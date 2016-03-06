#' Order Constrained Linear Optimization
#'
#' This function allows you to do silly things.
#' @param blah does blah things
#' @keywords oclo
#' @export
#' @examples
#' # Generate some data
#' n <- 25
#' x1 <- rnorm(n)
#' x2 <- nrorm(n)
#' y <- 10 + 5*x1 + 2*x2 + rnorm(n)
#' 
#' fit.oclo <- oclo(y~x1+x2)
#' 
oclo <- function(x, ...) {
  UseMethod("oclo", x)
}

#' 
#' @export
oclo.default <- function(x, data, ...) {
  stop("oclo() currently only accepts formulas")
}

#' 
#' @export
oclo.formula <- function(x, data=list(), ...) {
  # coerce to data.frame if a matrix
  if (!(class(data) %in% c("data.frame","data.table","matrix"))) { 
    stop("Must provide data of class of data.frame, data.table, or matrix")
  }
  if(class(data)=="matrix") {data <- data.frame(data)}

  # Remove intercept, if present
  cl <- match.call()  
  if(attr(terms(x),"intercept")==1) {
    x <- update(x, ~ . -1)
  }

  # create model.frame
  mframe <- model.frame(x, data)
  
  # create data
  X <- model.matrix(x,mframe)
  y <- as.matrix(mframe[,1, drop = FALSE])

  ocloData <- structure(list(y=y,
                             X=X,
                             model=mframe,
                             formula = cl), class="ocloData")  
  oclo(ocloData, ...)
}


bicFit <- function(tau,n,k,B) {
  knp <- sin(pi/2.0*tau*(n-k-1.0)/n)
  n * log(1.0 - knp^2) + (k-apply(B, 1, function(r) sum(r==0))) * log(n)  
}

crossover <- function(idx, ctpt, parent.A, parent.B, k) {
  c(parent.A[idx,1:ctpt],parent.B[idx,(ctpt+1):k])
}

### Genetic Algorithm Functions
# Initial betas
init.chromo <- function(k,
                        n.beta,
                        seed,
                        y,
                        X,
                        model.select) {

  if (!(is.null(X)) & !(class(X) %in% c("data.frame","data.table","matrix"))) { 
    warning(paste0("'seed' must be either 'random' or 'ols'. Defaulting to 'random'."))
  }

  n.K <- 2^(1:k)    
  n.K <- ceiling(n.K/sum(n.K)*n.beta)
  K <- sapply(1:k, function(kk) {
        c(rep(1,kk),rep(0,k-kk))
      })
  K <- K[,rep(1:ncol(K),n.K)]
  K <- apply(K,2,sample)
  lambda <- matrix(rnorm(ncol(X)*ncol(K)),nrow=ncol(K))
  if(seed=="ols") {
    lambda <- t(t(lambda) + coef(lm(y~X))[-1])
  }
  if(model.select) {lambda <- lambda*t(K)}
  return(lambda)
}

fitness <- function(y,X,lambda,n,k,surv.fun,repro.fun,n.elites) {
    ## survival fitness (tau)
    y.hat <- X%*%t(lambda)

    surv.fit <- abs(apply(y.hat,2,function(yhat) {surv.fun(y,yhat)}))

    ## reproductive fitness
    repro.fit <- repro.fun(surv.fit,n,k,lambda)

    # model size
    tracek <- apply(lambda,1,function(l)sum(l!=0))

    # best chromosomes 
    elites <- lambda[order(repro.fit,-surv.fit),][1:n.elites,]


    return(list(surv   = surv.fit,
                repro  = repro.fit,
                elites = elites,
                k      = tracek))
}

cull <- function(surv, n.beta, lambda) {
  # probability of surviving until reproduction
  p.survival <- surv/sum(surv)

  # Select mu - i.e. chromosomes that survive to reproduction
  mu.idx <- sample(1:length(surv),
                   ceiling(.75*n.beta),
                   replace=TRUE,
                   prob=p.survival)
  mu <- lambda[mu.idx,]
  rownames(mu) <- mu.idx
  return(mu)
}

initTrace <- function(trace.ga) {
  trace.ga$surv <- data.frame(mu.fit    = numeric(), 
                              sd.fit    = numeric(),
                              best.fit  = numeric())
  trace.ga$repro <- data.frame(mu.fit   = numeric(), 
                               sd.fit   = numeric(),
                               best.fit = numeric())
  trace.ga$k <- data.frame(k.mu = numeric(), 
                           k.sd = numeric())

  return(trace.ga)
}


addTrace <- function(trace.ga, fit) {
  trace.ga$k[nrow(trace.ga$k)+1,] <- c(mean(fit$k),
                                       sd(fit$k))
  trace.ga$surv[nrow(trace.ga$surv)+1,] <- c(mean(fit$surv),
                                             sd(fit$surv),
                                             max(fit$surv))
  trace.ga$repro[nrow(trace.ga$repro)+1,] <- c(mean(fit$repro),
                                               sd(fit$repro),
                                               min(fit$repro))
  return(trace.ga)
}

repro <- function(fit,mu.idx,n.beta,n.elites,k,mu,lambda,p.mutate,s.mutate,i) {
  p.reproduce <- rank(-fit$repro[mu.idx])
  p.reproduce <- p.reproduce/sum(p.reproduce)  
  lambda.idx <- sample(1:nrow(mu),n.beta-n.elites,replace=TRUE,prob=p.reproduce)

  # crossover
  crosspoint <- sample(1:(k-1),length(lambda.idx),replace=TRUE)
  parent.A <- mu[lambda.idx,]
  parent.B <- mu[sample(lambda.idx),]
  lambda <- do.call(rbind, 
                    Map(pryr::partial(crossover, k=k, 
                                parent.A=parent.A, 
                                parent.B=parent.B), 
                        1:length(lambda.idx), crosspoint))

  # mutate
  mutate.idx <- sample(1:prod(dim(lambda)),
                       rbinom(1,prod(dim(lambda)),p.mutate))
  if(i < 5) {
  lambda[mutate.idx] <- rnorm(length(mutate.idx),0,1)        
  } else if(i < 30) {
    lambda[mutate.idx] <- lambda[mutate.idx] + 
                          rnorm(length(mutate.idx),0,s.mutate)        
  } else {
          lambda[mutate.idx] <- lambda[mutate.idx] + 
                                rnorm(length(mutate.idx),0,s.mutate/2)
  } 
  
  # terminate unfit
  lambda <- lambda[!apply(lambda, 1, function(r) sum(r)==0),]
  
  # Add elites
  lambda <- rbind(fit$elites,lambda)

  return(lambda)
}

#' 
#' @export
oclo.ocloData <- function(gdata, ...,
                          oscale = TRUE, 
                          n.beta = 1000, 
                          n.gens=10,
                          model.select=TRUE,
                          surv.fun = pcaPP::cor.fk,
                          repro.fun = bicFit,
                          p.mutate = .01,
                          pdata = TRUE,
                          s.mutate = .25,
                          seed = "random",
                          n.elites = 5) {


  out <- structure(list(), class = "ocloFit")
  if(model.select) {
    class(out) <- c(class(out),"modelSelect")
  }
  trace.ga <- structure(list(), class="ocloTrace")

  y <- gdata$y
  X <- gdata$X
  k <- ncol(X)
  n <- length(y)

  # If not specified, default # of betas is 500 * # of predictors
  if(is.null(n.beta)) {
    n.beta <- 500*k
  }

  #### Genetic algorithm

  ## initialize chromosomes
  lambda <- init.chromo(k,n.beta,seed,y,X,model.select)

  # Initialize plotting variables
  if(pdata){
    trace.ga <- initTrace(trace.ga)
  }

  # Two-step GA for model selection
  if(model.select) {
    for(i in (1:n.gens)) {

      fit <- fitness(y,X,lambda,n,k,surv.fun,repro.fun,n.elites)

      # Plotting data
      if(pdata) {
        trace.ga <- addTrace(trace.ga, fit) 
      }

      # Survival Fitness
      mu <- cull(fit$surv, n.beta, lambda)

      # linear rank-selection by reproductive fitness
      lambda <- repro(fit,mu.idx=as.numeric(rownames(mu)),n.beta,n.elites,k,mu,lambda,p.mutate,s.mutate,i)

    }
    
    fit <- fitness(y,X,lambda,n,k,surv.fun,repro.fun,n.elites)
    if(pdata) { trace.ga <- addTrace(trace.ga, fit) }
    lambda <- lambda[order(fit$repro, -fit$surv),]


  } else {
    for(i in (1:n.gens)) {
      tracek <- apply(lambda,1,function(l)sum(l==0))

      ## survival fitness (tau)
      surv.fit <- abs(apply(X%*%t(lambda),2,function(yhat) {surv.fun(y,yhat)}))
      elites <- lambda[order(surv.fit),][1:5,]
  
      # Plotting data
      if(pdata) {
            out$tracefit$k <- rbind(out$tracefit$k,
                                    data.frame(k.mu=mean(k-tracek),
                                               k.sd=sd(k-tracek)))
            out$tracefit$surv <- rbind(out$tracefit$surv,
                                        data.frame(mu.fit=mean(surv.fit),best.fit=max(surv.fit)))
      }

      ##linear rank-selection by reproductive fitness
      p.reproduce <- rank(surv.fit)
      p.reproduce <- p.reproduce/sum(p.reproduce)  
      lambda.idx <- sample(1:nrow(lambda),n.beta-5,replace=TRUE,prob=p.reproduce)
  
      ## Perturb 
      # crossover

      crosspoint <- sample(1:(k-1),length(lambda.idx),replace=TRUE)
      parent.A <- lambda[lambda.idx,]
      parent.B <- lambda[sample(lambda.idx),]

      lambda <- do.call(rbind, Map(pryr::partial(crossover, 
                                           k=k, 
                                           parent.A=parent.A, 
                                           parent.B=parent.B), 
                                   1:length(lambda.idx), 
                                   crosspoint))
    
      # mutate
      mutate.idx <- sample(1:prod(dim(lambda)),
                           rbinom(1,prod(dim(lambda)),p.mutate))

# cow
# hard code minimum generations
# fix .25 / .125

      if(i < 10) {
        lambda[mutate.idx] <- lambda[mutate.idx] + rnorm(length(mutate.idx),0,.25)        
      } else {
              lambda[mutate.idx] <- lambda[mutate.idx] + rnorm(length(mutate.idx),0,.125)
      }
  
      # terminate unfit
      lambda <- lambda[!apply(lambda, 1, function(r) sum(r)==0),]
  
      # Add elites
      lambda <- rbind(elites,lambda)
    }
  
    surv.fit <- abs(apply(X%*%t(lambda),2,function(yhat) {surv.fun(y,yhat)}))
    repro.fit <- repro.fun(surv.fit,n,k,lambda)
    lambda <- lambda[order(repro.fit, surv.fit),]


  }

  # Unique solutions
  lambda <- unique(lambda)

  # rescale and add intercept to output (default)
  if (oscale) {
    scaling <- t(apply(lambda,1,function(beta) {
          coef(lm(y ~ X %*% beta))
        }))
    out$Beta <- cbind(`(Intercept)`=scaling[,1], lambda * scaling[,2])
    colnames(out$Beta)[-1] <- colnames(X)
    tau <- apply(cbind(1,X)%*%t(out$Beta),2,function(yhat) {pcaPP::cor.fk(y,yhat)})
    R.sq <- apply(cbind(1,X)%*%t(out$Beta),2,function(yhat) {cor(y,yhat)})^2
    bic <- bicFit(tau,n,k,lambda)
    out$Beta <- out$Beta[order(bic, -tau, -R.sq),]
    out$fits <- cbind(tau=tau[order(bic, -tau, -R.sq)],
                      R.sq=R.sq[order(bic, -tau, -R.sq)],
                      BIC=bic[order(bic, -tau, -R.sq)])
    out$coefficients <- out$Beta[1,]
    out$oscale <- data.frame(Intercept=scaling[,1], scale.fac=scaling[,2])
  } else {
  # add unit-scaled, non-oclo coefficients 
    out$Beta <- lambda
    out$Beta <- t(apply(out$Beta,1,function(b) b/sum(abs(b))))
    out$coefficients <- out$Beta[1,]
    names(out$coefficients) <- colnames(X)
  }

  out$call <- gdata$formula
  out$model <- gdata$model
  out$trace <- trace.ga

  out
}

#' 
#' @export
summary.ocloFit <- function(object, ...) {
  class(object) <- "summary.ocloFit"
  return(object)
}

#' 
#' @export
print.summary.ocloFit <- function(x, ...) {
  print.ocloFit(x, ...)
}

#' 
#' @export
print.ocloFit <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)

  # y.hat <- cbind(1,model.matrix(attr(x$model,"terms"),x$model))%*%coef(x)
  # r.sq <- 1 - sum((y.hat - x$model[,1])^2)/(sum((x$model[,1]-mean(x$model[,1]))^2))
  # tau <- cor.fk(y.hat,x$model[,1])
  # r <- cor(y.hat, x$model[,1])
  r.sq <- x$fits[,'R.sq'][1]
  r <- sqrt(x$fits[,'R.sq'][1])
  tau <- x$fits[,'tau'][1]
  bic <- x$fits[,'BIC'][1]

  cat("\nBIC\t\t\t:", round(bic,3),
      "\nPseudo-R^2\t\t:", round(r.sq,3),
      "\nPearson's r\t\t:", round(r,3),
      "\nKendall's tau-b\t\t:", round(tau,3),"\n")
}


#' 
#' @export
predict.ocloFit <- function(obj, X=NULL, ...) {
  if (!(is.null(X)) & !(class(X) %in% c("data.frame","data.table","matrix"))) { 
    stop("Data must be of class data.frame, data.table, or matrix")
  }

  if(length(X)==0) {    # Model matrix for fitted values
    mm <- cbind(1,model.matrix(attr(obj$model,"terms"),obj$model))
  } else {              # New model matrix from X
    mm <- cbind(1,model.matrix(attr(obj$model,"terms"),data.frame(X)))
  }

  fitted <- c(mm%*%coef(obj))
  names(fitted) <- 1:length(fitted)
  return(fitted)
}

#' 
#' @export
residuals.ocloFit <- function(obj, ...) {
  fitted <- predict(obj)
  return(obj$model[,1]-fitted)
}

#' 
#' @export
plot.ocloFit <- function(obj, ...) {
  if (length(obj$trace)==0) { 
    stop("Run oclo with pdata=TRUE to obtain plotting data.")
  }
  if(nrow(obj$trace$repro)==0) {
    pdata <- data.frame(obj$trace$surv)
    names(pdata) <- c("surv.mu","surv.best")
    par(mfrow=c(1,1))
    plot(x=1:nrow(pdata),y=pdata$surv.mu, type="l", col="red", ylim=c(0,1),
      xlab="Generation", ylab="Survival Fitness")
    lines(x=1:nrow(pdata),y=pdata$surv.best, col="blue")
  } else {
    pdata <- data.frame(cbind(obj$trace$surv,obj$trace$repro))
    names(pdata) <- c("surv.mu","surv.sd","surv.best","repro.mu","repro.sd","repro.best")
    par(mfrow=c(2,2))

    x <- 1:nrow(pdata)
    plot(x=1:nrow(pdata),y=pdata$surv.mu, type="l", col="red", ylim=c(0,1),
      xlab="Generation", ylab="Survival Fitness")
    polygon(c(x,rev(x)),c(pdata$surv.mu-pdata$surv.sd,rev(pdata$surv.mu+pdata$surv.sd)),col="lightgrey", border="darkgrey")
    lines(x=1:nrow(pdata),y=pdata$surv.mu, type="l", col="red")
    lines(x=1:nrow(pdata),y=pdata$surv.best, col="blue")

    plot(x=1:nrow(pdata),y=pdata$repro.mu, type="l", col="red",
      xlab="Generation", ylab="Reproductive Fitness")
    polygon(c(x,rev(x)),c(pdata$repro.mu-pdata$repro.sd,rev(pdata$repro.mu+pdata$repro.sd)),col="lightgrey", border="darkgrey")
    lines(x=1:nrow(pdata),y=pdata$repro.mu, type="l", col="red")    
    lines(x=1:nrow(pdata),y=pdata$repro.best, col="blue")

    pdata <- obj$trace$k
    x <- 1:nrow(pdata)
    plot(x=x,y=pdata$k.mu,type="l", ylim=c(0,ceiling(max(pdata$k.mu))), xlab="Generation", ylab="Average Model Size")
    polygon(c(x,rev(x)),c(pdata$k.mu-pdata$k.sd,rev(pdata$k.mu+pdata$k.sd)),col="lightgrey", border="darkgrey")
    lines(x=x,y=pdata$k.mu)
  }

  par(mfrow=c(1,1))
}


#' 
#' @export
jitterFit <- function(mod,            # Fitted oclo model
                      s     = .005,   # jitter sd
                      prob  = .25,    # probability of jitter
                      n     = 1e4,    # number of beta vectors
                      pdata = TRUE,   # Generate plotting data
                      ...) {

  out <- structure(list(), class="jitterFit")
  X <- model.matrix(attr(mod$model,"terms"),mod$model)
  y <- mod$model[,all.vars(attr(mod$model,"terms"))[1]]

  # Get all models with equivalent best fit
  best.idx <- which(mod$fits[,ncol(mod$fits)]==mod$fits[1,ncol(mod$fits)])
  oscale <- mod$oscale[best.idx,2]
  gbeta <- mod$Beta[best.idx,-1,drop=FALSE]/oscale

  # Create n new beta vectors
  m <- matrix(c(t(gbeta)),ncol=ncol(gbeta),nrow=n,byrow=TRUE)
  mut.idx <- which(m!=0)
  mut.idx <- sample(mut.idx,rbinom(1,length(mut.idx),prob))
  m[mut.idx] <- m[mut.idx] + rnorm(length(mut.idx),0,s)

  # Generate fitted values
  y.hat <- X%*%t(m)
  tau <- abs(apply(y.hat,2,function(yh) {pcaPP::cor.fk(y,yh)}))

  # Select unique fits that are better than the original
  idx1 <- which(tau>=mod$fits[1,'tau'])
  idx2 <- which(!duplicated(m[idx1,]))
  m <- m[idx1,][idx2,]
  tau <- tau[idx1][idx2]

  # oclo scaling
  scaling <- t(apply(m,1,function(beta) {
        coef(lm(y ~ X %*% beta))}))
  Beta <- cbind(`(Intercept)`=scaling[,1], m * scaling[,2])
  R.sq <- apply(y.hat[,idx1][,idx2],2,function(yh) {cor(y,yh)})^2
  ord <- order(-tau,-R.sq)

  # return best equivalent models
  fits <- cbind(tau,R.sq)
  fits <- fits[ord,]
  Beta <- Beta[ord,]
  best.idx <- which(fits[,2]==fits[1,2])

  out$fits <- fits[best.idx,]  
  out$beta <- Beta[best.idx,]
  out$Beta <- Beta

  return(out)
}

#' 
#' @export
summary.jitterFit <- function(object, ...) {
  class(object) <- "jitterFit.ocloFit"
  return(object)
}

#' 
#' @export
print.summary.jitterFit <- function(x, ...) {
  print.jitterFit(x, ...)
}

#' 
#' @export
print.jitterFit <- function(x, ...) {
  cat("\nCoefficients:\n")
  print(round(x$beta,5))
  cat("\nFit Metrics:\n")
  print(round(x$fits,5))

  # cat("\nBIC\t\t\t:", round(bic,3),
  #     "\nPseudo-R^2\t\t:", round(r.sq,3),
  #     "\nPearson's r\t\t:", round(r,3),
  #     "\nKendall's tau-b\t\t:", round(tau,3),"\n")
}

