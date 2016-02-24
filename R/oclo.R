library(pcaPP)
library(pryr)

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


#' 
#' @export
oclo.ocloData <- function(gdata, ...,
                          oscale = TRUE, 
                          n.beta = 1000, 
                          n.gens=10,
                          model.select=TRUE,
                          surv.fun = cor.fk,
                          repro.fun = bicFit,
                          p.mutate = .01,
                          pdata = FALSE,
                          s.mutate = .25,
                          seed = "random") {

  out <- structure(list(), class = "ocloFit")

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
    out$tracefit <- list()
    out$tracefit$surv <- data.frame(mu.fit=numeric(), best.fit=numeric())
    out$tracefit$repro <- data.frame(mu.fit=numeric(), best.fit=numeric())
    out$tracefit$k <- data.frame(k.mu=numeric(),k.sd=numeric())
  }

  # Two-step GA for model selection
  if(model.select) {
    for(i in (1:n.gens)) {
      tracek <- apply(lambda,1,function(l)sum(l==0))

      ## survival fitness (tau)
      surv.fit <- abs(apply(X%*%t(lambda),2,function(yhat) {surv.fun(y,yhat)}))

      repro.fit <- repro.fun(surv.fit,n,k,lambda)
      elites <- lambda[order(repro.fit,surv.fit),][1:5,]

      p.survival <- surv.fit/sum(surv.fit)
  
      # Select mu
      mu.idx <- sample(1:length(surv.fit),
                       ceiling(.75*n.beta),
                       replace=TRUE,
                       prob=p.survival)
      mu <- lambda[mu.idx,]
  
      ### reproductive fitness 
      # Plotting data
      if(pdata) {
            out$tracefit$k <- rbind(out$tracefit$k,
                                    data.frame(k.mu=mean(k-tracek),
                                               k.sd=sd(k-tracek)))
            out$tracefit$surv <- rbind(out$tracefit$surv,
                                        data.frame(mu.fit=mean(surv.fit),best.fit=max(surv.fit)))
            out$tracefit$repro <- rbind(out$tracefit$repro,
                                        data.frame(mu.fit=mean(repro.fit),best.fit=min(repro.fit)))
      }

      ##linear rank-selection by reproductive fitness
      p.reproduce <- rank(-repro.fit[mu.idx])
      p.reproduce <- p.reproduce/sum(p.reproduce)  
        lambda.idx <- sample(1:nrow(mu),n.beta-5,replace=TRUE,prob=p.reproduce)
  
      ## Perturb 
      # crossover

      crosspoint <- sample(1:(k-1),length(lambda.idx),replace=TRUE)
      parent.A <- mu[lambda.idx,]
      parent.B <- mu[sample(lambda.idx),]

      lambda <- do.call(rbind, Map(partial(crossover, k=k, parent.A=parent.A, parent.B=parent.B), 1:length(lambda.idx), crosspoint))
    
      # mutate
      mutate.idx <- sample(1:prod(dim(lambda)),rbinom(1,prod(dim(lambda)),p.mutate))
      if(i < 10) {
        lambda[mutate.idx] <- lambda[mutate.idx] + rnorm(length(mutate.idx),0,s.mutate)        
      } else if(i < 15) {
              lambda[mutate.idx] <- lambda[mutate.idx] + rnorm(length(mutate.idx),0,s.mutate/2)
      } else {
              lambda[mutate.idx] <- lambda[mutate.idx] + rnorm(length(mutate.idx),0,.001)
      }
  
      # terminate unfit
      lambda <- lambda[!apply(lambda, 1, function(r) sum(r)==0),]
  
      # Add elites
      lambda <- rbind(elites,lambda)
    }
    
    surv.fit <- abs(apply(X%*%t(lambda),2,function(yhat) {surv.fun(y,yhat)}))
    repro.fit <- repro.fun(surv.fit,n,k,lambda)
    lambda <- lambda[order(repro.fit, surv.fit),]


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

      lambda <- do.call(rbind, Map(partial(crossover, 
                                           k=k, 
                                           parent.A=parent.A, 
                                           parent.B=parent.B), 
                                   1:length(lambda.idx), 
                                   crosspoint))
    
      # mutate
      mutate.idx <- sample(1:prod(dim(lambda)),
                           rbinom(1,prod(dim(lambda)),p.mutate))

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
    tau <- apply(cbind(1,X)%*%t(out$Beta),2,function(yhat) {cor.fk(y,yhat)})
    R.sq <- apply(cbind(1,X)%*%t(out$Beta),2,function(yhat) {cor(y,yhat)})^2
    bic <- bicFit(tau,n,k,lambda)
    out$Beta <- out$Beta[order(bic, -tau, -R.sq),]
    out$fits <- cbind(tau=tau[order(bic, -tau, -R.sq)],
                      R.sq=R.sq[order(bic, -tau, -R.sq)],
                      BIC=bic[order(bic, -tau, -R.sq)])
    out$coefficients <- out$Beta[1,]
  } else {
  # add unit-scaled, non-oclo coefficients 
    out$Beta <- lambda
    out$Beta <- t(apply(out$Beta,1,function(b) b/sum(abs(b))))
    out$coefficients <- out$Beta[1,]
    names(out$coefficients) <- colnames(X)
  }

  out$call <- gdata$formula
  out$model <- gdata$model

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
  pdata <- obj$trace

  if (is.null(pdata)) { 
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
    names(pdata) <- c("surv.mu","surv.best","repro.mu","repro.best")
    par(mfrow=c(2,2))
    plot(x=1:nrow(pdata),y=pdata$surv.mu, type="l", col="red", ylim=c(0,1),
      xlab="Generation", ylab="Survival Fitness")
    lines(x=1:nrow(pdata),y=pdata$surv.best, col="blue")
    plot(x=1:nrow(pdata),y=pdata$repro.mu, type="l", col="red",
      xlab="Generation", ylab="Reproductive Fitness")
    lines(x=1:nrow(pdata),y=pdata$repro.best, col="blue")

    pdata <- obj$tracefit$k
    x <- 1:nrow(pdata)
    plot(x=x,y=pdata$k.mu,type="l", ylim=c(0,ceiling(max(pdata$k.mu))), xlab="Generation", ylab="Average Model Size")
    polygon(c(x,rev(x)),c(pdata$k.mu-pdata$k.sd,rev(pdata$k.mu+pdata$k.sd)),col="lightgrey", border="darkgrey")
    lines(x=x,y=pdata$k.mu)
  }

  par(mfrow=c(1,1))
}
