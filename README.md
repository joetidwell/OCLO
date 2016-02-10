gemmR
=====

The General Monotone Model (GeMM) in `R` with Order Constrained Linear Optimization (OCLO). This is an R version of the non-parametric regression algorithm proposed by Dougherty & Thomas (2012). The base OCLO algorithm works as follows:

1. Generate random weights with a genetic algorithm
2. Test the fit between the sum of the weighted predictors and an outcome using Bayesian Information Criterion based on tau-to-r transformation.
3. Select top ___ of candidate beta vectors, mutate and generate new weights, repeat for some number of generations.
4. Repeat entire process some number of times to check for convergence.

The resulting `ocloFit` object is an S3 class object. We're working on increasing usability and writing a proper vignette. Information on the basis for GeMM can be found in the [original paper](http://www.bsos.umd.edu/psyc/dougherty/pdf%20articles/DoughertyThomas2012Rev.pdf).

Installation
-----

You can install `oclo` via GitHub:

```r
library(devtools)
install_github("joetidwell/oclo")
```

Licensing
-----

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.