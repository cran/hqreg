\name{hqreg-package}
\alias{hqreg-package}
\docType{package}
\title{Regularization Paths for Lasso or Elastic-net Penalized Huber Loss Regression and Quantile Regression}
\description{Efficient algorithms for fitting regularization paths for lasso or elastic-net penalized regression models with 
  Huber loss, quantile loss or squared loss.}
\details{
\tabular{ll}{
Package: \tab hqreg\cr
Type: \tab Package\cr
Version: \tab 1.4-1\cr
Date: \tab 2024-09-23\cr
License: \tab GPL-3\cr
}
Very simple to use. Accepts \code{X,y} data for regression models, and
produces the regularization path over a grid of values for the tuning
parameter \code{lambda}. Also provides functions for plotting, prediction and parallelized cross-validation.
}
\references{Yi, C. and Huang, J. (2017) 
  \emph{Semismooth Newton Coordinate Descent Algorithm for 
  Elastic-Net Penalized Huber Loss Regression and Quantile Regression},
  \doi{10.1080/10618600.2016.1256816} \cr
  \emph{Journal of Computational and Graphical Statistics} \cr
}
\author{Congrui Yi <eric.ycr@gmail.com>}
\examples{
X = matrix(rnorm(1000*100), 1000, 100)
beta = rnorm(10)
eps = 4*rnorm(1000)
y = drop(X[,1:10] \%*\% beta + eps) 

# Huber loss
fit1 = hqreg(X, y)
coef(fit1, 0.01)
predict(fit1, X[1:5,], lambda = c(0.02, 0.01))
cv.fit1 = cv.hqreg(X, y)
plot(cv.fit1)

# Quantile loss
fit2 = hqreg(X, y, method = "quantile", tau = 0.2)
plot(fit2)

# Squared loss
fit3 = hqreg(X, y, method = "ls", preprocess = "rescale")
plot(fit3, xvar = "norm")
}
\keyword{models}
\keyword{regression}
