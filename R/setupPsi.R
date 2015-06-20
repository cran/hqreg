# Determine vector psi used to compute lambda.max
setupPsi <- function(yy, X, loss, gamma, c, penalty.factor)
{
  ind <- which(penalty.factor!=0)
  if(length(ind) == ncol(X)) {
    if(loss == "huber") {
      psi <- ifelse(abs(yy)>gamma, sign(yy), yy/gamma)
    } else if(loss == "quantile") {
      psi <- ifelse(abs(yy)>gamma, sign(yy), yy/gamma)+c
    } else {
      psi <- yy
    }
  } else {
    if(loss == "huber") {
      fit <- lm(yy~X[,-ind]+0)
      r <- fit$residuals
      psi <- ifelse(abs(r)>gamma, sign(r), r/gamma)
    } else if(loss == "quantile") {
      fit <- lm(yy~X[,-ind]+0)
      r <- fit$residuals
      psi <- ifelse(abs(r)>gamma, sign(r), r/gamma)+c
    } else {
      fit <- lm(yy~X[,-ind]+0)
      psi <- fit$residuals
    } 
  }
  psi
}

