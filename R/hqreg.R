hqreg <- function (X, y, loss = c("huber", "quantile", "squared"), gamma=quantile(abs(y-mean(y)), 0.05), tau = 0.5, alpha=1, nlambda=100, lambda.min = ifelse(nrow(X)>ncol(X), 0.001, 0.05), lambda, 
		  preprocess = c("standardize", "rescale", "none"),  screen = c("ASR", "SR", "none"), max.iter = 10000, eps = 1e-7, dfmax = ncol(X)+1, penalty.factor=rep(1, ncol(X)), message = FALSE) {
  
  # Error checking
  loss <- match.arg(loss)
  preprocess <- match.arg(preprocess)
  screen <- match.arg(screen)
  if (missing(lambda) && nlambda < 2) stop("nlambda should be at least 2")
  if (alpha < 0 || alpha > 1) stop("alpha should be between 0 and 1")
  if (alpha == 0) stop("alpha = 0 is not supported")
  if (loss == "huber" && gamma <= 0) stop("gamma should be positive for Huber loss")
  if (loss == "quantile" && (tau < 0 || tau > 1)) stop("tau should be between 0 and 1 for quantile loss")
  if (length(penalty.factor)!=ncol(X)) stop("the length of penalty.factor should equal the number of columns of X")

  call <- match.call()
  # Include a column for intercept
  n <- nrow(X)
  XX <- cbind(rep(1,n), X)
  penalty.factor <- c(0, penalty.factor) # no penalty for intercept term
  p <- ncol(XX)

  ym <- 0
  yq <- 0
  if (loss == "huber" || loss == "squared") {
    ym <- mean(y)
    yy <- y - ym
  }
  if (loss == "quantile") {
    yq <- quantile(y, tau)
    yy <- y - yq
    # initialize gamma to determine lambda sequence
    gamma <- quantile(abs(yy), 0.05)
    if (gamma < 0.00001) gamma = 0.00001
  }

  # Constant for quantile loss
  c <- 2*tau-1

  # For Huber loss, scale gamma accordingly when the response is standardized


  # Set up Psi if lambda = NULL
  psi <- 0
  user <- 0
  if (missing(lambda)) {
    lambda <- double(nlambda)
    psi <- setupPsi(yy, X, loss, gamma, c, penalty.factor)
  } else {
    nlambda <- length(lambda)
    user <- 1
  }
  
  # Flag for preprocessing and screening
  ppflag = switch(preprocess, standardize = 1, rescale = 2, none = 0)
  scrflag = switch(screen, ASR = 1, SR = 2, none = 0)
  # Fitting
  if (loss == "huber") {
    fit <- .C("huber", double(p*nlambda), integer(nlambda), as.double(lambda), integer(1), integer(1), as.double(XX), as.double(yy), as.double(psi), as.double(penalty.factor), 
              as.double(gamma), as.double(alpha), as.double(eps), as.double(lambda.min), as.integer(nlambda), as.integer(n), as.integer(p), as.integer(ppflag),
              as.integer(scrflag), as.integer(dfmax), as.integer(max.iter), as.integer(user), as.integer(message))
  } else if (loss == "quantile") {
    fit <- .C("quant", double(p*nlambda), integer(nlambda), as.double(lambda), integer(1), integer(1), as.double(XX), as.double(yy), as.double(psi), as.double(penalty.factor), 
              as.double(gamma), as.double(tau), as.double(alpha), as.double(eps), as.double(lambda.min), as.integer(nlambda), as.integer(n), as.integer(p), 
              as.integer(ppflag), as.integer(scrflag), as.integer(dfmax), as.integer(max.iter), as.integer(user), as.integer(message))
  } else {
    fit <- .C("squared", double(p*nlambda), integer(nlambda), as.double(lambda), integer(1), integer(1), as.double(XX), as.double(yy), as.double(psi), as.double(penalty.factor), 
              as.double(alpha), as.double(eps), as.double(lambda.min), as.integer(nlambda), as.integer(n), as.integer(p), as.integer(ppflag), as.integer(scrflag),
              as.integer(dfmax), as.integer(max.iter), as.integer(user), as.integer(message))
  }

  beta <- matrix(fit[[1]],nrow = p)
  iter <- fit[[2]]
  lambda <- fit[[3]]
  saturated <- fit[[4]]
  nv <- fit[[5]]
  # Eliminate saturated lambda values
  ind <- !is.na(iter)
  beta <- beta[, ind]
  iter <- iter[ind]
  lambda <- lambda[ind]
  
  # Intercept
  if (loss == "huber" || loss == "squared") beta[1,] <- beta[1,] + ym
  if (loss == "quantile") beta[1,] <- beta[1,] + yq
  
  # Names
  vnames <- colnames(X)
  if (is.null(vnames)) vnames=paste0("V",seq(p-1))
  vnames <- c("(Intercept)", vnames)
  dimnames(beta) <- list(vnames, round(lambda, 4))

  # Output
  structure(list(call = call,
                 beta = beta,
                 iter = iter,
                 saturated = saturated,
                 lambda = lambda,
                 alpha = alpha,
                 gamma = if (loss == "huber") gamma else NULL,
		 tau = if (loss == "quantile") tau else NULL,
                 penalty.factor = penalty.factor[-1],
		 loss = loss,
                 nv = nv),
            class = "hqreg")
}
