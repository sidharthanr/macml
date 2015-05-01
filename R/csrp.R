#' Cross sectional random parameters by MACML
#' 
#' \code{csrp} estimates a cross-sectional random parameters discrete choice 
#' model applying the MACML estimator.
#' 
#' @param dt a data frame with attribute measurements, etc.
#' @param n_obs number of observations. If \code{NULL}, defaults to 
#'   \code{nrow(dt)}
#' @param n_var number of attribute variables.
#' @param n_alts number of alternatives available to agents
#' @param initial_values intial parameters for the likelihood function. If 
#'   \code{NULL}, assume mean of zero and standard deviation of 1 for each 
#'   attribute.
#' @param seed10 random seed. If \code{NULL}, use 10,000 (value from original code)
#'
#' @export
csrp <- function(dt, n_obs, n_var, n_alts, initial_values, seed10){
  
  #check inputs ----
  if(is.na(n_obs)){
    n_obs <- nrow(dt)
  }
  
  if(is.na(initial_values)){
    initial_values <- c(rep(0, n_var), rep(1, n_var))
  }
  
  if(is.na(seed10){
    seed10 <- 10000
  }
  
  # manipulate dt into matrices ----
  attribute_mtx <- as.matrix(dt %>% select(x01:x25))
  choice_mtx <- as.matrix(dt %>% select(f1:f5))
  choice_index <- apply(dtaCh, 1, function(x) which(x == 1))
  
  appdel_tplt <- matrix(nrow = 0, ncol = n_alts * (n_alts - 1))
  for(i in 1:n_alts) {
    appdel_tplt_row = matrix(nrow = 1, ncol = 0)
    for(j in 1:n_alts) {
      if(j != i) {
        temp <- matrix(0, nrow = 1, ncol = n_alts)
        temp[1, j] <- 1
        appdel_tplt_row <- cbind(appdel_tplt_row, temp)
      }
    }
    appdel_tplt = rbind(appdel_tplt, appdel_tplt_row)
  }
  
  appdel <- t(kronecker(matrix(1, 1, n_alts - 1), choice_mtx))
  appdel1 <- t(appdel_tplt[choice_index, ])
  
  appdiff <- appdel -  appdel1 
  
  # Estimate first parameters ----
  ID <- 0.5 * diag(n_alts - 1) +  
    0.5 * matrix(1, nrow = n_alts - 1, ncol = n_alts - 1)
  outProd <- 0
  finalIter <- 0
  
  
  res1 <- optim(
    initial_values, 
    fn = lpr, gr = lgd, 
    method="BFGS",  hessian = TRUE,
    control = list(maxit = 1, trace = 3, fnscale = -1, REPORT = 1)
  )
  
}