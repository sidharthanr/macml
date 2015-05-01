#' lpr
#' 
#' @param x
#' 
#' @export
lpr <- function(x) {
  n_var = length(x) / 2
  beta = x[1:n_var]
  sigma = x[(n_var + 1):(2 * n_var)]
  v2 = as.vector(kronecker(matrix(1, nrow = n_alts, ncol = 1), beta))
  v2 = v2 * t(attribute_mtx)
  
  ut = matrix(nrow = n_alts, ncol = n_obs)
  for(i in 1:n_alts){
    ut[i,] = colSums(v2[((i - 1) * n_var + 1):(i * n_var), ])
  }
  
  dd = appdiff * (kronecker(matrix(1, nrow = n_alts - 1, ncol = 1), ut))
  dutil = matrix(nrow = n_alts - 1, ncol = n_obs)
  for(i in 1:(n_alts - 1)){
    dutil[i, ] =  colSums(dd[((i - 1) * n_alts + 1):(i * n_alts), ])
  }
  
  
  omega11 = diag(sigma) ^ 2
  
  ll = matrix(nrow = n_obs, ncol = 1)
  
  # TODO: This will be ***lots*** faster if we can vectorize it.
  for(i in 1:n_obs){
    xx = matrix(rep(attribute_mtx[i, ], n_alts - 1),
                nrow= n_alts * (n_alts - 1),
                ncol= n_var, byrow = TRUE)
    xx = (-appdiff[, i]) * xx 
    
    xx1 = matrix(nrow = n_alts - 1, ncol = n_var)
    for(j in 1:n_alts - 1) {
      xx1[j, ] = colSums(xx[((j - 1) * n_alts + 1):(j * n_alts), ])
    }
    
    omega = xx1 %*% omega11 %*% t(xx1) + ID
    sqrt_om = sqrt(diag(omega))
    
    kk1 = dutil[, i] / sqrt_om
    kk2 = t(omega / sqrt_om)/sqrt_om
    seed20 = seed10 
    retval = cdfmvna(kk1,kk2, seed20)
    seednext = retval[2]
    ll[i, 1] = retval[1]
  } 
  
  
  mean(log(ll))
  
}

#' Gradient of the CML
#'
#' \code{l_gradient} returns the gradient of the composite marginal likelihood
#' function.
#' 
#' @param x
#' 
#' @export
lgd <-function(x) {
  print(x);
  x_b = x[1:nvar];
  x_sig = x[(nvar+1):(2*nvar)];
  v2 = as.vector(kronecker(matrix(1,nrow=nc,ncol=1),x_b));
  v2 = v2 * t(dtaX);
  
  lowerSelMat = lower.tri(matrix(nrow=nc-1,ncol=nc-1),diag=F);
  
  
  ut = matrix(nrow=nc,ncol=nobs);
  for(inc in 1:nc){
    ut[inc,] = colSums(v2[((inc-1)*nvar+1):(inc*nvar),]);
  }
  dd = appdiff * (kronecker(matrix(1,nrow=nc-1,ncol=1),ut));
  dutil = matrix(nrow=(nc-1),ncol=nobs);
  for(inc in 1:(nc-1)){
    dutil[inc,] =  colSums(dd[((inc-1)*nc+1):(inc*nc),]);
  }
  
  seednext=seed10;
  omega11 = matrix(0,nrow=nvar,ncol=nvar);
  diag(omega11) <-  x_sig;
  omega11 = omega11^2;
  gg = matrix(nrow=nobs,ncol=length(x));
  for(iobs in (1:nobs)) {
    seed20=seednext;
    xx = matrix(rep(dtaX[iobs,],nc-1),nrow=nc*(nc-1),ncol=nvar,byrow = TRUE);
    xx = (-appdiff[,iobs]) * xx 
    
    xx1 = matrix(nrow=(nc-1),ncol=nvar);
    for(inc in 1:(nc-1)) {
      xx1[inc,] = colSums(xx[((inc-1)*nc+1):(inc*nc),]);
    }
    omega = xx1%*%omega11%*%t(xx1) + ID;
    om = diag(omega);
    sqrt_om = sqrt(om);
    kk1 = dutil[,iobs]/sqrt_om;
    kk2 = omega;kk2=t(kk2/sqrt_om)/sqrt_om;
    diag(kk2) <- rep(1,nc-1);
    
    retval = pdfmvna(kk1,kk2,seed20); 
    
    seednext = retval[length(retval)];
    zz = retval[1];
    ggg = matrix(retval[2:(length(retval)-1)],nrow=1); 
    
    
    kk2ext = kk2[lowerSelMat];
    
    # Start getting gradients with respect to parameters #
    
    ga = 1/(sqrt(om));
    
    ggw = -ggg[1:nc-1]%*%(ga*xx1);
    gr1 = ggg[nc:ncol(ggg)];
    gr2 = matrix(xx1,ncol=1);
    gr3 = gr2%*%t(gr2);
    grn = 0.5*diag(nvar);
    grn = kronecker(grn,matrix(1,nrow=nc-1,ncol=nc-1));
    grn = grn+(grn==0);
    
    gr3 = gr3*grn;
    cc1 = ncol(gr1);
    cc2 = (nvar*(nvar+1))/2;
    oms = combinations(nc-1,2);
    oms5 = om[oms[,1]];
    oms6 = om[oms[,2]];
    oms1 = oms5*oms6;
    
    
    k=1;j=1;grf = matrix(nrow=1,ncol=0);
    while(k <= nvar*(nc-1))
    {
      l = k;
      while(l <= (nc-1)*j)
      {
        ex = gr3[k:(k+nc-2),l:(l+nc-2)]; 
        ex = 2*ex;
        ex1 = diag(ex);  
        ex2 = ex[lowerSelMat];      
        ex3 = (((ex1[oms[,1]])/oms5)+((ex1[oms[,2]])/oms6))*(oms1);  
        grr1 = (ex2*sqrt(oms1)-0.5*kk2ext*(ex3))/oms1;
        grr2 = -0.5*((kk1)/(om))*(ex1);
        grr1 = ggg%*%(c(grr2,grr1));
        grf = cbind(grf,(2*grr1%*%x_sig[j])); 
        l=l+(nc-1);
      }    
      j=j+1;
      k=k+nc-1;  
    }
    gg[iobs,] = cbind(ggw,grf)/zz;        
  }
  if(finalIter==1){
    outProd <<- t(gg) %*% gg;
  }
  gg = apply(gg,2,mean); 
  
  return(gg);
}
