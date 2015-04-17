source("C:\\Research\\R\\cdfmvna.r")


data <- "C:/Research/R/data";
setwd(data);
dt<- read.csv("CS_uncorrelated.csv",header = TRUE);

dt = dt[1:5000,];
print(dim(dt));


nobs = 5000; 		#	Number of observations in data set
nvar = 5;			  #	Number of variables
nc = 5;			    #	Number of alternatives


# Utility Specification below

atrAlt1 <- cbind(dt$x01,dt$x02,dt$x03,dt$x04,dt$x05);
atrAlt2 <- cbind(dt$x06,dt$x07,dt$x08,dt$x09,dt$x10);
atrAlt3 <- cbind(dt$x11,dt$x12,dt$x13,dt$x14,dt$x15);
atrAlt4 <- cbind(dt$x16,dt$x17,dt$x18,dt$x19,dt$x20);
atrAlt5 <- cbind(dt$x21,dt$x22,dt$x23,dt$x24,dt$x25);


dtaX <- cbind(atrAlt1,atrAlt2,atrAlt3,atrAlt4,atrAlt5);
dtaCh <- cbind(dt$f1,dt$f2,dt$f3,dt$f4,dt$f5);            # The dummy variable indicating twhather the alternative is the chosen
dtaChNum = (apply(dtaCh,1,function(x) which(x==1)));            



appdel = kronecker(matrix(1,1,nc-1),dtaCh);
appdel_tplt =  matrix(nrow=0,ncol=nc*(nc-1));
for(inc in 1:nc)
{
  appdel_tplt_row = matrix(nrow=1,ncol=0);
  for(jnc in 1:nc)
  {
    if(jnc != inc)
    {
      temp = matrix(0,nrow=1,ncol=nc);
      temp[1,jnc] = 1;
      appdel_tplt_row = cbind(appdel_tplt_row,temp);
    }
  }
  appdel_tplt = rbind(appdel_tplt,appdel_tplt_row);
}
appdel1 = appdel_tplt[dtaChNum,];
appdel = t(appdel);
appdel1 = t(appdel1);
appdiff =  appdel -  appdel1;
  
ID = 0.5*diag(nc-1) + 0.5*matrix(1,nrow=nc-1,ncol=nc-1);

outProd = 0;  finalIter <<- 0;
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
   type <- match.arg(type)
   assign(".type", type, envir=baseenv())
   if(gcFirst) gc(FALSE)
   tic <- proc.time()[type]         
   assign(".tic", tic, envir=baseenv())
   invisible(tic)
}

toc <- function()
{
   type <- get(".type", envir=baseenv())
   toc <- proc.time()[type]
   tic <- get(".tic", envir=baseenv())
   print(toc - tic)
   invisible(toc)
}
lpr <-function(x)
{
  x_b = x[1:nvar];
  x_sig = x[(nvar+1):(2*nvar)];
  v2 = as.vector(kronecker(matrix(1,nrow=nc,ncol=1),x_b));
  v2 = v2 * t(dtaX);
  
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
  
  ll = matrix(nrow=nobs,ncol=1);

  for(iobs in (1:nobs))
  {
    seed20=seednext;
    xx = matrix(rep(dtaX[iobs,],nc-1),nrow=nc*(nc-1),ncol=nvar,byrow = TRUE);
    xx = (-appdiff[,iobs]) * xx 

    xx1 = matrix(nrow=(nc-1),ncol=nvar);
    for(inc in 1:(nc-1))
    {
      xx1[inc,] = colSums(xx[((inc-1)*nc+1):(inc*nc),]);
    }
    omega = xx1%*%omega11%*%t(xx1) + ID;
    sqrt_om = sqrt(diag(omega));
    kk1 = dutil[,iobs]/sqrt_om;
    kk2 = omega;kk2=t(kk2/sqrt_om)/sqrt_om;
    diag(kk2) <- rep(1,nc-1);
    retval = cdfmvna(kk1,kk2,seed20);
    seednext = retval[2];
    ll[iobs,1]= retval[1];
  } 
  ll = log(ll);
  ll = mean(ll);  
  
  return(ll);
}



lgd <-function(x)
{
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
  for(iobs in (1:nobs))
  {
    seed20=seednext;
    xx = matrix(rep(dtaX[iobs,],nc-1),nrow=nc*(nc-1),ncol=nvar,byrow = TRUE);
    xx = (-appdiff[,iobs]) * xx 
    
    xx1 = matrix(nrow=(nc-1),ncol=nvar);
    for(inc in 1:(nc-1))
    {
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


bInit = c(0,0,0,0,0,1,1,1,1,1);


seed10 = 10000;


finalIter <<- 0; tic(); 
res1 <- optim(bInit,lpr,lgd,method="BFGS",control = list(maxit=100, trace=3, fnscale=-1, REPORT=1), hessian = TRUE);
print(res1);  toc();

finalIter = 1;  tic();
res2 <- optim(res1$par,lpr,lgd,method="BFGS",control = list(maxit=1, trace=3, fnscale=-1, REPORT=1));
print(res2);   toc();

cov2 =  solve(res1$hessian) %*% outProd %*% solve(res1$hessian)/nobs/nobs  ;

estimate = res1$par;
stdErr = (sqrt(diag(cov2)));

print(cbind(estimate,stdErr));