

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
for(inc in 1:nc) {
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
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self")) {
   type <- match.arg(type)
   assign(".type", type, envir=baseenv())
   if(gcFirst) gc(FALSE)
   tic <- proc.time()[type]         
   assign(".tic", tic, envir=baseenv())
   invisible(tic)
}

toc <- function() {
   type <- get(".type", envir=baseenv())
   toc <- proc.time()[type]
   tic <- get(".tic", envir=baseenv())
   print(toc - tic)
   invisible(toc)
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
