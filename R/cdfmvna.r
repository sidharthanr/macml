##library(combinat)
library(gtools)
library(pbivnorm)

sys_randper = 1;

	m = length(a);fillerMatrix =lower.tri(diag(m),diag=F);  
	rRowCol <- which(fillerMatrix==T, arr.in=TRUE); rRow = rRowCol[,1]; rCol = rRowCol[,2];
#' The probability density function of the mvna
#' 
#' @param a 
#' @param r 
#' @param s
#' 
#' @export
#' 
pdfmvna <- function(a,r,s) {

	w = 0; s1 = 0;	
	if(sys_randper ==0){
		w <-  getComb(m,2);
		s1 = s;
	}
	else{
		w = matrix(nrow=0,ncol=m);
		set.seed(s);
		for(z in 1:sys_randper)
		{	
			aa = cbind(seq(1:m),runif(m));
			aa = aa[order(aa[,2])];
			w <- rbind(w,aa);
		}
		s1 = ceiling(runif(1)*1000000000 );
	}
	
	#w = matrix(1:m,nrow=1); ## TEMP OVERRIDE FOR SAME MACML PERMUTATION
	w = t(w); 
	n1 = ncol(w);
	p = 0;
	j = 1;

	a <- 	a * (a < 5.7) + 5.7*(a >= 5.7);
	ab = t(a);
	condpass=matrix(nrow=1,ncol=0);
  for(j in 1:n1)
	{
		x = ab[w[,j]];
		xvec = rep(x,m);
		rho = as.vector(r[w[,j],w[,j]]);		
		y = x[sapply(1:m,rep,m)];
		z = pbivnorm(xvec,y,rho);

		z = matrix(z,nrow=m,ncol=m);
		z3 = pnorm(x);

		z = z - diag(diag(z)) + diag(z3);

		z1 = matrix(pnorm(xvec) * pnorm(y),nrow=m,ncol=m);
		z2 = 1-z3;
		cond = 1;

		for(k in 3:m)
		{
			omega21 = z[k,1:k-1]-z1[k,1:k-1];
			omega11 = z[1:k-1,1:k-1]-z1[1:k-1,1:k-1];
			cm = chol2inv(chol(omega11));#chol2inv(chol/solve
			condk = z3[k] + omega21%*%(cm)%*%z2[1:k-1] ;   
			cond=cond*condk;
		}
		pcomb = z[1,2]*cond;
		condpass = cbind(condpass,cond);
		p = p+pcomb;  
	}
	
	gwfinal = matrix(0,nrow=1,ncol=m);
	grfinal = matrix(0,nrow=1,ncol=((m-1)*(m)/2));
	for(j in 1:n1)
	{
		x = ab[w[,j]];
		xvec = rep(x,m);
		rho = as.vector(r[w[,j],w[,j]]);		
		y = x[sapply(1:m,rep,m)];
		z = pbivnorm(xvec,y,rho);

		z = matrix(z,nrow=m,ncol=m);
		z3 = pnorm(x);     

		z = z - diag(diag(z)) + diag(z3);

		z1 = matrix(pnorm(xvec) * pnorm(y),nrow=m,ncol=m);
		z2 = 1-z3;
		
		y = matrix(y,nrow=m,ncol=m);
		rho = matrix(rho,nrow=m,ncol=m);
    		 
    rho1 <- rho;
    diag(rho1) <- rep(0,m);
    rho2 = sqrt(1-rho1^2);
    
    g3 = dnorm(x);
    g5 = (y-rho1*x)/rho2;
    g10 = g3*pnorm(g5);
    g11 = g3*pnorm(y);
    g12 = g10-g11;
    g13 = g12;   
    g14 = g3*(1-2*z3);
    g15 = g13;diag(g15) = g14;
    g20 = -g3;
  
    g25 = (1/rho2)*g3*dnorm(g5);
    diag(g25) <- rep(0,m);    
    g25Vec =  g25[fillerMatrix];
    
    rhovec = matrix(nrow=1,ncol=0);   
    for(idim in (1:(ncol(rho)-1)))
    {
      rhovec = cbind(rhovec,matrix(rho[idim,((idim+1):ncol(rho))],nrow=1));    
    }
    
    g30 = cbind(g10[1,2],g10[2,1],matrix(0,nrow=1,ncol=m-2));
    g53 = cbind(g25[1,2],matrix(0,nrow=1,ncol=(ncol(rhovec)-1)));
  
    gw1 = matrix(0,nrow=1,ncol=m);
    gr1 = matrix(0,nrow=1,ncol=ncol(rhovec));
    for(k in (3:m))
    {
      omega21 = z[k,1:k-1]-z1[k,1:k-1];
      omega11 = z[1:k-1,1:k-1]-z1[1:k-1,1:k-1];
      invomg11 = chol2inv(chol(omega11));
			condk = z3[k] + omega21%*%(invomg11)%*%z2[1:k-1] ;   
			g31 = matrix(0,nrow=1,ncol=m);
			g31[k] = g3[k];
      g40=matrix(nrow=1,ncol=0);
      g46=matrix(nrow=1,ncol=0);
      g51=matrix(nrow=1,ncol=0);
      g81 = z2[1:k-1];
      for(l in (1:(k-1)))
      {
        g35 = matrix(0,nrow=k-1,ncol=k-1); 
        g35[l,(1:(k-1))] = g15[l,(1:(k-1))];  
        g35=g35+t(g35);    
        diag(g35) = (diag(g35)/2);  
        g36 = -invomg11 %*% g35 %*% invomg11;
        g36 = omega21 %*% g36 %*% g81;  
        g40 = cbind(g40,g36);    
        
        g45 = matrix(0,nrow=1,ncol=k-1);
        g45[l] = g15[l,k];
        g46 = cbind(g46,(g45 %*% invomg11 %*% g81));
  
        g50 = matrix(0,nrow=k-1,ncol=1);
        g50[l]=g20[l];
        g51 = cbind(g51,(omega21 %*% invomg11 %*% g50));        
      }
      g40 = cbind(g40,matrix(0,nrow=1,ncol=m-(k-1)));
      
      g49 = g15[k,1:k-1];
      g46 = cbind(g46,(g49 %*% invomg11 %*% g81));
      
      g47 = matrix(0,nrow=1,ncol=m);
      g47[(1:ncol(g46))] = g46;
      
      g51 = cbind(g51,matrix(0,nrow=1,ncol=m-(k-1)));

      gw = (g31+g40+g47+g51) * as.vector((condpass[1,j]/condk)*z[1,2]);
      
      #####   Start here for gradients with respect to rho parameters   ###
      kk = ncol(rhovec);  
      g60 = matrix(0,nrow=1,ncol=0);
      g65 = matrix(0,nrow=1,ncol=0);
      
      for(l in (1:kk))
      {
        g57 = matrix(0,nrow=m,ncol=m);
                
        g57[rRow[l],rCol[l]]=1;
        g57[rCol[l],rRow[l]]=1;
        
               
        g59 = g57*g25;
        
        g59 = g59[1:(k-1),1:(k-1)];
        g58 = -invomg11 %*% g59 %*% invomg11;
        g58 = omega21 %*% g58 %*% z2[1:(k-1)];
        g60 = cbind(g60,g58);   
        
        g61 = g57[k,1:(k-1)]*g25[k,1:(k-1)];
        g62 = g61 %*% invomg11 %*% z2[1:(k-1)];
        g65 = cbind(g65,g62); 
      }
      gr = (g60+g65) * as.vector(((condpass[j])/condk)*z[1,2]);
      gw1 = gw+gw1;     
      gr1 = gr+gr1;  
    }
    
    gw2 = (g30*condpass[1,j])+gw1;
    gr2 = (g53*condpass[1,j])+gr1;
    
    sk=0;
    gr3 = matrix(0,nrow=m,ncol=m);
    for(mm in (1:(m-1)))
    {
      gr3[mm,(mm+1):m] = (gr2[(sk+1):(sk+m-mm)]);
      sk=sk+m-mm;
    }
    ####    commands below to resequence gradients based on permutation   #####
       
    aaa = cbind(w[,j],(1:m));
    aaa = aaa[sort.list(aaa[,1]),2];
    gwf = gw2[aaa];
    grf = gr3[aaa,aaa];  

    res = matrix(nrow=0,ncol=2);
    for(ir in (1:(m-1)))
    {
      jr = ir+1;
      temp = matrix((rep(w[ir,j],m-ir)),ncol=1);     
      temp = cbind(temp,matrix(w[jr:m,j],nrow=m-jr+1,ncol=1));          
      res = rbind(res,temp);                              
    }
    res2 = cbind(apply(res,1,min),apply(res,1,max));
    
    res2 = matrix(res2[,1] *(10^(floor (log(res2[,2])/log(10))+1)) + res2[,2],ncol=1);  

    aaa1 = cbind(res2,(1:nrow(res2)));
    aaa1 = matrix(aaa1[sort.list(aaa1[,1]),2],ncol=1);

    grf = gr2[aaa1];    
    gwfinal = gwfinal+gwf;  
    grfinal = grfinal+grf;         
    
	}
	if( (p/n1) > 0)
	{
		 return(cbind(p/n1,cbind(gwfinal,grfinal)/n1,s1));
	}
	else
	{
		newval <- pdfmvna(a,r,s1);
		return(newval);
	}    
}


#' The cumulative density function of the mvna
#' 
#' @param a
#' @param r
#' @param s
#' @export
cdfmvna <- function(a,r,s){	
	m = length(a);
	w = 0; s1 = 0;	
	if(sys_randper ==0){
		w <-  getComb(m,2) ;
		s1 = s;
	}
	else{
		w = matrix(nrow=0,ncol=m);
		set.seed(s);
		for(z in 1:sys_randper)
		{	
			aa = cbind(seq(1:m),runif(m));
			aa = aa[order(aa[,2])];
			w <- rbind(w,aa);
		}
		s1 = ceiling(runif(1)*1000000000 );
	}
  #w = matrix(1:m,nrow=1); ## TEMP OVERRIDE FOR SAME MACML PERMUTATION
	w = t(w);
	n1 = ncol(w);
	p = 0;
	j = 1;

	a <- 	a * (a < 5.7) + 5.7*(a >= 5.7);
	ab = t(a);

	for(j in 1:n1)
	{
		x = ab[w[,j]];
		xvec = rep(x,m);
		rho = as.vector(r[w[,j],w[,j]]);		
		y = x[sapply(1:m,rep,m)];
		z = pbivnorm(xvec,y,rho);

		z = matrix(z,nrow=m,ncol=m);
		z3 = pnorm(x);

		z = z - diag(diag(z)) + diag(z3);

		z1 = matrix(pnorm(xvec) * pnorm(y),nrow=m,ncol=m);
		z2 = 1-z3;
		cond = 1;

		for(k in 3:m)
		{
			omega21 = z[k,1:k-1]-z1[k,1:k-1];
			omega11 = z[1:k-1,1:k-1]-z1[1:k-1,1:k-1];
			cm = chol2inv(chol(omega11));#chol2inv(chol/solve
			condk = z3[k] + omega21%*%(cm)%*%z2[1:k-1] ;   
			cond=cond*condk;
		}
		pcomb = z[1,2]*cond;
		p = p+pcomb;  
	}
	if( (p/n1) > 0)
	{
		 return(cbind(p/n1,s1));
	}
	else
	{
		newval <- cdfmvna(a,r,s1);
		return(newval);
	}
}

#' Get combination

#' @param n
#' @param dim
#' 
#' @export
getComb  <- function(n,dim){
	combsALL <- combn(n,dim);
	combsALL = t(combsALL);
	permsall = matrix(nrow=0,ncol=n);


	numCombs <- nrow(combsALL);
	for(irow in 1:numCombs)
	{
		othDims <- matrix(nrow=1,ncol=0);
		for(idim in 1:n)
		{
			check = which(combsALL[irow,]==idim);
			if(length(check) == 0 ){
				othDims <- cbind(othDims,idim);
			}
		}
		othDims <- permutations(ncol(othDims),ncol(othDims),othDims);
		iniDims = kronecker(combsALL[irow,], matrix(1,ncol=nrow(othDims),nrow=1)  );
		iniDims =t(iniDims );
		permsall <- rbind(permsall,cbind(iniDims,othDims));
	}
	return(permsall);
}