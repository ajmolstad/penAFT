###
### interior point program for Gehan lasso
###
### Brent Johnson
### 25 Sept 2008
###
library(hqreg)
library(quantreg)


gehan.lasso <- function(x,y,delta,lambda,CENTER=TRUE) {
  	ynew <- 1000*(length(y))^2
  	dimnum <- dim(x) 
  	n1 <- dimnum[1]
  	n2 <- dimnum[2]
  
  	if(CENTER) {
		my.means <- apply(x,2,mean)
		my.sd <- apply(x,2,sd)
		x <- t((t(x)-my.means)/my.sd)
		}
	
  	yy0 <- rep(y,rep(n1,n1))
  	delta1 <- rep(delta,rep(n1,n1))
 	yy1 <- rep(y,n1)
  	yy <- delta1*(yy0-yy1)
  	xx0 <- matrix(rep(as.vector(x),rep(n1,n1*n2)),nrow=n1*n1)
  	xx1 <- t(matrix(rep(as.vector(t(x)),n1),nrow=n2))
  	xx <- xx0-xx1   

  	xxdif <- xx*delta1 
  	xnew <- apply(xxdif,2,sum)
  	lamvec <- rep(lambda,n2)
  	xnew <- rbind(xxdif,-xnew,diag(n1*n1*lamvec))
  	yynew <- c(yy,ynew,rep(0,n2))

  	fit <- rq.fit.fnb(xnew,yynew,tau=0.5)
  	bhat <- fit$coef
  	
  	list(b=bhat)
}

gehan.loss <- function(x,y,delta,beta) {
	res <- y - as.vector(x %*% beta)
	tmp.fun <- function(I,res,del) {
		if(del[I] < 0.5) out <- 0
		else {
			ei <- res[I]
			diff <- ei - res
			out <- -1 * sum(diff[diff <= 0])
			}
		out/length(res)
		}
	vec <- unlist(lapply(1:length(y),tmp.fun,res,delta))
	sum(vec)
	}


loc.gehanlasso.cvfun <- function(MVEC,x.tr,y.tr,del.tr,x.te,y.te,del.te,CENTER) {
	
	if(CENTER) {
		my.means <- apply(x.tr,2,mean)
		my.sd <- apply(x.tr,2,sd)
		x.tr <- t((t(x.tr)-my.means)/my.sd)
		my.means <- apply(x.te,2,mean)
		my.sd <- apply(x.te,2,sd)
		x.te <- t((t(x.te)-my.means)/my.sd)
		}
	tmp.cvfun <- function(J,LVEC,x.tr,y.tr,del.tr,x.te,y.te,del.te,CENTER) {
		LAMBDA <- LVEC[J]
		ft.tr <- gehan.lasso(x.tr,y.tr,del.tr,lambda=LAMBDA,CENTER=CENTER)
		b.tr <- ft.tr$b
		#cat("b.tr = ",round(b.tr,2),"\n")
		lss <- gehan.loss(x.te,y.te,del.te,b.tr)
		#cat("lss = ",round(lss,3),"\n")
		lss
		}
	loc.out <- unlist(lapply(1:length(MVEC),tmp.cvfun,MVEC,x.tr,y.tr,del.tr,x.te,y.te,del.te,CENTER))
	#cat("loc.out = ",round(loc.out,2),"\n")
	loc.out
	}

cv.gehan.lasso <- function(x,y,delta,lamvec,V=10,CENTER=TRUE) {
	if(!any(search()=="package:quantreg"))
    stop("cv.gehan.lasso:  The 'quantreg' package is not loaded.")

	#if(!any(search()=="package:survival"))
    #stop("cv.gehan.lasso:  The 'survival' package is not loaded.")
	
	n <- length(y)
	s <- sample(1:n,n,replace=FALSE)
	
	LOSS.MAT <- matrix(numeric(0),V,length(lamvec))
	for (J in 1:V) {
		cat("J = ",J,"\n")
		te.ind <- s[(1+(J-1)*floor(n/V)):(J*floor(n/V))]
		x.tr <- x[-te.ind,]
		y.tr <- y[-te.ind]
		del.tr <- delta[-te.ind]
		x.te <- x[te.ind,]
		y.te <- y[te.ind]
		del.te <- delta[te.ind]
		LOSS.MAT[J,] <- loc.gehanlasso.cvfun(lamvec,x.tr,y.tr,del.tr,x.te,y.te,del.te,CENTER=CENTER)
		}
	
	lvec <- apply(LOSS.MAT,2,mean)
	sdvec <- apply(LOSS.MAT,2,sd)
	out <- cbind(lamvec,lvec,sdvec)
	m.cv <- lamvec[order(lvec)[1]]
	list(out=out,mcv=m.cv)
	
	}
	
	