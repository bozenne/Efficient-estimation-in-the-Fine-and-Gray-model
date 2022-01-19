
library(mets)

##########################################################################
### NPMLE using stratified model for F_2 and FG for F_1 
##########################################################################

mledoubleFGR <- function(formula,data,offset=NULL,weights=NULL,X2=NULL,...) {# {{{
  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
  if (ncol(Y)==2) {
    exit <- Y[,1]
    entry <- NULL ## rep(0,nrow(Y))
    status <- Y[,2]
  } else {
    entry <- Y[,1]
    exit <- Y[,2]
    status <- Y[,3]
  }
  id <- strata <- NULL
  if (!is.null(attributes(Terms)$specials$cluster)) {
    ts <- survival::untangle.specials(Terms, "cluster")
    pos.cluster <- ts$terms
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  } else pos.cluster <- NULL
  if (!is.null(attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    pos.strata <- ts$terms
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name <- ts$vars
  }  else { strata.name <- NULL; pos.strata <- NULL}
###  if (!is.null(attributes(Terms)$specials$offset)) {
###    ts <- survival::untangle.specials(Terms, "offset")
###    pos.offset <- ts$terms
###    Terms  <- Terms[-ts$terms]
###    offset <- m[[ts$vars]]
###  }  else pos.offset <- NULL
  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
###  X2call <- X2;
###  if (!is.null(X2)) {
###	 if (X2[1]==0)  X2 <- matrix(nrow=0,ncol=0)
###	 else { X2 <-  X[,X2call,drop=FALSE];
###                X <- X[,-X2call,drop=FALSE];
###	 }
###  } else X2 <- X
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
###  if (ncol(X2)==0) X2 <- matrix(nrow=0,ncol=0)
  res <- c(mledoubleFG01R(X,entry,exit,status,id,strata,offset,weights,strata.name,X2call=X2call,...),
   list(call=cl,model.frame=m,formula=formula,strata.pos=pos.strata,cluster.pos=pos.cluster))
  class(res) <- c("doubleFG")

  res
}# }}}

mledoubleFG01R <- function(X,entry,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,
  strata.name=NULL,cumhaz=TRUE,
  beta,stderr=TRUE,method="nlm",no.opt=FALSE,Z=NULL,propodds=NULL,AddGam=NULL,
	     restrict=0,case.weights=NULL,X2call=NULL,...) {# {{{
  p <- p1 <- ncol(X);
  if (missing(beta))  beta <- NULL 
  if (p1==0) X <- cbind(rep(0,length(exit)))
  if (is.null(strata)) { strata <- rep(0,length(exit)); nstrata <- 1; strata.level <- NULL; } else {
	  strata.level <- levels(strata)
	  ustrata <- sort(unique(strata))
	  nstrata <- length(ustrata)
	  strata.values <- ustrata
      if (is.numeric(strata)) strata <-  fast.approx(ustrata,strata)-1 else  {
      strata <- as.integer(factor(strata,labels=seq(nstrata)))-1
    }
  }
  if (is.null(offset)) offset <- rep(0,length(exit))
  if (is.null(weights)) weights <- rep(1,length(exit))
  strata.call <- strata
  Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given
  if (!is.null(Z)) Zcall <- Z

  ## possible casewights to use for bootstrapping and other things
  if (is.null(case.weights)) case.weights <- rep(1,length(exit))

  trunc <- (!is.null(entry))
  if (!trunc) entry <- rep(0,length(exit))

  if (!is.null(id)) {
	  ids <- unique(id)
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(entry))-1;
   ## orginal id coding into integers
   id.orig <- id+1;

   dd <- .Call("FastCoxPrepStrata",entry,exit,status,X,id,trunc,strata,weights,offset,Zcall,case.weights,PACKAGE="mets")
###   print( table(dd$status))

   n1 <- sum(dd$status==1)
   n2 <- sum(dd$status==2)
   w1 <- which(dd$status==1)
   w2 <- which(dd$status==2)
   w0 <- which(dd$status==0)

   status <- c(dd$status)
   strata <- c(dd$strata)
   dd$nstrata <- nstrata
   if (is.null(beta)) 
   beta <- c(0,0,log(rep(1/n1,n1)),log(rep(1/n2,n2)))

###   print(c(n1,n2)); print(length(beta)); print(length(w1)); print(length(w2));

	obj <- function(pp,U=FALSE,all=FALSE) {# {{{

		beta <- pp[1:2]
		dLam1 <- rep(0,length(status))
		dLam2 <- rep(0,length(status))
		dLam1[w1] <- exp(pp[3:(3+n1-1)])
		dLam2[w2] <- exp(pp[-(1:(3+n1-1))])
                Lam1 <- cumsum(dLam1)
                Lam2 <- cumsumstrata(dLam2,strata,nstrata)
		Lam1tau <- tail(Lam1,1)
		lp <- X %*% beta
		RR <- exp(lp)

	ldF1 <- sum(log(dLam1[w1])+lp[w1]-Lam1[w1]*RR[w1])
	ldF2 <- sum(log(dLam2[w2])-Lam2[w2]-Lam1tau*RR[w2])

psurv <- (1-(1-exp(-Lam1[w0]*RR[w0]))-(1-exp(-Lam2[w0]))*(exp(-Lam1tau*RR[w0])))
###		print(summary(psurv))
	lpsurv <- sum(log(psurv))
###		print(c(ldF1,ldF2,lpsurv))
ploglik <-ldF1+ldF2+lpsurv
###	print(ploglik)
return(-ploglik)
}# }}}

opt <- nlm(obj,beta,...)
opt$method <- "nlm"

opt$coef <- opt$estimate[1:2]
opt$dLam1 <- exp(opt$estimate[3:(3+n1-1)])
opt$Lam1 <- cumsum(opt$dLam1)
opt$dLam2 <- exp(opt$estimate[-(1:(3+n1-1))])
opt$strata <- dd$strata[w2]
     opt$Lam2 <- cumsumstrata(opt$dLam2,opt$strata,nstrata)
  return(opt)

}# }}}

##########################################################################

###First we simulate some competing risks data using some utility functions.
###
###We simulate data with two causes based on the Fine-Gray model:
###\begin{align}
###F_1(t,X) &  = P(T\leq t, \epsilon=1|X)=( 1 - exp(-\Lambda_1(t) \exp(X^T \beta_1))) \\
###F_2(t,X) & = P(T\leq t, \epsilon=2|X)= ( 1 - exp(-\Lambda_2(t) \exp(X^T \beta_2))) 
###              \cdot (1 - F_1(\infty,X))  
###\end{align}
###where the baselines are given as $\Lambda_j(t) = \rho_j (1- exp(-t/\nu_j))$ for $j=1,2$, and the 
###$X$ being two independent binomials. 

### simulation of data with  binary covariates 
library(mets)
depcens <- 0
set.seed(100)
rho1 <- 0.2; rho2 <- 5.9
n <- 800
beta=c(0.0,-0.1,-0.5,0.3)
dep=0
rc <- 0.5
###
dats <- simul.cifs(n,rho1,rho2,beta,depcens=0,bin=1,rc=rc, rate=c(1,7))
cif1 <- cifreg(Event(time,status)~Z1+Z2,dats,propodds=NULL)
###
dsort(dats) <- ~time

npmle1 <- mledoubleFGR(Event(time,status)~Z1+Z2+strata(Z1,Z2),dats,method="nlm")
npmle1$coef

######################################################################################

library(doMC)
cc=detectCores()
registerDoMC(cc)

### simulations and fitting of FG-model, augmented model, and FG-CM (with strataified censoring)

onerun <- function(i,n,mmle=0,ddfg=0,depcens=1,rc=0.5,bin=1,res=0,no.opt=FALSE,lbeta=NULL) {# {{{

if (i%%500==0) print(i)

if (length(bin)==1) bin <- rep(bin,2)
if (is.null(lbeta)) lbeta <- beta 

 dats <- simul.cifs(n,rho1,rho2,beta,depcens=depcens,bin=bin,rc=rc)
 dsort(dats) <- ~time

 if (bin[1]==0) dcut(dats,breaks=2) <- Z1g~Z1 else dats$Z1g <- dats$Z1
 if (bin[2]==0) dcut(dats,breaks=2) <- Z2g~Z2 else dats$Z2g <- dats$Z2
if (sum(dats$status==1)>=10) {
 ### possibly biased due to censoring dependence 
 fg <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL,beta=lbeta[1:2],no.opt=no.opt)

 fgcm <- tryCatch( cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL,
			  cens.model=~strata(Z1g,Z2g)) ,error=function(x) NULL) 
 cr <- phreg(Surv(time,status==0)~+1,data=dats)
dfgwF2 <- dfg <- fgaug2 <-  fgaug <- NULL
if (ddfg==1)
 dfg  <- tryCatch( doubleFGR(Event(time,status)~Z1+Z2,data=dats,restrict=res),error=function(x) NULL) 

mle <- fgsaug2 <-  fgsaug <- fgsaug3 <-  NULL
if (mmle==1)
mle <- tryCatch( mledoubleFGR(Event(time,status)~Z1+Z2+strata(Z1,Z2),dats,method="nlm"),error=function(x) NULL) 

if (!is.null(mle)) { mle$coef <- mle$coef;  mle$se.coef <- rep(NA,2)}

if (!is.null(fg)) fgsaug <-tryCatch( 
    FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1g,Z2g),data=dats,E=fg$E,cause=1) 
    ,error=function(x) NULL) 

 if (!is.null(fgsaug)) fgsaug2 <- tryCatch(
    FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1g,Z2g),data=dats,E=fgsaug$E,cause=1) 
    , error=function(x) NULL) 

if (!is.null(fgsaug)) fgsaug3 <- tryCatch(
    FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1g,Z2g),data=dats,E=fgsaug2$E,cause=1) 
    , error=function(x) NULL) 


 if (is.null(fgcm)) fgcm <- list(coef=rep(NA,2),se.coef=rep(NA,2))
 if (is.null(mle)) mle <- list(coef=rep(NA,2),se.coef=rep(NA,2))
 if (is.null(dfg)) dfg <- list(coef=rep(NA,2),se.coef=rep(NA,2))
 if (is.null(dfgwF2)) dfgwF2 <- list(coef=rep(NA,2),se.coef=rep(NA,2))
 if (is.null(fgaug)) fgaug <- list(coef=rep(NA,2),se.coef=rep(NA,2))
 if (is.null(fgaug2)) fgaug2 <- list(coef=rep(NA,2),se.coef=rep(NA,2))
 if (is.null(fgsaug)) fgsaug <- list(coef=rep(NA,2),se.coef=rep(NA,2))
 if (is.null(fgsaug2)) fgsaug2 <- list(coef=rep(NA,2),se.coef=rep(NA,2))
 if (is.null(fgsaug3)) fgsaug3 <- list(coef=rep(NA,2),se.coef=rep(NA,2))

 res <- list(
	       fg=fg$coef      ,se.fg   =    fg$se.coef     ,
	       se.fg1=fg$se1.coef,
	       dfg=dfg$coef    ,se.dfg   =  dfg$se.coef     ,
	     fgcm=fgcm$coef    ,se.fgcm =  fgcm$se.coef     ,
	      aug=fgaug$coef   ,se.aug  = fgaug$se.coef     ,
	     aug2=fgaug2$coef  ,se.aug2 = fgaug$se.coef     ,
	     saug=fgsaug$coef  ,se.saug =fgsaug$se.coef     ,
	    saug2=fgsaug2$coef ,se.saug2=fgsaug$se.coef     ,
	    saug3=fgsaug2$coef ,se.saug3=fgsaug$se.coef     ,
	      mle=mle$coef[1:2],se.mle  =   mle$se.coef[1:2],
	       dfgwF2=dfgwF2$coef[1:2])
} else res <- list(
	       fg=rep(NA,2),se.fg   =rep(NA,2),
	       se.fg1=rep(NA,2),
	       se.dfg=rep(NA,2),
	     fgcm=rep(NA,2),se.fgcm =rep(NA,2),
	      aug=rep(NA,2),se.aug  =rep(NA,2),
	     aug2=rep(NA,2),se.aug2 =rep(NA,2),
	     saug=rep(NA,2),se.saug =rep(NA,2),
	    saug2=rep(NA,2),se.saug2=rep(NA,2),
	    saug3=rep(NA,2),se.saug3=rep(NA,2),
	      mle=rep(NA,2),se.mle  =rep(NA,2),
	      dfg=rep(NA,2),dfgwF2  =rep(NA,2))


 return(res)
}# }}}

onerun(1,800,mmle=1)

## functions to  handle output 

ana <- function(res) {# {{{

fg <- do.call("rbind",lapply(res,function(x) x$fg))
mle <- do.call("rbind",lapply(res,function(x) x$mle))
###dfg <- do.call("rbind",lapply(res,function(x) x$dfg))
###dfgwF2 <- do.call("rbind",lapply(res,function(x) x$dfgwF2))
fgcm <- do.call("rbind",lapply(res,function(x) x$fgcm))
###aug <- do.call("rbind",lapply(res,function(x) x$aug))
###aug2 <- do.call("rbind",lapply(res,function(x) x$aug2))
saug <- do.call("rbind",lapply(res,function(x) x$saug))
saug2 <- do.call("rbind",lapply(res,function(x) x$saug2))
saug3 <- do.call("rbind",lapply(res,function(x) x$saug3))
if (is.null(saug3)) saug3 <- matrix(NA,10,2)
bres <- cbind(
apply(fg,2,mean,na.rm=TRUE),
apply(mle,2,mean,na.rm=TRUE),
###apply(dfg,2,mean,na.rm=TRUE),
###apply(dfgwF2,2,mean,na.rm=TRUE),
apply(fgcm,2,mean,na.rm=TRUE),
###apply(aug,2,mean,na.rm=TRUE),
###apply(aug2,2,mean,na.rm=TRUE),
apply(saug,2,mean,na.rm=TRUE),
apply(saug2,2,mean,na.rm=TRUE),
apply(saug3,2,mean,na.rm=TRUE),
beta[1:2])
colnames(bres) <- c("fg","mle",##"dfg","dfgwF2",
		    "fgmc",
###		    "fg-aug","fg-aug2",
		    "fg-Saug",
		    "fg-Saug2",
		    "fg-Saug3",
		    "truth")
###bres

###fgscore <- do.call("rbind",lapply(res,function(x) x$fgscore))
###augaug <- do.call("rbind",lapply(res,function(x) x$augaug1))
###augaug2 <- do.call("rbind",lapply(res,function(x) x$augaug2))
###augaug2 <- augaug
###augaug <- cbind( 
###apply(fgscore,2,mean,na.rm=TRUE),
###apply(augaug,2,mean,na.rm=TRUE),
###apply(augaug2,2,mean,na.rm=TRUE))
###colnames(augaug) <- c("fgscore","augmentation1","augmentation2")

nabres <- cbind(
apply(is.na(fg),2,mean,na.rm=TRUE),
apply(is.na(mle),2,mean,na.rm=TRUE),
###apply(is.na(dfg),2,mean,na.rm=TRUE),
apply(is.na(fgcm),2,mean,na.rm=TRUE),
###apply(is.na(aug),2,mean,na.rm=TRUE),
###apply(is.na(aug2),2,mean,na.rm=TRUE),
apply(is.na(saug),2,mean,na.rm=TRUE),
apply(is.na(saug2),2,mean,na.rm=TRUE),
apply(is.na(saug3),2,mean,na.rm=TRUE),
nrow(fg))
colnames(nabres) <- c("fg","mle",##"dfg",
		      "fgcm",##"fg-aug","fg-aug2",
		      "fg-saug", "fg-saug2","fg-saug3","truth")
# nabres

sds <- cbind(
   apply(fg,2,sd,na.rm=TRUE),
###  apply(dfg,2,sd,na.rm=TRUE),
 apply(fgcm,2,sd,na.rm=TRUE),
###  apply(aug,2,sd,na.rm=TRUE),
###  apply(aug2,2,sd,na.rm=TRUE),
  apply(saug,2,sd,na.rm=TRUE),
  apply(saug2,2,sd,na.rm=TRUE),
   apply(mle,2,sd,na.rm=TRUE),
  apply(saug3,2,sd,na.rm=TRUE)
)
colnames(sds) <- c("fg",## "dfg",
		   "fgcm",
###		   "fg-aug","fg-aug2",
		   "fg-Saug","fg-Saug2",
		   "mle","fg-Saug3")
###print(sds)
sds

varrel <- cbind(sds[,1]/sds[,3],sds[,1]/sds[,4],sds[,1]/sds[,5])^2
colnames(varrel) <- c("fg/fg-Saug","fg/fg-Saug2","fg/mle")

rsefg <- do.call("rbind",lapply(res,function(x) x$se.fg))
rsedfg <- do.call("rbind",lapply(res,function(x) x$se.dfg))
if (is.null(rsedfg)) rsedfg <- matrix(NA,10,2)
rsefg1 <- do.call("rbind",lapply(res,function(x) x$se.fg1))
rsefgcm <- do.call("rbind",lapply(res,function(x) x$se.fgcm))
rsemle <- do.call("rbind",lapply(res,function(x) x$se.mle))
rsefgs <- do.call("rbind",lapply(res,function(x) x$se.saug))
rsefgs2 <- do.call("rbind",lapply(res,function(x) x$se.saug2)) 
###    if (length(x$se.saug2)==1) y <- rep(x$se.aug2,2) else y  <-  x$se.aug; print(y); return(y)} ))
ese <- cbind(
    apply(rsefg,2,mean,na.rm=TRUE),
    apply(rsefg1,2,mean,na.rm=TRUE),
###    apply(rsedfg,2,mean,na.rm=TRUE),
    apply(rsefgcm,2,mean,na.rm=TRUE),
    apply(rsemle,2,mean,na.rm=TRUE),
    apply(rsefgs,2,mean,na.rm=TRUE),
    apply(rsefgs2,2,mean,na.rm=TRUE)
   )
colnames(ese) <- c("se-fg","se-fg-naive","se-dfg",
		   "se-fgcm",##"se-mle",
		   "se-fgtS","se-fgS2")
ese

varrel <- cbind(varrel,(ese[,2]/ese[,1]))
colnames(varrel)[4] <- "(ese-fg-naive/ese-fg)"


covf <- function(mu,semu,truth=beta[1:2]) 
{
 up <- mu+1.96*semu
 lo <- mu-1.96*semu
 icovv <- (t(up) > truth) & (t(lo) < truth)
 icovv <- apply(icovv,1,mean,na.rm=TRUE)
 return(icovv)
}

covv <- cbind(
covf(fg,rsefg),
covf(fg,rsefg1),
covf(fgcm,rsefgcm),
covf(mle,rsemle),
covf(saug,rsefgs),
covf(saug2,rsefgs2))
colnames(covv) <- c("cov-fg","cov-fg-naive","cov-fgcm","cov-mle","cov-fgtS","cov-fgS2")

list(bias=round(bres,2),na=nabres,sds=round(sds,4),
     varrel=round(varrel,4),ese=round(ese,4),covv=round(covv,3))

}# }}}

rell <- function(restot,ns=c(200,400,800),rho1s=c(0.2,0.4),rho2s=c(1,10))
{# {{{
#

tot <- c()
k <- 0
for (n in ns)
for (rho1 in rho1s)
for (rho2 in rho2s) 
for (rc in seq(crate)) {
cens <- c(10,25,40)[rc]
k <- k+1
ccc <- restot[[k]]$mm$status[1]
tot <- rbind(tot,c(restot[[k]]$n,cens,ccc,rho1,rho2,
		   restot[[k]]$ana$varrel[1:2,5],restot[[k]]$ana$varrel[1:2,3]))
}
tot <- data.frame(tot)
names(tot) <- c("n","cens","censo","rho1","rho2","NPMPLE-1","NPMLE-2","AUG-1","AUG-2")
dsort(tot) <- ~cens+rho1+rho2+n
rownames(tot) <- NULL
tot

}# }}}

coverage <- function(restot) 
{# {{{
totbias <- c()
k <- 0
for (n in c(200,400,800))
for (rho1 in c(0.2,0.4))
for (rho2 in c(1,10)) 
for (rc in seq(crate)) {
cens <- c(10,25,40)[rc]
k <- k+1
ccc <- restot[[k]]$mm$status[1]
totbias <- rbind(totbias,c(restot[[k]]$n,cens,ccc,rho1,rho2,
		   restot[[k]]$ana$covv[1:2,1], restot[[k]]$ana$covv[1:2,3],
		   restot[[k]]$ana$covv[1:2,5]))
}
totbias
###
totbias <- data.frame(totbias)
names(totbias) <- c("n","cens","censo","rho1","rho2",
		paste("FG",1:2,sep=""), paste("FG-CM",1:2,sep=""), paste("FG-AUG",1:2,sep="")
		)
dsort(totbias) <- ~cens+rho1+rho2+n
###
totbias[,-3]
}# }}}

bias <- function(restot) 
{# {{{
totbias <- c()
k <- 0
for (n in c(200,400,800))
for (rho1 in c(0.2,0.4))
for (rho2 in c(1,10)) 
for (rc in seq(crate)) {
cens <- c(10,25,40)[rc]
k <- k+1
ccc <- restot[[k]]$mm$status[1]
totbias <- rbind(totbias,c(restot[[k]]$n,cens,ccc,rho1,rho2,
		   restot[[k]]$ana$bias[1:2,1]+c(0,0.1),
		   restot[[k]]$ana$bias[1:2,2]+c(0,0.1),
		   restot[[k]]$ana$bias[1:2,5]+c(0,0.1),
		   restot[[k]]$ana$bias[1:2,8]+c(0,0.1),
		   restot[[k]]$ana$bias[1:2,9]+c(0,0.1),
		   restot[[k]]$ana$bias[1:2,10]+c(0,0.1)
		   ))
}
totbias
###
totbias <- data.frame(totbias)
names(totbias) <- c("n","cens","censo","rho1","rho2",
		paste("FG",1:2,sep=""),
		paste("NPMLE",1:2,sep=""),
		paste("FG-CM",1:2,sep=""),
		paste("FG-AUG",1:2,sep=""),
		paste("FG-AUG2",1:2,sep="")
		)
dsort(totbias) <- ~cens+rho1+rho2+n
###
totbias[,-3]
}# }}}

### simulation test
nsim <- 100
rho1 <- 0.4; rho2 <- 8; 
res400d02 <- foreach (i=0:n.sim) %dopar% onerun(i,800,mmle=1,rc=1,depcens=0,res=3)
ana(res400d02)


####################################################################################
### running all models with independent censoring and storing output ###############
###  with different censorings rates, by crate                       ############### 
####################################################################################

rate <- c(0.035,0.12,0.55,0.03,0.1,0.6)
crate <- c(0.035,0.12,0.3)

## independent censoring situation
nsim <- 100
restot <- list()
k <- 0
for (n in c(200,400,800))
for (rho1 in c(0.2,0.4))
for (rho2 in c(1,10)) 
for (rc in crate) {
k <- k+1
print(c(k,rho1,rho2,rc))
rrc <- (rho1+rho2)*rc
dats <- simul.cifs(n,rho1,rho2,beta,depcens=0,bin=1,rc=rrc)
mm <- dtable(dats,status~I(time<6),prop=1)$"TRUE"$table
res <- foreach (i=0:nsim) %dopar% onerun(i,n,rc=rrc,depcens=0,res=3,ddfg=0,mmle=1)
anao=ana(res)
print(c(k,rho1,rho2,rc))
print(anao)
restot[[k]] <- list(ana=anao,n=n,rho1=rho1,rho2=rho2,rc=rc,mm=mm)
}



## dependent censoring situation
restot <- list()
k <- 0
for (n in c(200,400,800))
for (rho1 in c(0.2,0.4))
for (rho2 in c(1,10)) 
for (rc in crate) {
k <- k+1
print(c(k,rho1,rho2,rc))
rrc <- (rho1+rho2)*rc
if (rho2==10 & rho1==0.2 & rc==0.035) rrc <- 1.5*rrc
dats <- simul.cifs(n*100,rho1,rho2,beta,depcens=0,bin=1,rc=rrc)
###mm1 <- dtable(dats,~status,prop=1)$table
mm <- dtable(dats,status~I(time<6),prop=1)$"TRUE"$table
###print(mm$status)
res <- foreach (i=0:nsim) %dopar% onerun(i,n,rc=rrc,depcens=1,res=3,ddfg=0,mmle=1)
anao=ana(res)
print(c(k,rho1,rho2,rc))
print(anao)
restot[[k]] <- list(ana=anao,n=n,rho1=rho1,rho2=rho2,rc=rc,mm=mm)
}


