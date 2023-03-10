############################################################################################
bivlocsc2stg <- function(formula1, formula2, formula3, weights=NULL, xstrata=NULL, data, 
                           obstype.name="obstype", fit=TRUE, xs.includes=FALSE, off.set=NULL,
			   errdistrn="normal", errmodpars=6, start=NULL, Qstart=NULL, 
                           control=mlefn.control(...), control.inner=mlefn.control.inner(...), ...)
############################################################################################
{
# Process formula -------------------------------------------------------------------
# formula1 of the form y1~y2+x1+x2..., formula2 of the form y2~x1+x2..., formula3 of the form ~x1+x2...
 
mf <- call <- match.call()
mf[[1]] <- as.name("model.frame")
names(mf)[1] <- "model"
mf$na.action <- as.name("na.keep")
NULL -> mf$xstrata -> mf$obstype.name -> mf$fit -> mf$xs.includes -> mf$off.set
NULL -> mf$errdistrn -> mf$errmodpars -> mf$start -> mf$Qstart 
NULL -> mf$control -> mf$control.inner


# First formula processing ----------------------------------------------------------
mf1 <- mf[c("model","formula1","data","weights")]
mf1$na.action <- as.name("na.keep")
resp1 <- as.character(mf1$formula[2])
names(mf1)[2]<- "formula"
mf1 <- eval(mf1, sys.frame(sys.parent()))
Terms <- y1vname <- attr(mf1,"terms")
X1 <- model.matrix(Terms,mf1)
if(!is.matrix(X1)) X1 <- matrix(X1)
X1Orig <- X1
w <- model.extract(mf1, weights)
terms1 <- attr(terms(formula1),"term.labels")
order1 <- attr(terms(formula1),"order")
assign1 <- attr(X1,"assign")
fnames1<-names(attr(X1,"contrasts"))
y1 <- model.extract(mf1,response)
names(y1) <- NULL



# Second formula processing ---------------------------------------------------------
mf2 <- mf[c("model","formula2","data")]
mf2$na.action <- as.name("na.keep")
resp2 <- as.character(mf2$formula[2])
names(mf2)[2] <- "formula"
mf2 <- eval(mf2, sys.frame(sys.parent()))
Terms <- y2vname <- attr(mf2,"terms")
X2 <- model.matrix(Terms,mf2)
if(!is.matrix(X2)) X2 <- matrix(X2)
X2Orig <- X2
terms2 <- attr(terms(formula2),"term.labels")
order2 <- attr(terms(formula2),"order")
assign2 <- attr(X2,"assign")
fnames2<-names(attr(X2,"contrasts"))
y2 <- model.extract(mf2,response)
names(y2) <- NULL


# Third formula processing ---------------------------------------------------------
mf3 <- mf[c("model","formula3","data")]
mf3$na.action <- as.name("na.keep")
names(mf3)[2] <- "formula"
mf3$formula <- as.formula(paste(resp2,"~",mf3$formula[2]))
mf3 <- eval(mf3, sys.frame(sys.parent()))
Terms <- attr(mf3,"terms")
X3 <- model.matrix(Terms,mf3)
if(!is.matrix(X3)) X3 <- matrix(X3)
terms3 <- attr(terms(formula3),"term.labels")
order3 <- attr(terms(formula3),"order")
assign3 <- attr(X3,"assign")
fnames3 <- names(attr(X3,"contrasts"))	
 

# -----------------------------------------------------------------------------------
if (!is.R())  xtabs <- xtabsforS  ### So can run in Splus ###

n <- if (is.matrix(y2)) dim(y2)[1] else length(y2)
if (is.null(w)) w <- rep(1,n)


# Input checking ---------------------------------------------------------------------

if (is.na(match(obstype.name,names(data)))) 
    stop(paste("Dataframe did not have a column named",obstype.name))
obstype <- data[,match(obstype.name,names(data))]
if (any(is.na(data[,obstype.name])))
    stop(paste("Missing values not permitted in obstype variable", obstype.name))
    
prospective <- all(obstype %in% c("uncond","y|x"))
# for S may need prospective <- all(!is.na(match(obstype,c("uncond","y|x"))))

if (any(is.na(w))) stop("Weights should not contain NAs")


# Form the off.set -----------------------------------------------------------------
if (is.null(off.set)) off.set <- matrix(0, n, 3)  
else {
   off.set <- as.matrix(off.set)
   if (dim(off.set)[1] != n) 
   	stop("off.set should have the same row.length as n !")
   if (dim(off.set)[2] != 3)
   	stop("off.set should have three columns !")
}


# Set up X4-matrix with only those variables with interaction with Y2 --------------
# An 1 vector is used for Y2-variable itself

prosp <- function(formula, data) {
  mf <- call <- match.call()
  mf[[1]] <- as.name("model.frame")
  names(mf)[1] <- "model"
  mf$na.action <- as.name("na.keep")
  names(mf)[2] <- "formula"
  mf <- eval(mf, sys.frame(sys.parent()))
  Terms <- attr(mf,"terms")
  X <- model.matrix(Terms,mf)
  if(!is.matrix(X)) X <- matrix(X)
  X 
}

y2id <- match(resp2, names(data))
noNA <- apply(!is.na(data),1,sum)
tempdata <- data[noNA==dim(data)[2],]

tempdata[,y2id] <- rep(max(y2,na.rm=TRUE)+1,dim(tempdata)[1])
y2posiny1 <- NULL
X1temp <- prosp(formula1,data=tempdata)
for (j in 1:dim(X1temp)[2]) {
   if (any((X1temp[,j]-X1[noNA==dim(data)[2],j])!=0))
      y2posiny1 <- c(y2posiny1, j)
}

tempdata <- data
tempdata[,y2id] <- rep(1, dim(data)[1])
X4 <- prosp(formula1, data=tempdata)
X4 <- X4[,y2posiny1,drop=FALSE]


# Form the x-strata -----------------------------------------------------------------

if (!is.null(xstrata)) { ## Set Null xstrata if there is only one stratum exists !
  vx <- as.matrix(data[,match(xstrata,names(data))])
  vx1 <- vx[which(apply(!is.na(vx),1,sum)==dim(vx)[2]),,drop=FALSE]
  if (dim(unique(vx1))[1]==1) xstrata <- NULL }

xStrat <- if (is.null(xstrata)) rep(1,n) else {
    xstratform <- as.formula(paste("~","-1 + ",paste(xstrata,sep="",collapse=":")))
    data1 <- as.data.frame(lapply(data, factor))
    xstratmat <- model.matrix(xstratform, model.frame(xstratform,data1,na.action=function(x)x))
    for(i in 1:ncol(xstratmat)) xstratmat[,i] <- i * xstratmat[,i]
    apply(xstratmat,1,sum)
}
names(xStrat) <- NULL
nStrat <- max(xStrat)


# Form the y1-vector -------------------------------------------------------
matrixY1 <- FALSE

if ((is.vector(y1)|is.factor(y1))&!is.list(y1)) {
    y1fac<-factor(y1)
    if (length(levels(y1fac))>1){
        y1 <- 1 + 1*((y1fac)==levels(y1fac)[1])
        y1name <- if (length(levels(y1fac))==2)  as.character(levels(y1fac))[2]
                  else paste("Not.",as.character(levels(y1fac))[1],sep="")
    }
    else if (length(levels(y1fac))==1){
        y1 <- 1*((y1fac)==levels(y1fac)[1])
        y1name <- c(as.character(levels(y1fac)),NA)
    }
    else stop("illegal Y1-data")
}
else if (is.matrix(y1)) {#in form of (case,control), expand out into vector form
    matrixY1 <- TRUE
    y1name <- colnames(y1)[1]
    w <- as.vector(y1*w)
    X1 <- rbind(X1,X1)
    X2 <- rbind(X2,X2)
    X3 <- rbind(X3,X3)
    X4 <- rbind(X4,X4)

    y1 <- rep(c(1,2), rep(n,2))
    y2 <- rep(y2, 2)
    off.set <- rbind(off.set, off.set)
    obstype <- rep(obstype, 2)
    xStrat <- rep(xStrat, 2)
    n <- 2*n
}
else stop("Illegal Y1-data")


# Form y-matrix and y1strata ---------------------------------------------
y <- cbind(2-y1, y2) ###make Y1 values to be 1/0 !
dimnames(y) <- list(NULL, c(resp1, resp2))

yStrat <- y1; nyStrat <- max(y1) 


# Remove random missing values --------------------------------------------------
mat <- cbind(y,X1,X2,X3,xStrat)
miss1 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="uncond"]
miss2 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="retro"]
miss3 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="y|x"]

mat <- cbind(X1,X2,X3,xStrat)
miss4 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="xonly"]

mat <- cbind(y[,1],xStrat) 
miss5 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="strata"]

toremove <- c(miss1,miss2,miss3,miss4,miss5)
missReport <- NULL
if(length(toremove)>0){
   #cat("\nObservations deleted due to missing data:\n")
   if ((length(miss1)>0))
   	missReport <- rbind(missReport, paste("uncond:",length(miss1),"rows relating to", 
		      	sum(w[miss1]),"observations"))
   if ((length(miss2)>0))
	missReport <- rbind(missReport, paste("retro:",length(miss2),"rows relating to",
			sum(w[miss2]),"observations"))
   if ((length(miss3)>0))
	missReport <- rbind(missReport, paste("y|x:",length(miss3),"rows relating to",
			sum(w[miss3]),"observations"))
   if ((length(miss4)>0))
	missReport <- rbind(missReport,paste("xonly:",length(miss4),"rows relating to",
			sum(w[miss4]),"observations"))
   if ((length(miss5)>0))
	missReport <- rbind(missReport,paste("strata:",length(miss5),"rows relating to",
			sum(w[miss5]),"observations"))
   dimnames(missReport) <- list(rep("",nrow(missReport)),rep("",ncol(missReport)))

   w <- w[-toremove]
   obstype <- factor(as.character(obstype[-toremove]))
   xStrat <- xStrat[-toremove]
   yStrat <- yStrat[-toremove]
   y <- y[-toremove,,drop=FALSE]
   X1 <- X1[-toremove,,drop=FALSE]
   X2 <- X2[-toremove,,drop=FALSE]
   X3 <- X3[-toremove,,drop=FALSE]
   X4 <- X4[-toremove,,drop=FALSE]
   off.set <- off.set[-toremove,,drop=FALSE]
   n <- length(w)
}
yStratfac <- factor(yStrat)
xStratfac <- factor(xStrat)
obstypefac <- factor(as.character(obstype))


# Stratum Counts Report -------------------------------------------------------------
#cat("\nStrata Counts Report:\n")

StrReport <- xStrReport <- key <- NULL
if (!matrixY1)
     yKey <- c(as.character(y1vname[[2]]), y1name, as.character(y2vname[[2]]))
else yKey <- c("y1", y1name, as.character(y2vname[[2]]))

sub <- (1:n)[obstype=="xonly"]
if(length(sub)>0){
    df <- data.frame(w, obstype=obstypefac, yStrat=yStratfac, xStrat=xStratfac)[-sub,]
    df1 <- data.frame(w, obstype=obstypefac, yStrat=yStratfac, xStrat=xStratfac)[sub,]
}
else df <- data.frame(w, obstype=obstypefac, yStrat=yStratfac, xStrat=xStratfac) 

StrReport <- ftable(xtabs(w~obstype+yStrat+xStrat,data=df))
if(exists("df1")) xStrReport <-  xtabs(w~obstype+xStrat,data=df1)

if (xs.includes == TRUE) {
    oblevel0 <- levels(factor(as.character(obstype)))
    oblevel <- rep(oblevel0, rep(nyStrat,length(oblevel0)))  

    report <- ftable(xtabs(w~obstype+yStrat+xStrat, data=df[df$obstype=="strata",]))
    strep <- xtabs(w~yStrat+xStrat, data=df[df$obstype!="strata",])
    report[oblevel=="strata"] <- strep
    StrReport <- StrReport-report

    if (exists("df1")) {
        report1 <- xtabs(w~obstype+xStrat, data=df1[df1$obstype=="strata",])
        strep1 <- xtabs(w~xStrat, data=df1[df1$obstype!="strata",])
        report1[oblevel0=="strata"] <- strep1
        xStrReport <- xStrReport-report1
    }
}

if (!(is.null(xstrata))) {
    key <- dimnames(xstratmat)[[2]]
    names(key) <- 1:length(key)
}

if (!fit) {
    ans <- list(missReport=missReport, StrReport=StrReport, xStrReport=xStrReport, key=key, 
       	      fit=fit, call=call, assign1=assign1, assign2=assign2, assign3=assign3, 
              fnames1=fnames1, fnames2=fnames2, fnames3=fnames3, terms1=terms1, terms2=terms2,
              terms3=terms3, order1=order1, order2=order2, order3=order3, v1=length(terms1)+1,
              v2=length(terms2)+1,  v3=length(terms3)+1)
    class(ans) <- "bivlocsc2stg"
    return(ans) # RETURN AT THIS POINT IF fit IS FALSE
}


# ------------------------------------------------------------------------------------
rmat <- Qmat <- NULL

if (!prospective) {
   havestrata <- length(obstype[obstype=="strata"])>0
   rk <- rep(0,n)
   if ((havestrata & xs.includes) | !havestrata) 
      rk[obstype=="retro"] <- -1 * w[obstype=="retro"]
   if  (havestrata & xs.includes) 
      rk[obstype=="uncond"] <- -1 * w[obstype=="uncond"]
   rk[obstype=="strata"] <- 1 * w[obstype=="strata"]

   rmat <- xtabs(rk ~ yStratfac + xStratfac)
   rmat <- matrix(rmat,dim(rmat)[1],dim(rmat)[2],dimnames=dimnames(rmat))
   if (dim(rmat)[1] == 1){ # only one level of y1 represented in the data
      if (dimnames(rmat)[[1]]=="1") rmat <- rbind(rmat,rep(0,dim(rmat)[2]))
       	 else rmat <- rbind(rep(0,dim(rmat)[2]),rmat) }

   if (is.null(Qstart)){
     ind <- if (havestrata & !xs.includes) 
                 (1:n)[obstype=="uncond"|obstype=="retro"|obstype=="strata"]
            else (1:n)[obstype=="strata"]
     temp <- xtabs(w~yStrat+xStrat,data=data.frame(yStrat=yStratfac,xStrat=xStratfac,w)[ind,])

     if (is.matrix(temp) && (dim(temp) == c(2,nStrat)) && min(temp)>0)
       	    Qmat <- sweep(temp,2,apply(temp,2,sum),"/")
     else  stop("Attempt to construct missing Qstart failed")
   }
   else { # Qstart is present
      if(is.matrix(Qstart) && (dim(Qstart) == c(2,nStrat)) &&  # matrix with correct dimensions
         (round(apply(Qstart,2,sum),5) == rep(1,dim(Qstart)[2]))) # cols sum to 1
            Qmat <- Qstart
      else stop("illegal Qstart value")
   }
}


# Create start if needed ------------------------------------------------------------
if (is.null(start)){ # use glm() to obtain starting values
   ### calculate the weights (N/n)
   wmult <- rep(1,n) 
   ind <- 1:n
   if (!prospective) {
      if (all(!is.na(match(obstype,c("retro","xonly")))))  # all retro or xonly
 	          stop("Cannot construct starting values from this data")
   	
      ind <- (1:n)[obstype=="retro" | obstype=="y|x" | obstype=="uncond"]
      if (length(levels(factor(yStrat[ind])))<2 | sum(w[ind]) < 50)
     	        stop("Cannot construct starting values from this data")

      if (!is.na(match("retro",obstype))){ # adjust values of wmult if have retro obsns
           nMat1 <- xtabs(w~yStrat+xStrat,data=data.frame(yStrat=yStratfac,xStrat=xStratfac,
                        w,obstype)[obstype=="uncond" | obstype=="retro",])
           nMat <- xtabs(w~yStrat+xStrat,data=data.frame(yStrat=yStratfac,xStrat=xStratfac,
			w,obstype)[obstype=="retro",])
   	   NMat <- xtabs(w~yStrat+xStrat,data=data.frame(yStrat=yStratfac,xStrat=xStratfac,
			w,obstype)[obstype=="strata",])
   	   if (!xs.includes) NMat <- NMat + nMat1

     	   Nn <- (NMat+1)/(nMat+1)
           Nn <- sweep(Nn,2,apply(Nn,2,FUN=function(x)mean(x,na.rm=TRUE)),"/")
   	   Nnvec <- as.vector(Nn) 
   	   xystrat <- (xStrat-1)*nyStrat+yStrat
           wmult <- Nnvec[xystrat] # multiply by corresponding values in Nn
           if (any(is.na(wmult))) stop("Missing values in constructed wmult vector")
      }
   }
   
   start1 <- glm(y[,1][ind] ~ X1[ind,]-1, weights=w[ind]*wmult[ind],
	   	family=binomial)$coefficients

   gg <- lm(y[,2][ind] ~ X2[ind,]-1, weights=w[ind]*wmult[ind])
   start2 <- gg$coefficients

   scalef <- 1
   if (errdistrn=="logistic") scalef <- 1.8
   if (errdistrn=="t") {dfs <- errmodpars[1]; scalef <- sqrt(dfs/(dfs-2))}
	start3 <- c(log(summary(gg)$sigma/scalef),rep(0,ncol(X3)-1))
	
   start <- c(start1,start2,start3)
}
else { # Match up any values in start with names in Xs, fill in holes with 0's.
   newstart <- rep(0,dim(X1)[2]+dim(X2)[2]+dim(X3)[2])
   if (length(start) != length(newstart)) {
 	vname <- c(dimnames(X1)[[2]],dimnames(X2)[[2]],dimnames(X3)[[2]])
   	posns <- match(names(start),vname)
        newstart[posns] <- start
   	start <- newstart
   }
}


# Form xarray ----------------------------------------------------------------
hvalue <- y[obstype=="strata",1]
hxStrat <- xStrat[obstype=="strata"]
Cmult <- w[obstype=="strata"]

X1 <- X1[obstype != "strata",,drop=FALSE]
X2 <- X2[obstype != "strata",,drop=FALSE]
X3 <- X3[obstype != "strata",,drop=FALSE]
X4 <- X4[obstype != "strata",,drop=FALSE]
w <- w[obstype != "strata"]
y <- y[obstype != "strata",,drop=FALSE]
off.set <- off.set[obstype != "strata",,drop=FALSE]
xStrat <- xStrat[obstype != "strata"]
yStrat <- yStrat[obstype!="strata"]
obstype <- obstype[obstype != "strata"]
n <- length(w)

ntheta <- dim(X1)[2]+dim(X2)[2]+dim(X3)[2]
nthetaOrig <- c(dim(X1)[2], dim(X2)[2], dim(X3)[2])
thetaname <- c(dimnames(X1)[[2]],dimnames(X2)[[2]],dimnames(X3)[[2]])
nX <- 3

X <- array(0,c(n,ntheta,nX))
dimnames(X) <- list(NULL,thetaname,c("Y1mod","Y2mod1","Y2mod2"))
X[,1:dim(X1)[2],1] <- X1
X[,(dim(X1)[2]+1):(dim(X1)[2]+dim(X2)[2]),2] <- X2
X[,(dim(X1)[2]+dim(X2)[2]+1):ntheta,3] <- X3

Aposn <- (1:n)[obstype=="uncond"|obstype=="retro"|obstype=="y|x"]
Acounts <- w[Aposn]
if (!prospective) {
	Bposn <- (1:n)[obstype=="uncond"|obstype=="retro"|obstype=="xonly"]
	Bcounts <- w[Bposn]
}
else Bposn <- Bcounts <- numeric(0)


# ----------------------FUNCTION CALL----------------------------------------------
if (errdistrn=="logistic") errdist <- logisterr
else if (errdistrn=="normal") errdist <- normerr
else if (errdistrn=="t") errdist <- terr
else stop(paste("error distribution",errdistrn,"not implemented"))

modelfn <- lspml2locsc; stratfn <- NULL
ProspModInf <- MEtaProspModInf; StratModInf <- MEtaStratModInf.spml2locsc 

res <- mlefn(theta=start,MLInf,y=as.matrix(y),x=X,ProspModInf=ProspModInf,
	     StratModInf=StratModInf, modelfn=modelfn, stratfn=stratfn, errdist=errdist, 
             errmodpars=errmodpars, Aposn=Aposn,Acounts=Acounts,Bposn=Bposn,Bcounts=Bcounts,
 	     rmat=rmat,Qmat=Qmat,xStrat=as.numeric(xStrat), off.set=off.set, 
             X4=X4, y2posiny1=y2posiny1, nthetaOrig=nthetaOrig, strata=as.numeric(xStrat),
	     control=control, control.inner=control.inner)   

if (res$error != 1) {
  names(res$theta) <- names(res$score) <- dimnames(X)[[2]]
  dimnames(res$inf) <- list(dimnames(X)[[2]],dimnames(X)[[2]])
  #covmat <- solve(res$inf) # change to avoid scale problems with inverse
  temp <- 1/sqrt(diag(res$inf))
  covmat <- diag(temp)%*%solve(diag(temp)%*%res$inf%*%diag(temp))%*%diag(temp)
  dimnames(covmat) <- dimnames(res$inf)
  dd <- sqrt(diag(covmat))
  correlation <- covmat/outer(dd, dd)
  pred.Y2 <- X2Orig %*% res$theta[(dim(X1)[2]+1):(dim(X1)[2]+dim(X2)[2])]
}
else pred.Y2 <- covmat <- correlation <- NULL

# ------------------------------------------------------------------------------------
ans <- list(missReport=missReport, StrReport=StrReport, xStrReport=xStrReport, key=key, 
	    yKey=yKey, fit=fit, error=res$error, coefficients=res$theta, loglk=res$loglk, 
	    score=res$score, inf=res$inf, fitted.Y2=pred.Y2, cov=covmat, cor=correlation, 
	    Qmat = res$extra$Qmat, call=call, assign1=assign1, assign2=assign2, assign3=assign3, 
	    fnames1=fnames1, fnames2=fnames2, fnames3=fnames3, terms1=terms1, terms2=terms2, 
	    terms3=terms3,order1=order1, order2=order2, order3=order3, v1=length(terms1)+1,
            v2=length(terms2)+1,  v3=length(terms3)+1)
class(ans) <- "bivlocsc2stg"
ans
}


#####################################################################################
print.bivlocsc2stg <- function (x, digits = max(3, getOption("digits") - 3), ...) 
#####################################################################################    
{
cat("\nCall:\n", deparse(x$call), "\n", sep = "")

if (!is.null(x$missReport)){
	   cat("\nObservations deleted due to missing data:")
    print(x$missReport,quote = FALSE, row.names=FALSE)
}
cat("\nStratum Counts Report:\n")
print(x$StrReport,quote = FALSE, row.names=FALSE)
if (!is.null(x$xStrReport)) {
    cat("\nObservations of obstype==xonly\n")
    print(x$xStrReport)
}

cat("\nModels for prob of joint distribution of ", x$yKey[1]," and ",x$yKey[3],
	" given covariates\n",sep="")
if (!is.null(x$key)) {
    cat("\nKey to x-Strat:\n")
    print.default(x$key,quote=FALSE, row.names=FALSE)
}

if (x$fit) {
   if (x$error != 0) print("WARNING, FIT UNSUCCESSFUL")
    cat("\nCoefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nloglikelihood =",x$loglk, "  using", length(coef(x)), "parameters\n\n")

#   cat("Covariance:\n")
#   print.default(format(x$cov, digits = digits), print.gap = 2, quote = FALSE)
#   cat("\n")
}

invisible(x)
}


#####################################################################################
summary.bivlocsc2stg <- function(object) {
#####################################################################################

# Works for R only so far

z <- object
sigma <- z$cov
thetas <- z$coefficients    
assign1 <- z$assign1; assign2 <- z$assign2; assign3 <- z$assign3
terms1 <- z$terms1; terms2 <- z$terms2; terms3 <- z$terms3
order1 <- z$order1; order2 <- z$order2; order3 <- z$order3
fnames1 <- z$fnames1;  fnames2 <- z$fnames2;  fnames3 <- z$fnames3
 

# calculate wald statistics ---------------------------------------------------

dowald <- function(assign, terms, thetas, sigma) {
    # find dfs for each variable
    num <- -1
    pos <- waldchis <- pvals <- coef.table <- NULL

    for (i in 1:length(assign)) {
    if (assign[i] > num) { 
       pos <- c(pos,i) 
       num <- assign[i] 
       }
    }
    pos1 <- c(pos[-1],length(assign)+1)
    values <- assign[pos]    
    dfs <- lengths <- pos1 - pos

    # do wald calculations
    if(length(pos)>1){
    for (i in 2:length(pos)) {
        thetai <- thetas[pos[i]:(pos[i]+lengths[i]-1)]
        covi <- sigma[pos[i]:(pos[i]+lengths[i]-1),pos[i]:(pos[i]+lengths[i]-1)]
        waldchi <- t(thetai) %*% solve(covi,thetai)
        waldchis <- c(waldchis,waldchi)
        pval <- (1-pchisq(waldchi,df=dfs[i]))
        pvals <- c(pvals,pval) 
    }
    coef.table <- cbind(dfs[-1], waldchis, pvals)
    dimnames(coef.table)[[1]] <- terms[values]
    dimnames(coef.table)[[2]] <- c("Df","Chi", "Pr(>Chi)")
    }
    coef.table
}

l1 <- length(assign1); l2 <- length(assign2)
lt <- length(thetas)
ct1 <- dowald(assign1,terms1,thetas[1:l1],sigma[1:l1, 1:l1])
ct3 <- dowald(assign2,terms2,thetas[(l1+1):(l1+l2)], sigma[(l1+1):(l1+l2), (l1+1):(l1+l2)])
ct5 <- dowald(assign3,terms3,thetas[(l1+l2+1):lt], sigma[(l1+l2+1):lt, (l1+l2+1):lt])


 # calculate t statistics for non factors ---------------------------------------------
if(is.null(dim(sigma))) s.err <- sqrt(sigma) else s.err <- sqrt(diag(sigma))

zstats <- thetas/s.err
pvals <- (1 - pnorm(abs(zstats)))*2
coef.table2 <- cbind(thetas, s.err, zstats, pvals)
dimnames(coef.table2)[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

ct2 <- ct4 <- ct6 <- NULL

if (length(assign1) > 0) ct2 <- coef.table2[1:l1, ]
if (length(assign2) > 0) ct4 <- coef.table2[(l1+1):(l1+l2), ]
if (length(assign3) > 0) ct6 <- coef.table2[(l1+l2+1):lt, ]

ans <- list(missReport=z$missReport, call=z$call,StrReport=z$StrReport, xStrReport=z$xStrReport,
       key=z$key, yKey=z$yKey, fit=z$fit, error=z$error, loglk=z$loglk, coefficients=z$coefficients,
       coef.table1=ct1, coef.table2=ct2, coef.table3=ct3, coef.table4=ct4, coef.table5=ct5,
       coef.table6=ct6)
class(ans) <- "summary.bivlocsc2stg"
ans
}    


#####################################################################################
print.summary.bivlocsc2stg <- function (x, digits = max(3, getOption("digits") - 3), ...) {
#####################################################################################    
cat("\nCall:\n", deparse(x$call), "\n", sep = "")

if (!is.null(x$missReport)){
	   cat("\nObservations deleted due to missing data:")
    print(x$missReport,quote = FALSE, row.names=FALSE)
}
cat("\nStratum Counts Report:\n")
print(x$StrReport,quote = FALSE, row.names=FALSE)
if (!is.null(x$xStrReport)) {
    cat("\nObservations of obstype==xonly\n")
    print(x$xStrReport)
}

cat("\nModels for prob of joint distribution of ", x$yKey[1]," and ",x$yKey[3],
	" given covariates\n",sep="")
if (!is.null(x$key)) {
    cat("\nKey to x-Strat:\n")
    print.default(x$key,quote=FALSE, row.names=FALSE)
}

if (x$fit) {
   if (x$error != 0) print("WARNING, FIT UNSUCCESSFUL")
    cat("\nloglikelihood =",x$loglk, "  using", length(coef(x)), "parameters\n\n")  
  
    if(!is.null(x$coef.table1)){
        cat("\nY1-Model\n") 
        cat("Wald Tests:\n")
        print(x$coef.table1,digits=4)
        cat("\n")
    }
    if(!is.null(x$coef.table2)){
        cat("Coefficients:\n")
        print(x$coef.table2,digits=4)
    }   

    if(!is.null(x$coef.table3)){
        cat("\nY2-Model (Location)\n") 
        cat("Wald Tests:\n")
        print(x$coef.table3,digits=4)
        cat("\n")
    }   
    if(!is.null(x$coef.table4)){
        cat("Coefficients:\n")
        print(x$coef.table4,digits=4)
    }   

    if(!is.null(x$coef.table5) | !is.null(x$coef.table6))
        cat("\nY2-Model (Scale)\n")
    if(!is.null(x$coef.table5)) {
        cat("Wald Tests:\n")
        print(x$coef.table5,digits=4)
        cat("\n")
    }   
    if(!is.null(x$coef.table6)){
        cat("Coefficients:\n")
        print(x$coef.table6,digits=4) 
    }   
}
else cat("\nCall requested that model not be fitted.\n")

invisible(x)
}

