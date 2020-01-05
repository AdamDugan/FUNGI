variant<-function(x)
{
  length(unique(x))>1
}

flip.check <- function(v1, v2){
    NA.cols <- 0
    ## check if cols have "NA" values
    c1 <- tryCatch(as.numeric(v1),error=function(e) e, warning=function(w) w)
    c2 <- tryCatch(as.numeric(v2),error=function(e) e, warning=function(w) w)
    if(is(c1,"warning")) NA.cols <- NA.cols+1
    if(is(c2,"warning")) NA.cols <- NA.cols+1
        
    ## make sure both vectors are numeric
    v1 <- suppressWarnings(as.numeric(v1))
    v2 <- suppressWarnings(as.numeric(v2))
    
    dif = v1-v2
    flip.pos = which(abs(dif)==2)
    # if there are flips, do 
    if(length(flip.pos)!=0){
        zero.pos = which(v2==0)
        two.pos = which(v2==2)
        zero.noflip = setdiff(zero.pos,intersect(flip.pos,zero.pos)) #v1=0 or v1=1; v2=0 SNPs
        two.noflip = setdiff(two.pos,intersect(flip.pos,two.pos))    #v1 =2 or v1 = 1; v2 = 2
        prev.zero = which(v1==0)        # Sneha says v2==0
        zerozero = intersect(prev.zero,zero.noflip) #v1=0;v2=0
        prev.two = which(v2==2)
        twotwo = intersect(prev.two,two.noflip)
        no.flip = length(zerozero)+length(twotwo)
        # if re-coding does not create more flips in other positions, do
        if(length(flip.pos)>no.flip){
            v2[flip.pos] <- v2[flip.pos]+2*sign(dif[flip.pos])
            v2[zero.noflip] <- v2[zero.noflip]+2
            v2[two.noflip] <- v2[two.noflip]-2
        }
        else{
            v2 <- v2
        }
    }
    # if no flips, do nothing
    else{
        v2 <- v2
    }
    result <- list(v2, NA.cols)
    return(result)
}

trait.simu2<-function(geno,func.frq,mean,sd,miu,
                     cov=NULL, dist=c("normal","cauchy","poisson","binary"),Ntrait=2)
{
  dist <- match.arg(dist);
  
  idx.true<-sample(1:ncol(geno),ncol(geno)*func.frq);
  geno<-as.matrix(geno[,idx.true]);
  
  #beta<-rnorm(ncol(geno),mean,sd)*abs(log10(maf));
  beta<-rnorm(ncol(geno),mean,sd);
  #beta<-runif(ncol(geno),mean-sd,mean+sd)
  
  if(is.null(cov)){
    y.tmp<-geno%*%beta
  }else{
    alpha<-rep(c(-1,1),ncol(cov)/2);
    y.tmp<-geno%*%beta+cov%*%alpha;
  }
  
  Y<-NULL;
  
  for(i in 1:Ntrait){
    if(dist=="normal"){
      y<-miu+y.tmp+rnorm(nrow(geno),0,1)
    }else if(dist=="cauchy"){
      y<-miu+y.tmp;
      y<-rcauchy(length(y),y,)
    }else if(dist=="poisson"){
      y<-exp(-miu-y.tmp);
      y<-rzipois(length(y), y, pstr0 =0.6)
      #y<-rpois(length(y), y)
    }else if(dist=="binary"){
      y<-1/(1+exp(-miu-y.tmp));
      y<-as.numeric(y>runif(nrow(geno)))
    }
    Y<-cbind(Y,y)
  }
  Y
}

trait.simu<-function(geno,func.frq,mean,sd,miu,
                     cov=NULL, dist=c("normal","cauchy","poisson","binary"))
{
  dist <- match.arg(dist);
  
  idx.true<-sample(1:ncol(geno),ncol(geno)*func.frq);
  geno<-as.matrix(geno[,idx.true]);
  
  #beta<-rnorm(ncol(geno),mean,sd)*abs(log10(maf));
  beta<-rnorm(ncol(geno),mean,sd);
  #beta<-runif(ncol(geno),mean-sd,mean+sd)
  
  if(is.null(cov)){
    y.tmp<-geno%*%beta
  }else{
    alpha<-rep(c(-1,1),ncol(cov)/2);
    y.tmp<-geno%*%beta+cov%*%alpha;
  }
    
  if(dist=="normal"){
    y<-miu+y.tmp+rnorm(nrow(geno),0,1)
  }else if(dist=="cauchy"){
    y<-miu+y.tmp;
    y<-rcauchy(length(y),y,)
  }else if(dist=="poisson"){
    y<-exp(-miu-y.tmp);
    y<-rzipois(length(y), y, pstr0 =0.6)
    #y<-rpois(length(y), y)
  }else if(dist=="binary"){
    y<-1/(1+exp(-miu-y.tmp));
    y<-as.numeric(y>runif(nrow(geno)))
  }
  as.vector(y)
}


## Using truncated power basis
smooth.construct.tr.smooth.spec<-function(object,data,knots)
## a truncated power spline constructor method function
## object$p.order = null space dimension
{ m <- object$p.order[1]
  if (is.na(m)) m <- 2 ## default
  if (m<1) stop("silly m supplied")
  if (object$bs.dim<0) object$bs.dim <- 10 ## dim(H1) or number of terms to be penalized
  nk<-object$bs.dim-m-1 ## number of knots
  if (nk<=0) stop("k too small for m")
  x <- data[[object$term]] ## the data
  x.shift <- mean(x) # shift used to enhance stability
  k <- knots[[object$term]] ## will be NULL if none supplied
  if (is.null(k)) # space knots through data
      { n<-length(x)
        k<-quantile(x[2:(n-1)],seq(0,1,length=nk+2))[2:(nk+1)]
    }
  if (length(k)!=nk) # right number of knots?
      stop(paste("there should be ",nk," supplied knots"))
  x <- x - x.shift # basis stabilizing shift
  k <- k - x.shift # knots treated the same!
  X<-matrix(0,length(x),object$bs.dim)
  for (i in 1:(m+1)) X[,i] <- x^(i-1)
  for (i in 1:nk) X[,i+m+1]<-(x-k[i])^m*as.numeric(x>k[i])
  object$X<-X # the finished model matrix
  if (!object$fixed) # create the penalty matrix
      { object$S[[1]]<-diag(c(rep(0,m+1),rep(1,nk)))
    }
  object$rank<-nk # penalty rank
  object$null.space.dim <- m+1 # dim. of unpenalized space
  ## store "tr" specific stuff ...
  object$knots<-k;object$m<-m;object$x.shift <- x.shift
  object$df<-ncol(object$X) # maximum DoF (if unconstrained)
  class(object)<-"tr.smooth" # Give object a class
  object
}

library(mgcv)
library(nlme)
##library(vegan)
##library(refund)

## check and flip genotype codeing to minimize genotype value
## turbulance so fitted curves are smoother.
flip.xtong <- function(v1, v2, level = 2L)
{
### v1 ---- the reference variant to be compared against
### v2 ---- the target variant to be check and flipped.
### NA is tolarated by ignoring these indivudials' contribution
### to terbulence.
### level ---- the top genotype value, 2 for dosage genotype
### coded as 0, 1, and 2.

    ## non-NA mask
    msk <- which(!is.na(v1+v2));

    ## compliment of v2 <=> flipped v2
    v2.c <- level - v2;

    ## compare the before and after flipping terbulence between
    ## two variants across all individual, if NA turns up at an
    ## individual of either varaint, this individual contribute
    ## no turbulence.
    ## d is in monotune with sum(|v1-v2|-|v1-v2.c|)
    d <- crossprod(level-v1[msk]*2L, level-v2.c[msk]*2L);

    if(d > 0L)
    {
        return (v2.c);
    }
    else
    {
        return (v2);
    }
}

## flip genotype codeing to reduce turbulance
## for a group of individuals
## gmx ---- the genotype matrix
## mrg ---- which dimension of the matrix represents
## the variants(so the other is for idividuals)?
##    1 --> row variant, col individual
##    2 --> col variant, row individual
gflp.xtong <- function(gmx, mrg=1L)
{
    ## allocate the output - flipped genotype matrix;
    out <- matrix(data=0L, nrow=nrow(gmx), ncol=ncol(gmx));

    if(mrg==1L) # row variant, column individual
    {
        ## the first variant doesn't require flipping
        out[1L,] <- gmx[1L,];

        ## deal with the rest, use previous(i-1 th.) variant as a
        ## reference, check and flip current variant if necessary
        for(i in 2L:nrow(gmx))
        {
            out[i,] <- flip.xtong(out[i-1L,], gmx[i,]);
        }
    }
    else        # row individual, column variant.
    {
        ## the first variant doesn't require flipping
        out[,1L] <- gmx[,1L];

        ## deal with the rest, use previous(i-1 th.) variant as a
        ## reference, check and flip current variant if necessary
        for(i in 2L:ncol(gmx))
        {
            out[,i] <- flip.xtong(out[,i-1L], gmx[,i]);
        }
    }
    out;
}

## Fit genotype function so the genotypes can be smoothed.
## gmx ---- genotype matrix.
## pos ---- genotype position.
## mrg ---- margin of variant in the matrix.
##     1: row major -- row variant, col individual
##     2: col major -- col variant, row individual
## return smoothed genotype matrix.
gfit.apply <- function(gmx, pos, mrg=1L)
{
    ## integrety check
    ndv <- dim(gmx)[3L-mrg];     # number of individuals
    ngv <- dim(gmx)[mrg];        # number of variants
    stopifnot(length(pos)==ngv);

    ## standarize variant position to [0,1]
    pos = (pos - pos[1L])/(pos[ngv] - pos[1L]);

    ## fit function from original genotype matrix
    fmx<-apply(X = gmx, MARGIN = 3L-mrg, FUN=function(idv)
    {
        ## fit genotype function for the j th. individual
        f <- try(gam(idv ~ s(pos, fx = FALSE, k = -1L, bs = 'cr')));
        if (inherits(f, "try-error"))
        {
            cat("Subject", j ,"has too many NA values to do anything meaningful",'\n')
            rep_int(x = NA, times = ngv);
        }
        else
        {
            f$fitted.values;
        }
    });
    if(mrg==2L)
    {
        fmx <- t(fmx);
    }
    fmx;
}
