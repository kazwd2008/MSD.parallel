################################################################################
#  Modified Stahel-Donoho Estimators for Multivariate Outlier Detection  
#      Ver.3.1 2013/06/15
#      Ver.3.2 2021/11/15  Outlier detection step added
################################################################################
#   by WADA, Kazumi 
################################################################################

RMSD.p <- function(inp, nb=0, sd=0, pt=0.999, dv=10000) {
 inp.d <- ncol(inp)        			# number of variables 
 inp.n <- nrow(inp)        			# number of observations 

#--------------------------------
#    parallelization 
#--------------------------------
require(doParallel)
require(foreach)
  type <- if (exists("mcfork", mode="function")) "FORK" else "PSOCK"
    cores <- getOption("mc.cores", detectCores())
    cl <- makeCluster(cores, type=type)    
    registerDoParallel(cl)
  RNGkind("L'Ecuyer-CMRG")     
    if (sd != 0) set.seed(sd)   

#--------------------------------
#     create orthogonal bases
#--------------------------------
## set number of bases so that it can be devided by number of cores
## "dv" is the max number of elements in a chunk

if (nb == 0) bb.n <- trunc(exp(2.1328+0.8023*inp.d) / inp.d) else bb.n <- nb
dv.cr <- ceiling(dv/cores/(inp.d^2))  #  Number of bases in a chunk
bb.cr <- ceiling(bb.n / dv.cr)	 #  Number of chunks
rn.cr <- dv.cr * inp.d^2         #  Number of elements which consists of bases in a chunk
kijun 	<- qchisq(0.95, inp.d)   #  reference for trimming
clusterExport(cl, "orthonormalization")

#-----------------------------------------------------------
#  projection, residual and weights computation in parallel
#-----------------------------------------------------------
bwt.cr <- foreach(cr=1:bb.cr, .combine='c') %dopar% { 
  res 	<- array(0, c(inp.n, inp.d, dv.cr))     #  residuals
  wt 	<- array(0, c(inp.n, inp.d, dv.cr))     #  weights	

  Fprj	<- function(pj) t(pj %*% t(inp))	# projection

  basis <- array(runif(rn.cr), c(inp.d, inp.d, dv.cr)) 	
  #basis  <- apply(basis, 3, gso)
  basis   <- apply(basis, 3, orthonormalization)
    basis   <- array(basis, c(inp.d, inp.d, dv.cr)) 

  prj <- apply(basis, 3, Fprj)	
    prj <- array(prj, c(inp.n, inp.d, dv.cr))

  medi <- apply(prj, c(2, 3), median)        # median
  madx <- apply(prj, c(2, 3), mad)      # MAD * 1.4826 (MAD / 0.674)

  for (l in 1:dv.cr) {          #  robust standardization of residuals
      res[,,l] <- t(abs(t(prj[,,l]) - medi[,l]) / madx[,l])
  }      

  # trimming weight
  k0	   <- which(res <= sqrt(kijun))	
  k1	   <- which(res > sqrt(kijun))	
  wt[k0]  <- 1				
  wt[k1]  <- kijun / (res[k1]^2)	

  wts <- apply(wt, c(1,3), prod)
  apply(wts, 1, min)     #  selecting the smallest weight for each observation
}      

stopCluster(cl)
#-----------------------------------------------------------

bwt.cr <- array(bwt.cr, c(inp.n, cores))  
bwt    <- apply(bwt.cr, 1, min)    #  selecting the smallest weight through chunks

### initial robust covariance matrix
u1 <- apply(inp * bwt, 2, sum) / sum(bwt)
V1 <- t(t(t(inp) - u1) * bwt) %*% (t(t(inp) - u1) * bwt) / sum(bwt^2)	

### avoiding NaN error
u1 <- ifelse(is.nan(u1), 0, u1)
V1 <- ifelse(is.nan(V1), 0, V1)

### robust PCA (LAPACK)
eg	<- eigen(V1, symmetric=TRUE)
ctb	<- eg$value / sum(eg$value) 	# contribution ratio

##############################
# projection pursuit (PP)
##############################

res2	<- array(0, c(inp.n, inp.d))     # residuals
wt2	<- array(0, c(inp.n, inp.d))         # weight by observations x variables
wts2 	<- array(0, inp.n)               # final weight by observations

prj2  <- t(eg$vector %*% (t(inp) - u1))  # projection
medi2 <- apply(prj2, 2, median)          # median and 
madx2 <- apply(prj2, 2, mad)             #  MAD for standadization
res2 <- t(abs(t(prj2) - medi2) / madx2)  # standardized residuals

### trimming
k0	   <- which(res2 <= sqrt(kijun))
k1	   <- which(res2 > sqrt(kijun))	
wt2[k0]  <- 1			
wt2[k1]  <- kijun / (res2[k1]^2)	
wts2 <- apply(wt2, 1, prod)  
wts2 <- pmin(wts2, bwt) 

##############################
# final mean vector and covariance matrix
##############################
 
u2 <- apply(inp * wts2, 2, sum) / sum(wts2)			
V2 <- t(t(t(inp) - u2) * wts2) %*% (t(t(inp) - u2) * wts2) / sum(wts2^2)

#return(list(u1=u1, V1=V1, bwt=bwt, u2=u2, V2=V2, wts2=wts2, eg=eg, ctb=ctb))

##############################
# outlier detection 
##############################

# compute mahalanobis distance of each observation for outlier detection
mah <- mahalanobis(inp, u2, V2)
FF <- mah * ( inp.n - inp.d )* inp.n /(( inp.n^2 - 1 )* inp.d)	# test statistics
cf <- qf(pt, inp.d, inp.n - inp.d) 	# 99.9 percentile point of F distribution
ot <- rep(1, inp.n)			# outlier flag (1: OK;   2: outlier)
ot[which(FF > cf)] <- 2 		# outliers are those exceed the threshold

return(list(u=u2, V=V2, wt=wts2, mah=mah,fs=cf, ot=ot))
}

#################################################################################
# orthonormalization: Gram-Schmidt Orthonormalization contained in "far" package 
#################################################################################
# A set of unit vectors is returned in case of collinearlity.

orthonormalization <- function (u = NULL, basis = TRUE, norm = TRUE) {
#    if (is.null(u)) 
#        return(NULL)
#    if (!(is.matrix(u))) 
#        u <- as.matrix(u)
    p <- nrow(u)
    n <- ncol(u)
#    if (prod(abs(La.svd(u)$d) > 1e-08) == 0) 
#        stop("colinears vectors")
#    if (p < n) {
#        warning("too much vectors to orthogonalize.")
#        u <- as.matrix(u[, 1:p])
#        n <- p
#    }
#    if (basis & (p > n)) {
#        base <- diag(p)
#        coef.proj <- crossprod(u, base)/diag(crossprod(u))
#        base2 <- base - u %*% matrix(coef.proj, nrow = n, ncol = p)
#        norm.base2 <- diag(crossprod(base2))
#        base <- as.matrix(base[, order(norm.base2) > n])
#        u <- cbind(u, base)
#        n <- p
#    }
    if (prod(abs(La.svd(u)$d) > 1e-08) == 0) {		# added by wada
        warning("collinears vectors")			# added by wada
        v <- matrix(0, nr=p, nc=p)			# added by wada
        diag(v) <- 1					# added by wada
        return(v)					# added by wada
    }    						# added by wada
    v <- u
    if (n > 1) {
        for (i in 2:n) {
            coef.proj <- c(crossprod(u[, i], v[, 1:(i - 1)]))/diag(crossprod(v[, 
                1:(i - 1)]))
            v[, i] <- u[, i] - matrix(v[, 1:(i - 1)], nrow = p) %*% 
                matrix(coef.proj, nrow = i - 1)
        }
    }
    if (norm) {
        coef.proj <- 1/sqrt(diag(crossprod(v)))
        v <- t(t(v) * coef.proj)
    }
    return(v)
}
###################################################################