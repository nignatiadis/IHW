# below code is copied from deprecated CRAN 'pi0' package
# https://github.com/gitlongor/pi0/blob/master/R/CBUM.R
# Long Qu modified the code from Markitsis and Lai (2010)
# This file is licensed under GPL (>= 2.0)

CBUM <- function(p, start.pi0=0.5, thresh.censor=0.05, eps=1e-5, niter=Inf, verbose=FALSE) {
    verbose=isTRUE(verbose)
    a <- NA
    G <- length(p)
    pi1 <- 1
    if(thresh.censor<=0) thresh.censor=min(p)/2
    stopifnot(thresh.censor>0)
    #indices for which the pi are <thresh.censor 
    idx.uncencor <- which(p>=thresh.censor)
    n.censor=G-length(idx.uncencor)
    sum.z.censor=(1-start.pi0)*n.censor
    z.uncensor=rep(1-start.pi0, G-n.censor)
    sum.z.uncensor=sum(z.uncensor)

    log.p.uncensor=log(p[idx.uncencor])
    log.thresh.censor=log(thresh.censor)

    iter=1
    repeat {
        pi_temp <- (sum.z.censor + sum(z.uncensor)) / G
        a_temp <- - sum.z.uncensor / ( log.thresh.censor*sum.z.censor + sum(log.p.uncensor * z.uncensor) )

        z.censor.num = pi_temp * thresh.censor^(a_temp)
        sum.z.censor <- n.censor * ( z.censor.num / ( (1-pi_temp) * thresh.censor + z.censor.num ) )

        z.uncensor.num = pi_temp * a_temp * exp(log.p.uncensor * (a_temp-1)) 
        z.uncensor <-  z.uncensor.num / ( 1-pi_temp + z.uncensor.num) 
        sum.z.uncensor = sum(z.uncensor)

        if(verbose)
            cat('iter',iter,'\tgamma=',1-pi_temp, '\talpha=',a_temp, '\tmax.diff=',max(abs(pi1-pi_temp), abs(a-a_temp)), fill=TRUE)
        if( iter>=niter || (abs(pi1-pi_temp) < eps && isTRUE(abs(a-a_temp)<eps)) ) break
        iter=iter+1
        pi1 <- pi_temp
        a <- a_temp

    }

    ans = ( (1-pi1)+pi1*a )
    attr(ans, 'converged')=iter<niter
    attr(ans, 'iter')=iter
    attr(ans, 'alpha')=a_temp
    attr(ans, 'lfdr')=if(thresh.censor<min(p)) ans/((1-pi1)+pi1*a*p^(a-1)) else NULL
    attr(ans, 'thresh.censor')=thresh.censor
    attr(ans, 'call')=match.call()
    class(ans) = 'CBUM'
    ans
}