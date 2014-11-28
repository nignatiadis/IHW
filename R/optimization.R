# functions that help with the optimization procedures along the n-dimensional simplex!



# use subplex as default for now
simplex_optimization <- function(npars, fn){
	par <- w2phi(rep(1,npars))
	opt  <-  subplex(
        par = par,
        fn = fn,
        control = list(trace=0))
	phi2w(opt$par)
}

