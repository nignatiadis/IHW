# small suite of functions to check if input from user is correct etc.
# also useful for unit tests!

check_uniform_groups <- function(groups){
	diff(range(table(groups))) <= 1
}
