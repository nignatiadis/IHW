test_that("number of elements in each group should differ by at most 1", {
    expect_uniform_groups <- function(groups){
        expect_true(diff(range(table(groups))) <= 1)
    }
    set.seed(1)
    x1 <- runif(1000)
	expect_uniform_groups(groups_by_filter(x1,7))

    x2 <- c(rep(0,100), runif(100))
    # b2 <- cut(x2, quantile(x, probs=seq(0, 1, length.out=8+1))) #error
    expect_uniform_groups(groups_by_filter(x2, 8))

    x3 <- c(rep(0,10), runif(100))
    # b3 <- cut(x3, quantile(x, probs=seq(0, 1, length.out=8+1))) #NAs introduced
    expect_uniform_groups(groups_by_filter(x3, 8))


})