test_that("number of elements in each group should differ by at most 1", {
    check_uniform_groups <- function(groups){
        diff(range(table(groups))) <= 1
    }

 	groups <- groups_by_filter(runif(1000),7)
	expect_true(check_uniform_groups(groups))
})
