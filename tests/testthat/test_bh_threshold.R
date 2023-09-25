test_that("applying BH via p.adjust or get_bh_threshold should yield the same result",{
	set.seed(1)
	pv <- c(runif(10000),rbeta(2000,.5,7))
	expect_equal(sum(p.adjust(pv, method="BH") <= 0.01), sum(pv <= get_bh_threshold(pv, 0.01)))
})
