test_that("phi2w and w2phi are inverse functions", {
	ws <- random_weights(15)
	expect_equal(phi2w(w2phi(ws)), ws)
	phis <- runif(20, pi/6, pi/4)
	expect_equal(w2phi(phi2w(phis)), phis)
})