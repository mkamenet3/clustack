########################################
#Tests for Functions
########################################
test_that("Likelihood weights correct", {
    lik <- c(0.5, 0.75, 0.3, 0.4, 0.01)
    expect_that(round(likweights(lik),3), 
                equals(c(0.188, 0.166, 0.208, 0.198, 0.240)))
    
})