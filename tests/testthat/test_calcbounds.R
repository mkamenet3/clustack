########################################
#Tests for Functions
########################################
test_that("calcbounds switches from cluster RR to cell RR", {
    lik <- c(0.5, 0.75, 0.3, 0.4, 0.01)
    expect_that(round(likweights(lik),3), 
                equals(c(0.188, 0.166, 0.208, 0.198, 0.240)))
    
})


calcbounds <- function(cellrates=FALSE, cellsix=NULL){
    if(is.null(cellsix)){
        cellsix <- TRUE
    } 
    print(cellrates)
    print(cellsix)
}

calcbounds()
calcbounds(cellrates = TRUE, cellsix = c(1,2,3))




