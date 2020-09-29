########################################
#Tests for Functions
########################################
#Buckland toy example: change
foo <- function(a,b, cell=FALSE){
    if(cell==TRUE){
        print("works")
    } else {
        print("nonlog")
    }
}
foo(2,4,cell=TRUE)
foo(2,4, cell=FALSE)
foo(2,4)