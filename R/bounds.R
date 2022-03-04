
# bucklandbounds.cells2 <- function(thetaa, res, w, outObs,outExp,IC,transform="none",tsparsemat,overdisp.est, cellsix_out, conf.level) {
#     thetai <- res$Lambda_dense
#     critval <- qnorm(1-(1-conf.level)/2)
#     #print(critval)
#     if(IC=="aic"){
#         thetaa <- res$wLambda[res$selection.aic,][cellsix_out]
#     } else {
#         thetaa <- res$wLambda[res$selection.bic,][cellsix_out]
#     }
#     #print(str(thetaa))
#     
#
#THIS WORKS2!! 
test0 <- t(sparsemat)%*%Yx #these are all the total counts in each PC. These then need to be expanded out to
#each cell, so each cell will have the same number
test1a <- hadamard.prod(as.matrix(t(sparsemat)), matrix(test0@x)%*%matrix(data = 1, nrow = 66870, ncol = 1))
#sparsemat*matrix(as.vector(test0), ncol=1)
test1a <- as.matrix(sweep(sparsemat, MARGIN=2, test0@x, `*`))
test1a[test1a==0] <- Inf
test2a<- 1/test1a

test3 <- (log(thetai[,cellsix[5]]) - log(thetaa[5]))^2
test3a <- (sqrt(test2a + test3))%*%w
test3a[cellsix_out[5]]
LBa = exp(log(thetaa[5])-critval*(test3a[cellsix_out[5]]))
UBa = exp(log(thetaa[5])+critval*(test3a[cellsix_out[5]]))


# ################################################################
# 
# #THIS WORKS2!! Need to figure out the outside cluster bounds, should be 0
# test0 <- t(sparsemat)%*%Yx #these are all the total counts in each PC. These then need to be expanded out to
# #each cell, so each cell will have the same number
# test1a <- hadamard.prod(as.matrix(t(sparsemat)), matrix(test0@x)%*%matrix(data = 1, nrow = 66870, ncol = 1))
#     #sparsemat*matrix(as.vector(test0), ncol=1)
# test1a <- as.matrix(sweep(sparsemat, MARGIN=2, test0@x, `*`))
# test1a[test1a==0] <- Inf
# test2a<- 1/test1a
# 
# test3 <- (log(thetai[,cellsix[3]]) - log(thetaa[3]))^2
# test3a <- (sqrt(test2a + test3))%*%w
# test3a[cellsix_out[5]]
# LBa = exp(log(thetaa[5])-critval*(test3a[cellsix_out[5]]))
# UBa = exp(log(thetaa[5])+critval*(test3a[cellsix_out[5]]))
# 
# 
# 
# ################################################################    
# #    #THIS WORKS-ISH, need to figure out zeros
#     if(transform=="log")#{
#         #print("transform")
#         if(!is.null(overdisp.est)){
#             varthetai_1 <- as.matrix(as.vector(outYx)*sparsemat)
#             varthetai_1[varthetai_1==0] <- Inf
#             varthetai <- overdisp.est*(1/varthetai_1)
#         } else {
#             test0 <- t(sparsemat)%*%Yx
#             test1 <- colSums(as.matrix(as.vector(test0)*sparsemat))
#             test1[test1==0] <- Inf
#             varthetai<- 1/test1
#             
#             
#             
#             
#             varthetai_1 <- as.matrix(as.vector(outYx)*sparsemat)
#             varthetai_1[varthetai_1==0] <- Inf
#             varthetai <- 1/varthetai_1
#         }
#         test3 <- (log(thetai[,cellsix[3]]) - log(thetaa[3]))^2
#         test4 <- sqrt(varthetai + test3)%*%matrix(w, ncol = 1)
#         LBa = exp(log(thetaa[3])-critval*(test4))
#         UBa = exp(log(thetaa[3])+critval*(test4))
#         
#         
#         #bb <- sapply(1:nrow(res$Lambda_dense), function(k) log(res$Lambda_dense[k,]) - log(round(res$wLambda[1,],3)))
#         se_thetaa <- as.vector(sqrt(varthetai + bb^2)%*%matrix(w, ncol = 1))
#         LBa = exp(round(log(res$wLambda[1,],3)) -critval*(se_thetaa))
#         UBa = exp(round(log(res$wLambda[1,],3)) +critval*(se_thetaa))
# 
# #         
# #     } else {
# #         #browser()
# #         if(!is.null(overdisp.est)){
# #             varthetai_1 <- as.matrix(as.vector(outExp)*sparsemat)
# #             varthetai_1[varthetai_1==0] <- Inf
# #             varthetai <- overdisp.est*(thetai/t(varthetai_1))
# #         } else {
# #             varthetai <- (thetai/t(varthetai_1))#sapply(1:nrow(tsparsemat), function(k) thetai[k,]/outExp_out )
# #         }
# #         withintheta <- sapply(1:length(cellsix_out), function(j) (thetai[,cellsix_out][,j] - thetaa[j]))^2
# #         wtnbtn <- sapply(1:nrow(tsparsemat), function(k) sqrt(varthetai[k,cellsix_out] + withintheta[k,]))
# #         varthetas_w <- matrix(w, nrow = 1)%*%t(wtnbtn)
# #         var_thetaa <- (as.vector(varthetas_w))^2
# #         UBa = as.vector(thetaa) + critval*sqrt(var_thetaa)
# #         LBa = as.vector(thetaa) - critval*sqrt(var_thetaa)
# #     }
# #     return(list(buckland.LB = LBa,
# #                 clusterstack= thetaa,
# #                 buckland.UB = UBa))
# # }



#'@title bucklandbounds.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on Buckland et al. for specific cells (given by \code{cellsix_out}) 
#'@param thetaa Stacked relative risk estimate for the cluster(s)
#'@param res Resultant object from \code{detectclusters()}.
#'@param w Likelihood-weights for potential clusters.
#'@param outObs Observed counts for each potential cluster.
#'@param outExp Expected counts for each potential cluster.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC).
#'@param transform If a transformation is used. Default is \code{"none"}. Other transformation currently available is \code{"log"}.
#'@param tsparsemat Transpose of large sparse matrix where columns are potential clusters and rows are space-time locations.
#'@param overdisp.est Estimate of overdispersion (or underdispersion).
#'@param cellsix_out Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate for each of the cells in \code{cellsix_out}
#'@export
#'
# bucklandbounds.cells <- function(thetaa, res, w, outObs,outExp,IC,transform="none",tsparsemat,overdisp.est, cellsix_out, conf.level) {
#     thetai <- res$Lambda_dense
#     critval <- qnorm(1-(1-conf.level)/2)
#     #print(critval)
#     if(IC=="aic"){
#         thetaa <- res$wLambda[res$selection.aic,][cellsix_out]
#     } else {
#         thetaa <- res$wLambda[res$selection.bic,][cellsix_out]
#     }
#     #print(str(thetaa))
#     if(transform=="log"){
#         #print("transform")
#         if(!is.null(overdisp.est)){
#             varthetai_1 <- as.matrix(as.vector(outYx)*sparsemat)
#             varthetai_1[varthetai_1==0] <- Inf
#             varthetai <- overdisp.est*(1/varthetai_1)
#         } else {
#             varthetai_1 <- as.matrix(as.vector(outYx)*sparsemat)
#             varthetai_1[varthetai_1==0] <- Inf
#             varthetai <- 1/varthetai_1
#         }
#         withintheta <- sapply(1:length(cellsix_out), function(j) (log(thetai[,cellsix_out][,j]) - log(thetaa[j])))^2
#         wtnbtn <- sapply(1:length(cellsix_out), function(k) sqrt(varthetai[k] + withintheta[,k]))
#         varthetas_w <- matrix(w, nrow = 1)%*%wtnbtn
#         var_thetaa <- (as.vector(varthetas_w))^2
#         UBa = exp(as.vector(log(thetaa)) + critval*sqrt(var_thetaa))
#         LBa = exp(as.vector(log(thetaa)) - critval*sqrt(var_thetaa))
#         
#     } else {
#         #browser()
#         if(!is.null(overdisp.est)){
#             varthetai_1 <- as.matrix(as.vector(outExp)*sparsemat)
#             varthetai_1[varthetai_1==0] <- Inf
#             varthetai <- overdisp.est*(thetai/t(varthetai_1))
#         } else {
#             varthetai <- (thetai/t(varthetai_1))#sapply(1:nrow(tsparsemat), function(k) thetai[k,]/outExp_out )
#         }
#         withintheta <- sapply(1:length(cellsix_out), function(j) (thetai[,cellsix_out][,j] - thetaa[j]))^2
#         wtnbtn <- sapply(1:nrow(tsparsemat), function(k) sqrt(varthetai[k,cellsix_out] + withintheta[k,]))
#         varthetas_w <- matrix(w, nrow = 1)%*%t(wtnbtn)
#         var_thetaa <- (as.vector(varthetas_w))^2
#         UBa = as.vector(thetaa) + critval*sqrt(var_thetaa)
#         LBa = as.vector(thetaa) - critval*sqrt(var_thetaa)
#     }
#     return(list(buckland.LB = LBa,
#                 clusterstack= thetaa,
#                 buckland.UB = UBa))
# }

bucklandbounds.cells <- function(thetaa, res, w, outObs,outExp,IC,transform="none",tsparsemat,overdisp.est, cellsix_out, conf.level) {
    thetai <- res$Lambda_dense
    critval <- qnorm(1-(1-conf.level)/2)
    #print(critval)
    if(IC=="aic"){
        thetaa <- res$wLambda[res$selection.aic,][cellsix_out]
    } else {
        thetaa <- res$wLambda[res$selection.bic,][cellsix_out]
    }
    #print(str(thetaa))
    if(transform=="log"){
        #print("transform")
        if(!is.null(overdisp.est)){
            varthetai_1 <- as.matrix(as.vector(outYx)*sparsemat)
            varthetai_1[varthetai_1==0] <- Inf
            varthetai <- overdisp.est*(1/varthetai_1)
        } else {
            varthetai_1 <- as.matrix(as.vector(outYx)*sparsemat)
            varthetai_1[varthetai_1==0] <- Inf
            varthetai <- 1/varthetai_1
        }
        bb <- sapply(1:nrow(res$Lambda_dense), function(k) log(res$Lambda_dense[k,]) - log(round(res$wLambda[1,],3)))
        se_thetaa <- as.vector(sqrt(varthetai + bb^2)%*%matrix(w, ncol = 1))
        #
        #tmp[,j]=log(as.matrix(model_estimates)[,j])-log(as.matrix(model_estimates)%*%c(example_clusters$wgt2))
        #
        # withintheta <- sapply(1:length(cellsix_out), function(j) (log(thetai[,cellsix_out][,j]) - log(thetaa[j])))^2
        # wtnbtn <- sapply(1:length(cellsix_out), function(k) sqrt(varthetai[k] + withintheta[,k]))
        # varthetas_w <- matrix(w, nrow = 1)%*%wtnbtn
        # var_thetaa <- (as.vector(varthetas_w))^2
        LBa = exp(round(log(res$wLambda[1,],3)) -critval*(se_thetaa))
        UBa = exp(round(log(res$wLambda[1,],3)) +critval*(se_thetaa))
        # UBa = exp(as.vector(log(thetaa)) + critval*(se_thetaa))
        # LBa = exp(as.vector(log(thetaa)) - critval*(se_thetaa))
        
    } else {
        #browser()
        if(!is.null(overdisp.est)){
            varthetai_1 <- as.matrix(as.vector(outExp)*sparsemat)
            varthetai_1[varthetai_1==0] <- Inf
            varthetai <- overdisp.est*(thetai/t(varthetai_1))
        } else {
            varthetai <- (thetai/t(varthetai_1))#sapply(1:nrow(tsparsemat), function(k) thetai[k,]/outExp_out )
        }
        withintheta <- sapply(1:length(cellsix_out), function(j) (thetai[,cellsix_out][,j] - thetaa[j]))^2
        wtnbtn <- sapply(1:nrow(tsparsemat), function(k) sqrt(varthetai[k,cellsix_out] + withintheta[k,]))
        varthetas_w <- matrix(w, nrow = 1)%*%t(wtnbtn)
        var_thetaa <- (as.vector(varthetas_w))^2
        UBa = as.vector(thetaa) + critval*sqrt(var_thetaa)
        LBa = as.vector(thetaa) - critval*sqrt(var_thetaa)
    }
    return(list(buckland.LB = LBa,
                clusterstack= thetaa,
                buckland.UB = UBa))
}



#'@title ba2.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on Burnham and Anderson for specific cells (given by \code{cellsix_out}) 
#'@param thetaa Stacked relative risk estimate for the cluster(s).
#'@param res Resultant object from \code{detectclusters()}.
#'@param w Likelihood-weights for all potential clusters.
#'@param outObs Observed counts for each potential cluster.
#'@param outExp Expected counts for each potential cluster.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC).
#'@param transform If a transformation is used. Default is \code{"none"}. Other transformation currently available is \code{"log"}.
#'@param tsparsemat Transpose of large sparse matrix where columns are potential clusters and rows are space-time locations.
#'@param overdisp.est Estimate of overdispersion (or underdispersion).
#'@param cellsix_out Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate for each of the cells in \code{cellsix_out}
#'@export
#'
ba2.cells <- function(thetaa,res, w, outObs,outExp, IC, transform,tsparsemat, overdisp.est, cellsix_out, conf.level=0.95){
    thetai <- res$Lambda_dense
    critval <- qnorm(1-(1-conf.level)/2)
    if(IC=="aic"){
        thetaa <- res$wLambda[res$selection.aic,][cellsix_out]
    } else {
        thetaa <- res$wLambda[res$selection.bic,][cellsix_out]
    }
    if(transform=="log"){
        if(!is.null(overdisp.est)){
            varthetai_1 <- as.matrix(as.vector(outYx)*sparsemat)
            varthetai_1[varthetai_1==0] <- Inf
            varthetai <- overdisp.est*(1/varthetai_1)
        } else {
            varthetai_1 <- as.matrix(as.vector(outYx)*sparsemat)
            varthetai_1[varthetai_1==0] <- Inf
            varthetai <- (1/varthetai_1)
        }
        withintheta_A <- sapply(1:length(cellsix_out), function(j) (log(thetai[, cellsix_out][,j]) - log(thetaa[j])))^2
        wtnbtn_A <- sapply(1:nrow(tsparsemat), function(k) (varthetai[cellsix_out,k] + withintheta[k,]))
        #sqrt(varthetai[k] + withintheta[,k]))
            #sapply(1:length(cellsix_out), function(k) varthetai + withintheta[k,])
        var_thetaa_A <- as.vector(matrix(w, nrow=1)%*%wtnbtn)
        UBa_A = exp(as.vector(log(thetaa)) + critval*sqrt(var_thetaa))
        LBa_A = exp(as.vector(log(thetaa)) - critval*sqrt(var_thetaa))
    } else{
        if(!is.null(overdisp.est)){
            varthetai <- sapply(1:nrow(tsparsemat), function(k) overdisp.est*thetai[k,]/outExp_out)
        } else {
            varthetai <- sapply(1:nrow(tsparsemat), function(k) thetai[k,]/outExp_out)
        }
        withintheta <- sapply(1:length(cellsix_out), function(j) (thetai[,cellsix_out][,j] - thetaa[j]))^2
        wtnbtn <- sapply(1:nrow(tsparsemat), function(k) sqrt(varthetai[cellsix_out,k] + withintheta[k,]))
        var_thetaa <- as.vector(matrix(w, nrow=1)%*%t(wtnbtn))
        UBa = as.vector(thetaa) + critval*sqrt(var_thetaa)
        LBa = as.vector(thetaa) - critval*sqrt(var_thetaa)
    }
    return(list(ba2.LB = LBa,
                clusterstack= thetaa,
                ba2.UB = UBa))
}

#' @title mata_tailareazcore
#' @description Calculate MATA tail area weighted z-score based on Turek et al.
#' @param thetai Each individual cluster relative risk.
#' @param thetaa Stacked relative risk estimate.
#' @param sd.thetai Standard deviation of estimated variance for each individual cluster relative risk.
#' @param w Likelihood-weights for all Q<K potential clusters.
#' @param alpha \eqn{\alpha}, probability of rejecting the null hypothesis when it is true.
mata_tailareazscore <- function(thetai, thetaa, sd.thetai, w, alpha){
    thetai <- as.vector(thetai)
    zval <- (thetaa - thetai)/sd.thetai
    zpnorm <- pnorm(zval)
    w_zpnorm <- sum((w*zpnorm))-alpha
}


#'@title matabounds.cells
#'@description Wrapper function to calculate confidence bounds for stacked cluster relative risk estimates based on MATA (model-averaged tail area) intervals based on Turek et al. for specific cells (given by \code{cellsix_out}).
#'@param thetaa Stacked relative risk estimate for the cluster(s).
#'@param res Resultant object from \code{detectclusters()}.
#'@param w Likelihood-weights for all Q<K potential clusters.
#'@param outObs_out Observed counts for each potential cluster.
#'@param outExp_out Expected counts for each potential cluster.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC).
#'@param transform If a transformation is used. Default is \code{"none"}. Other transformation currently available are \code{"log"} and \code{"sqrt"}.
#'@param tsparsemat Transpose of large sparse matrix where columns are potential clusters and rows are space-time locations.
#'@param overdisp.est Estimate of overdispersion (or underdispersion).
#'@param cellsix_out Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate
#'@export
#'
matabounds.cells <- function(thetaa,res, w, outObs_out,outExp_out, IC, transform=c("none","log", "sqrt"),tsparsemat,overdisp.est, cellsix_out, conf.level=0.95) {
    if(is.null(transform)){
        none = matabounds_none.cells(thetaa, res, w,outObs_out, outExp_out, IC, tsparsemat, overdisp.est, cellsix_out, conf.level)
    }
    switch(transform,
           none = matabounds_none.cells(thetaa, res, w, outObs_out, outExp_out, IC, tsparsemat, overdisp.est, cellsix_out, conf.level),
           log = matabounds_log.cells(thetaa, res, w, outObs_out, outExp_out, IC, tsparsemat, overdisp.est, cellsix_out, conf.level),
           sqrt = matabounds_sqrt.cells(thetaa, res, w, outObs_out, outExp_out, IC, tsparsemat, overdisp.est, cellsix_out, conf.level))
}



#'@title matabounds_none.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on MATA (model-averaged tail area) intervals based on Turek et al. for specific cells (given by \code{cellsix_out}) with no transformation(\code{transform="none"})..
#'@param thetaa Stacked relative risk estimate for the cluster(s).
#'@param res Resultant object from \code{detectclusters()}.
#'@param w Likelihood-weights for all Q<K potential clusters.
#'@param outObs_out Observed counts for each potential cluster.
#'@param outExp_out Expected counts for each potential cluster.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC).
#'@param tsparsemat Transpose of large sparse matrix where columns are potential clusters and rows are space-time locations.
#'@param overdisp.est Estimate of overdispersion (or underdispersion).
#'@param cellsix_out Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate
#'@export
#'
matabounds_none.cells <- function(thetaa, res, w,outObs_out,outExp_out, IC, tsparsemat, overdisp.est,cellsix_out, conf.level=0.95) {
    thetai <- res$Lambda_dense
    alpha <- (1-conf.level)/2
    if(IC=="aic"){
        thetaa <- res$wLambda[res$selection.aic,][cellsix_out]
    } else {
        thetaa <- res$wLambda[res$selection.bic,][cellsix_out]
    }
    if(!is.null(overdisp.est)){
        varthetai <- sapply(1:nrow(tsparsemat), function(k) overdisp.est*(thetai[k,]/outExp_out))
    } else {
        varthetai <- sapply(1:nrow(tsparsemat), function(k) thetai[k,]/outExp_out)
    }
    mataLB <- sapply(1:length(cellsix_out),
                     function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                         thetai= thetai[,cellsix_out[k]],
                                         sd.thetai=sqrt(varthetai[cellsix_out[k],]),
                                         w=w, alpha=alpha, tol=1e-8)$root)
    
    mataUB <- sapply(1:length(cellsix_out),
                     function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                         thetai= thetai[,cellsix_out[k]],
                                         sd.thetai=sqrt(varthetai[cellsix_out[k],]),
                                         w=w, alpha=1-alpha, tol=1e-8)$root)
    
    return(list(mata.LB = mataLB,
                clusterstack= thetaa,
                mata.UB = mataUB))
}



#'@title matabounds_log.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on MATA (model-averaged tail area) intervals based on Turek et al. for specific cells (given by \code{cellsix_out}) with natural log transformation(\code{transform="log"}).
#'@param thetaa Stacked relative risk estimate for the cluster(s).
#'@param res Resultant object from \code{detectclusters()}.
#'@param w Likelihood-weights for all Q<K potential clusters.
#'@param outObs_out Observed counts for each potential cluster.
#'@param outExp_out Expected counts for each potential cluster.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC).
#'@param tsparsemat Transpose of large sparse matrix where columns are potential clusters and rows are space-time locations.
#'@param overdisp.est Estimate of overdispersion (or underdispersion).
#'@param cellsix_out Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate
#'@export
#'
matabounds_log.cells<- function(thetaa, res, w,outObs_out,outExp_out, IC, tsparsemat, overdisp.est,cellsix_out, conf.level=0.95) {
    thetai <- res$Lambda_dense
    alpha <- (1-conf.level)/2
    if(IC=="aic"){
        thetaa <- res$wLambda[res$selection.aic,][cellsix_out]
    } else {
        thetaa <- res$wLambda[res$selection.bic,][cellsix_out]
    }
    if(!is.null(overdisp.est)){
        #logTvarthetai <- sapply(1:nrow(tsparsemat), function(k) overdisp.est*(1/outObs_out))
        varthetai_1 <- as.matrix(as.vector(outYx)*sparsemat)
        varthetai_1[varthetai_1==0] <- Inf
        logTvarthetai  <- overdisp.est*(1/varthetai_1)
    } else{
        #logTvarthetai <- sapply(1:nrow(tsparsemat), function(k) 1/(thetai[k,]*outObs_out))
        varthetai_1 <- as.matrix(as.vector(outYx)*sparsemat)
        varthetai_1[varthetai_1==0] <- Inf
        logTvarthetai  <- (1/varthetai_1)
    }
    
    
    mataLB <- sapply(1:length(cellsix_out),
                     function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                         thetai= log(thetai[,cellsix_out[k]]),
                                         sd.thetai=sqrt(logTvarthetai[cellsix_out[k],]),
                                         w=w, alpha=alpha, tol=1e-8)$root)
    mataUB <- sapply(1:length(cellsix_out),
                     function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                         thetai= log(thetai[,cellsix_out[k]]),
                                         sd.thetai=sqrt(logTvarthetai[cellsix_out[k],]),
                                         w=w, alpha=1-alpha, tol=1e-8)$root)
    return(list(matalog.LB = exp(mataLB),
                clusterstack= thetaa,
                matalog.UB = exp(mataUB)))
}






#' @title selectuniqRR
#' @description Helper function for \code{calcbounds()}. Identifies unique relative risks across potential clusters that are not equal to 1.
#' @param uniqRRs Matrix of unique relative risks for each potential cluster and background.
selectuniqRR <- function(uniqRRs){
    clusterRR_i <- rep(NA, nrow(uniqRRs))
    flag1 <- sapply(1:nrow(uniqRRs), function(k) length(unique(uniqRRs[k,])))
    non1s <- uniqRRs[which(flag1==2),]
    clusterRR_i[which(flag1==2)] <- sapply(1:nrow(non1s), function(k) non1s[k,which(non1s[k,]!=1)])
    clusterRR_i[which(flag1==1)] <- 1    
    return(unlist(clusterRR_i))
    
}


#'@title calcbounds
#'@description Calculates lower and upper confidence bounds.
#'@param id_ic The number of clusters identified either by BIC (QBIC) or AIC (QAIC).
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC). This should match \code{id_ic} argument.
#'@param res Resultant object from \code{detectclusters()}.
#'@param byloc Boolean. Use \code{TRUE} when stacking by location. Use \code{FALSE} when stacking by potential cluster.
#'@param Ex Expected counts (unstandardized) for each cell.
#'@param Yx Observed counts for each cell.
#'@param cellsix_out Indices of the cells to calculate bounds for. 
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@details The confidence bounds methods returned are:
#'\itemize{
#'\item \code{buck}: Stacked confidence bounds based on Buckland et al.
#'\item \code{buckTlog}: Stacked confidence bounds based on Buckland et al.; natural log-transformed.
#'\item \code{ba}: Model-average weighted confidence bounds based on Burnham and Anderson
#'\item \code{baTlog}: Model-average weighted confidence bounds based on Burnham and Anderson.
#'\item \code{mata}: Model-average tail area (MATA) intervals based on Turek et al.
#'\item \code{mataTlog}: Model-average tail area (MATA) intervals based on Turek et al.; natural log-transformed.
#'
#'}
calcbounds <- function(id_ic, IC, res, byloc, Ex, Yx, cellsix=NULL, sparsemat, conf.level=0.95){
    IC <- tolower(IC)
    if(!is.null(cellsix) & is.null(sparsemat)){
        stop("You must provide the sparsemat when calculating rates for each cell.")
    }
    if(is.null(cellsix)){
        stop("For cell-wise estimates, you must provide the index of the cells of interest")
    } 
    thetai_uniq <- sapply(1:nrow(res$Lambda_dense), function(k) unique(res$Lambda_dense[k,]))
    thetai_uniqi <- as.matrix(do.call(rbind, thetai_uniq), ncol=2)
    thetai <- selectuniqRR(thetai_uniqi)
    if(IC=="aic"){
        out <- vector(mode = "list", length = res$selection.aic)
        w <- matrix(res$wtMAT[,1:res$selection.aic], ncol=res$selection.aic)
        for(i in 1:res$selection.aic){
            thetaa <- sum(thetai*w[,i])
            out[[i]] <- calcbounds.cells(id_ic, IC, res, byloc, Ex, Yx,as.vector(w[,i]), thetaa,thetai, sparsemat, cellsix, conf.level)
            
        }
    } else {
        out <- vector(mode = "list", length = res$selection.bic)
        w <- matrix(res$wtMAT[,1:res$selection.bic], ncol=res$selection.bic)
        for(i in 1:res$selection.bic){
            thetaa <- sum(thetai*w[,i])
            out[[i]] <-calcbounds.cells(id_ic, IC, res, byloc, Ex, Yx,as.vector(w[,i]), thetaa,thetai, sparsemat, cellsix, conf.level)
            
        }
        
    }
    return(out)
}





#'@title calcbounds.cells
#'@description Helper function for \code{calcbounds()}.
#'@param id_ic The number of clusters identified either by BIC (QBIC) or AIC (QAIC).
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC). This should match \code{id_ic} argument.
#'@param res Resultant object from \code{detectclusters()}.
#'@param byloc Boolean. Use \code{TRUE} when stacking by location. Use \code{FALSE} when stacking by potential cluster.
#'@param Ex Expected counts (unstandardized) for each cell.
#'@param Yx Observed counts for each cell.
#'@param w Likelihood-weights for all potential clusters.
#'@param thetaa Stacked relative risk estimate for the cluster(s)
#'@param thetai Individual relative risk estimates for each potential cluster.
#'@param sparsemat Large sparse matrix where rows are potential clusters and columns are space-time locations.
#'@param cellsix Indices of the cells to calculate bounds for. 
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return Returns large list of confidence bounds and stacked estimates in addition to timings for each of the confidence bounds methods.    
calcbounds.cells <- function(id_ic, IC, res, byloc=FALSE, Ex, Yx,w, thetaa,thetai, sparsemat, cellsix, conf.level=0.95){
    browser()
    if(is.null(conf.level)){conf.level <- 0.95}
    if(id_ic==0){
        # emptynonstack<- vector(mode = "list", length = 2)
        # emptynonstack[[1]] <- 0
        # emptynonstack[[2]] <- matrix(rep(0, length(cellsix)*3), ncol=3)
        # outnonstack<- outnonstackTlog <- outnonstack_asymp <- outnonstack_asympTlog <- emptynonstack
        # names(outnonstack) <- c("nonstack.theta.time", "nonstack.theta")
        # names(outnonstackTlog) <- c("nonstack.theta.time", "nonstack.theta")
        # names(outnonstack_asymp) <- c("nonstack_asymp.theta.time", "nonstack_asymp.theta")
        # names(outnonstack_asympTlog) <- c("nonstack_asymp.theta.time", "nonstack_asymp.theta")
        return(list(
            # outnonstack= outnonstack,
            # outnonstackTlog = outnonstackTlog,
            # outnonstack_asymp = outnonstack_asymp,
            # outnonstack_asympTlog = outnonstack_asympTlog,
            outbuck.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            outbuckTlog.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            outba2.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            outba2Tlog.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            outmata.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            outmataTlog.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            
            # outnonstack.time = 000,
            # outnonstackTlog.time = 000,
            # outnonstack_asymp.time = 000,
            # outnonstack_asympTlog.time = 000,
            outbuck.theta.time = 000,
            outbuckTlog.theta.time = 000,
            outba2.theta.time = 000,
            outba2Tlog.theta.time = 000,
            outmata.theta.time = 000,
            outmataTlog.theta.time = 000
        ))
        
        
    } else {
        cellrisk_wt_out <- rep(NA, length(cellsix))
        cellsix_out <- cellsix
        if(byloc==TRUE){
            #print(paste0("byloc", byloc))
            outExp <- t(sparsemat)%*%Ex
            outYx <- t(sparsemat)%*%Yx
            #outExp_out <- Ex[cellsix]
            #outYx_out <- Yx[cellsix]
            for(i in 1:length(cellsix)){
                cellsixvec <- rep(0,dim(res$wLambda)[2])
                cellsixvec[cellsix[i]] <-1
                overlapid <- matrix(cellsixvec, nrow=1)%*%t(res$Lambda_dense)
                cellrisk_wt <- overlapid%*%w 
                cellrisk_wt_out[[i]] <- cellrisk_wt
            }
        }
        else if (byloc==FALSE){
            #print(paste0("byloc: ", byloc))
            outExp <- t(sparsemat)%*%Ex #66870 x 1
            outYx <- t(sparsemat)%*%Yx #66870 x 1
            #outExp_out <- rep(outExp@x[res$maxid[id_ic]], times=length(cellsix))
            #outObs_out <- rep(outYx@x[res$maxid[id_ic]], times=length(cellsix))
            cellrisk_wt_out <- NULL
        }    
        
     
        
        outbuck.theta.time <- system.time(outbuck.theta <- bucklandbounds.cells(thetaa,
                                                                                res,
                                                                                w=w,
                                                                                #Ex,
                                                                                outObs,
                                                                                outExp,
                                                                                IC=IC,
                                                                                transform="none",
                                                                                tsparsemat=t(sparsemat),
                                                                                overdisp.est, cellsix_out,
                                                                                conf.level))
        outbuckTlog.theta.time <- system.time(outbuckTlog.theta <- bucklandbounds.cells(thetaa,
                                                                                        res,
                                                                                        w=w,
                                                                                        #Ex,
                                                                                        outObs_out,
                                                                                        outExp_out,
                                                                                        IC=IC,
                                                                                        transform="log",
                                                                                        tsparsemat=t(sparsemat),
                                                                                        overdisp.est, cellsix_out,
                                                                                        conf.level))
        message("Buckland bounds finished")
        outba2.theta.time <- system.time(outba2.theta <- ba2.cells(thetaa,
                                                                      res,
                                                                      w=w,
                                                                      #Ex,
                                                                   outObs_out,
                                                                   outExp_out,
                                                                      IC,
                                                                      transform="none",
                                                                      tsparsemat=t(sparsemat),
                                                                      overdisp.est,
                                                                      cellsix_out,
                                                                      conf.level))
        outba2Tlog.theta.time <- system.time(outba2Tlog.theta <- ba2.cells(thetaa,
                                                                              res,
                                                                              w=w,
                                                                              #Ex,
                                                                           outObs_out,
                                                                           outExp_out,
                                                                              IC,
                                                                              transform="log",
                                                                              tsparsemat=t(sparsemat),
                                                                              overdisp.est,
                                                                              cellsix_out,
                                                                              conf.level))
        
        message("Burnham & Anderson bounds finished")
        outmata.theta.time <- system.time(outmata.theta <- matabounds.cells(thetaa,
                                                                            res,
                                                                            w=w,
                                                                            #Ex,
                                                                            outObs_out,
                                                                            outExp_out,
                                                                            IC,
                                                                            transform="none",
                                                                            tsparsemat=t(sparsemat),
                                                                            overdisp.est,
                                                                            cellsix_out,
                                                                            conf.level))
        outmataTlog.theta.time <- system.time(outmataTlog.theta <- matabounds.cells(thetaa,
                                                                                    res,
                                                                                    w=w,
                                                                                    #Ex,
                                                                                    outObs_out,
                                                                                    outExp_out,
                                                                                    IC,
                                                                                    transform="log",
                                                                                    tsparsemat=t(sparsemat),
                                                                                    overdisp.est,
                                                                                    cellsix_out,
                                                                                    conf.level))
        message("MATA bounds finished")
        
        
        return(list(
            # outnonstack= outnonstack,
            # outnonstackTlog = outnonstackTlog,
            # outnonstack_asymp = outnonstack_asymp,
            # outnonstack_asympTlog = outnonstack_asympTlog,
            outbuck.theta = outbuck.theta,
            outbuckTlog.theta = outbuckTlog.theta,
            outba2.theta = outba2.theta,
            outba2Tlog.theta = outba2Tlog.theta,
            outmata.theta = outmata.theta,
            outmataTlog.theta = outmataTlog.theta,
            
            outnonstack.time = outnonstack.time[[3]],
            outnonstackTlog.time = outnonstackTlog.time[[3]],
            outnonstack_asymp.time = outnonstack_asymp.time[[3]],
            outnonstack_asympTlog.time = outnonstack_asympTlog.time[[3]],
            outbuck.theta.time = outbuck.theta.time[[3]],
            outbuckTlog.theta.time = outbuckTlog.theta.time[[3]],
            outba2.theta.time = outba2.theta.time[[3]],
            outba2Tlog.theta.time = outba2Tlog.theta.time[[3]],
            outmata.theta.time = outmata.theta.time[[3]],
            outmataTlog.theta.time = outmataTlog.theta.time[[3]]
        ))
    }
    
}

