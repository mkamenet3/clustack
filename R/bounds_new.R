library(matrixcalc)
####################################################
selectuniq <- function(uniqRRs, background){
    clusterRR_i <- rep(NA, nrow(uniqRRs))
    flag1 <- sapply(1:nrow(uniqRRs), function(k) length(unique(uniqRRs[k,])))
    non1s <- uniqRRs[which(flag1==2),]
    clusterRR_i[which(flag1==2)] <- sapply(1:nrow(non1s), function(k) non1s[k,which(non1s[k,]!= background)])
    clusterRR_i[which(flag1==1)] <-  background    
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
#'@param cellsix Indices of the cells to calculate bounds for. 
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@details The confidence bounds methods returned are:
#'\itemize{
#'\item \code{buckTlog}: Stacked confidence bounds based on Buckland et al.; natural log-transformed.
#'\item \code{baTlog}: Model-average weighted confidence bounds based on Burnham and Anderson.
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
    thetai <- selectuniq(thetai_uniqi)
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
   # browser()
    if(is.null(conf.level)){conf.level <- 0.95}
    critval <- qnorm(1-(1-conf.level)/2)
    if(id_ic==0){
        return(list(
            outbuckTlog.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            outba2Tlog.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            outmataTlog.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),

            outbuckTlog.theta.time = 000,
            outba2Tlog.theta.time = 000,
            outmataTlog.theta.time = 000
        ))
    } else {
        cellrisk_wt_out <- rep(NA, length(cellsix))
        if(byloc==TRUE){
            outExp <- t(sparsemat)%*%Ex
            outYx <- t(sparsemat)%*%Yx
            for(i in 1:length(cellsix)){
                cellsixvec <- rep(0,dim(res$wLambda)[2])
                cellsixvec[cellsix[i]] <-1
                overlapid <- matrix(cellsixvec, nrow=1)%*%t(res$Lambda_dense)
                cellrisk_wt <- overlapid%*%w 
                cellrisk_wt_out[[i]] <- cellrisk_wt
            }
        }
        else if (byloc==FALSE){
            cellrisk_wt_out <- NULL
        }    
        ###########################
        #Calculate key matrices
        ###########################
        thetai <- res$Lambda_dense
        thetai_uniqvals <- apply(thetai,1, function(x) unique)
        thetai_uniqmat <- data.matrix(as.data.frame(thetai_uniqvals))
        thetai_uniq <- selectuniqRR(t(thetai_uniqmat), 1) #66870
        
        
        thetaa <- res$wLambda[res$selection.bic,][cellsix]
        Yx_pc <- t(sparsemat)%*%Yx 
        Yx_pcloc <- as.matrix(sweep(sparsemat, MARGIN=2, Yx_pc@x, `*`)) 
        Yx_pcloc[Yx_pcloc==0] <- Inf
        if(!is.null(overdisp.est)){
            var_thetai <- overdisp.est*(1/Yx_pcloc)
        }else{
            var_thetai <- 1/Yx_pcloc    
        }
        
        var_thetai_uniqvals <- apply(var_thetai, MARGIN=2, unique)
        var_thetai_uniqmat <- data.matrix(as.data.frame(moo2))
        var_thetai_uniq <- selectuniq(t(var_thetai_uniqmat),0)
        
        
        withinvar <- sapply(1:length(cellsix), function(i) (log(thetai_uniq) - log(thetaa[i]))^2)
        
        
        
        ###########################
        #Calculate Bounds by Method
        ###########################
        outbuckTlog.theta.time <- system.time(outbuckTlog.theta <- bucklandbounds.cells(thetaa, 
                                                                                        var_thetai_uniq, 
                                                                                        withinvar, 
                                                                                        w, 
                                                                                        sparsemat,
                                                                                        cellsix,
                                                                                        critval))
 
        message("Buckland bounds finished")
        outba2Tlog.theta.time <- system.time(outba2Tlog.theta <- ba2.cells(thetaa, 
                                                                           var_thetai_uniq, 
                                                                           withinvar, 
                                                                           w, 
                                                                           sparsemat,
                                                                           cellsix, 
                                                                           critval))
        
        message("Burnham & Anderson bounds finished")
        outmataTlog.theta.time <- system.time(outmataTlog.theta <- matabounds_log.cells(log(thetaa),
                                                                                        log(thetai),
                                                                                        var_thetai_uniq, 
                                                                                        w, 
                                                                                        sparsemat,
                                                                                        cellsix, 
                                                                                        conf.level))
        message("MATA bounds finished")
        
        ###########################
        #Return Results
        ###########################
        return(list(
            outbuckTlog.theta = outbuckTlog.theta,
            outba2Tlog.theta = outba2Tlog.theta,
            outmataTlog.theta = outmataTlog.theta,
            
            outbuckTlog.theta.time = outbuckTlog.theta.time[[3]],
            outba2Tlog.theta.time = outba2Tlog.theta.time[[3]],
            outmataTlog.theta.time = outmataTlog.theta.time[[3]]
        ))
    }
    
}



#'@title bucklandbounds.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on Buckland et al. for specific cells (given by \code{cellsix}) 
#'@param thetaa Stacked relative risk estimates for each cell.
#'@param var_thetai_uniq Unique values of variance for each unique individual cell relative risk estimate.
#'@param withinvar Within estimate variablility (between each unique thetai and the stacked estimate for the cell)
#'@param w Likelihood-weights for potential clusters.
#'@param sparsemat Large sparse matrix where columns are potential clusters and rows are space-time locations.
#'@param cellsix Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate for each of the cells in \code{cellsix}
#'@export
bucklandbounds.cells <- function(thetaa, var_thetai_uniq, withinvar, w, sparsemat,cellsix, critval){
    combo <- sapply(1:length(cellsix), function(i) var_thetai_uniq+withinvar[,i])
    var_est<- sapply(1:length(cellsix), function(i) t(t(sparsemat)*combo[,i]))
    var_estw <- sapply(1:length(cellsix), function(i) sqrt(var_est[[i]])%*%w)
    LBa = sapply(1:length(cellsix), function(i) exp(log(thetaa[i])-critval*(var_estw[[i]][cellsix[i]])))
    UBa = sapply(1:length(cellsix), function(i) exp(log(thetaa[i])+critval*(var_estw[[i]][cellsix[i]])))
    
    return(list(buckland.LB = LBa,
                clusterstack= thetaa,
                buckland.UB = UBa))
}

#'@title ba2.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on Burnham and Anderson for specific cells (given by \code{cellsix}) 
ba2.cells <- function(thetaa, var_thetai_uniq, withinvar, w, sparsemat,cellsix, critval){
    combo <- sapply(1:length(cellsix), function(i) var_thetai_uniq+withinvar[,i])
    var_est<- sapply(1:length(cellsix), function(i) t(t(sparsemat)*combo[,i]))
    var_estw <- sapply(1:length(cellsix), function(i) var_est[[i]]%*%w)
    LBa = sapply(1:length(cellsix), function(i) exp(log(thetaa[i])-critval*(sqrt(var_estw[[i]][cellsix[i]]))))
    UBa = sapply(1:length(cellsix), function(i) exp(log(thetaa[i])+critval*(sqrt(var_estw[[i]][cellsix[i]]))))
    
    return(list(ba2.LB = LBa,
                clusterstack= thetaa,
                ba2.UB = UBa))
}

###########################

#'@title matabounds_log.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on MATA (model-averaged tail area) intervals based on Turek et al. for specific cells (given by \code{cellsix}) with natural log transformation(\code{transform="log"}).
#'@param thetaa Stacked relative risk estimate for the cluster(s).
#'@param res Resultant object from \code{detectclusters()}.
#'@param w Likelihood-weights for all Q<K potential clusters.
#'@param outObs_out Observed counts for each potential cluster.
#'@param outExp_out Expected counts for each potential cluster.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC).
#'@param tsparsemat Transpose of large sparse matrix where columns are potential clusters and rows are space-time locations.
#'@param overdisp.est Estimate of overdispersion (or underdispersion).
#'@param cellsix Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate
#'@export
#'
matabounds_log.cells<- function(thetaa, thetai,var_thetai_uniq, w, sparsemat,cellsix, conf.level) {
    alpha <- (1-conf.level)/2
    mataLB <- sapply(1:length(cellsix),
                     function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                         thetai= thetai[,cellsix[k]],
                                         sd.thetai=sqrt(var_thetai_uniq),
                                         w=w, alpha=alpha, tol=1e-8)$root)
    mataUB <- sapply(1:length(cellsix),
                     function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                         thetai= thetai[,cellsix[k]],
                                         sd.thetai=sqrt(var_thetai_uniq),
                                         w=w, alpha=1-alpha, tol=1e-8)$root)
    return(list(matalog.LB = exp(mataLB),
                clusterstack= exp(thetaa),
                matalog.UB = exp(mataUB)))
}

#' @title mata_tailareazcore
#' @description Calculate MATA tail area weighted z-score based on Turek et al.
#' @param thetaa Stacked relative risk estimate.
#' @param thetai Each individual cluster relative risk.
#' @param sd.thetai Standard deviation of estimated variance for each individual cluster relative risk.
#' @param w Likelihood-weights for all Q<K potential clusters.
#' @param alpha \eqn{\alpha}, probability of rejecting the null hypothesis when it is true.
mata_tailareazscore <- function(thetaa,thetai, sd.thetai, w, alpha){
    zval <- (thetaa - thetai)/sd.thetai
    zpnorm <- pnorm(zval)
    w_zpnorm <- sum((w*zpnorm))-alpha
}

