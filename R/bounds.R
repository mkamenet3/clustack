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
calcbounds <- function(id_ic, IC, res, byloc, Ex, Yx, cellsix=NULL, sparsemat, conf.level=0.95, overdisp.est=NULL){
    #browser()
    if(is.null(overdisp.est)){overdisp.est<-NULL}
    IC <- tolower(IC)
    if(!is.null(cellsix) & is.null(sparsemat)){
        stop("You must provide the sparsemat when calculating rates for each cell.")
    }
    if(is.null(cellsix)){
        stop("For cell-wise estimates, you must provide the index of the cells of interest")
    } 
    thetai_uniq <- sapply(1:nrow(res$Lambda_dense), function(k) unique(res$Lambda_dense[k,]))
    thetai_uniqi <- as.matrix(do.call(rbind, thetai_uniq), ncol=2)
    thetai <- selectuniq(thetai_uniqi,1)
    if(IC=="aic"){
        select.aic <- max(1, res$selection.aic)
        out <- vector(mode = "list", length = select.aic)
        if(id_ic==0){
            out <- poisLik(Ex, Yx, t(sparsemat))
            Lik <- out$Lik
            w <- matrix(likweights(Lik),ncol=select.aic)
        } else{
            w <- matrix(res$wtMAT[,1:select.aic], ncol=select.aic)  
        }
        for(i in 1:select.aic){
            thetaa <- sum(thetai*w[,i])
            out[[i]] <- calcbounds.cells(id_ic, IC, res, byloc, Ex, Yx,as.vector(w[,i]), thetaa,thetai, sparsemat, cellsix, conf.level,overdisp.est)
            
        }
    } else {
        select.bic <- max(1, res$selection.bic)
        out <- vector(mode = "list", length = select.bic)
        if(id_ic==0){
            out <- poisLik(Ex, Yx, t(sparsemat))
            Lik <- out$Lik
            w <- matrix(likweights(Lik), ncol=select.bic)
        }else{
            w <- matrix(res$wtMAT[,1:select.bic], ncol=select.bic) 
        }
        for(i in 1:select.bic){
            thetaa <- sum(thetai*w[,i])
            out[[i]] <-calcbounds.cells(id_ic, IC, res, byloc, Ex, Yx,as.vector(w[,i]), thetaa,thetai, sparsemat, cellsix, conf.level,overdisp.est)
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
calcbounds.cells <- function(id_ic, IC, res, byloc=FALSE, Ex, Yx,w, thetaa,thetai, sparsemat, cellsix, conf.level=0.95, overdisp.est){
    #browser()
    if(is.null(conf.level)){conf.level <- 0.95}
    critval <- qnorm(1-(1-conf.level)/2)
    alpha2 <-(1-conf.level)/2
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
        
        ###########################
        #Calculate key matrices
        ###########################
        thetai <- res$Lambda_dense
        thetai_uniqvals <- apply(thetai,1, unique)
        #browser()
        thetai_uniqmat <- data.matrix(as.data.frame(thetai_uniqvals))
        thetai_uniq <- selectuniq(t(thetai_uniqmat), 1) #66870
        message("Extracted unique thetai's")
        
        thetaa <- res$wLambda[res$selection.bic,][cellsix]
        Yx_pcloc <- t(sparsemat)%*%Yx 
        
        if(!is.null(overdisp.est)){
            se_thetai <- sqrt(overdisp.est*(1/Yx_pcloc))
        }else{
            se_thetai <- sqrt(1/Yx_pcloc)
        }
        #replace Inf with 0
        ixinf <- which(is.infinite(se_thetai))
        se_thetai[ixinf] <-0
        
        
        withinvar <- sapply(1:length(cellsix), function(i) (log(thetai_uniq) - log(thetaa[i])))
        ses <- sapply(1:length(cellsix), function(i) (se_thetai)^2 + (withinvar[,i])^2)
        var_est<- sapply(1:length(cellsix), function(i) sparsemat %*% Diagonal(length(as.vector(ses[[i]])),
                                                                               as.vector(ses[[i]])))
        
        
        message("Finished calculating SEs")
        
        #browser()
        ###########################
        #Calculate Bounds by Method
        ###########################
        #thetaa, ses, w, sparsemat,cellsix, critval
        outbuckTlog.theta.time <- system.time(outbuckTlog.theta <- bucklandbounds.cells(thetaa, 
                                                                                        var_est, 
                                                                                        w, 
                                                                                        sparsemat,
                                                                                        cellsix,
                                                                                        critval))
        
        message("Buckland bounds finished")
        outba2Tlog.theta.time <- system.time(outba2Tlog.theta <- ba2.cells(thetaa, 
                                                                           var_est, 
                                                                           w, 
                                                                           sparsemat,
                                                                           cellsix, 
                                                                           critval))
        
        message("Burnham & Anderson bounds finished")
        outmataTlog.theta.time <- system.time(outmataTlog.theta <- matabounds_log.cells(var_est, thetai, thetaa, w,
                                                                                        cellsix,alpha2))
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
bucklandbounds.cells <- function(thetaa, var_est, w, sparsemat,cellsix, critval){
    var_estw <- sapply(1:length(cellsix), function(i) sqrt(var_est[[i]])%*%w)
    LBa = sapply(1:length(cellsix), function(i) exp(log(thetaa[i])-critval*(var_estw[[i]][cellsix[i]])))
    UBa = sapply(1:length(cellsix), function(i) exp(log(thetaa[i])+critval*(var_estw[[i]][cellsix[i]])))
    
    return(list(buckland.LB = LBa,
                clusterstack= thetaa,
                buckland.UB = UBa))
}

#'@title ba2.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on Burnham and Anderson for specific cells (given by \code{cellsix}) 
ba2.cells <- function(thetaa, var_est, w, sparsemat,cellsix, critval){
    var_estw <- sapply(1:length(cellsix), function(i) var_est[[i]]%*%w)
    LBa = sapply(1:length(cellsix), function(i) exp(log(thetaa[i])-critval*(sqrt(var_estw[[i]][cellsix[i]]))))
    UBa = sapply(1:length(cellsix), function(i) exp(log(thetaa[i])+critval*(sqrt(var_estw[[i]][cellsix[i]]))))
    
    return(list(ba2.LB = LBa,
                clusterstack= thetaa,
                ba2.UB = UBa))
}

#########MATA Bounds
mata_ZLB <- function(theta, sesi,thetaii,w){
    resLB <-sum(w*(1-pnorm((log(thetaii)-theta)/(sesi+(log(thetaii)-theta==0)))))
    
}

mata_ZUB <- function(theta,sesi, thetaii,w){
    1-mata_ZLB(theta,sesi, thetaii,w)
}
f_LB<-function(theta,sesi,thetaii,w, alpha2){
    mata_ZLB(theta, sesi,thetaii,w)-alpha2  
} 
f_UB<-function(theta,sesi,thetaii,w,alpha2){
    mata_ZUB(theta,sesi, thetaii,w)-alpha2  
} 

matabounds_log.cells <- function(var_est,thetai, thetaa,w, cellsix,alpha2){
    #browser()
    se_est <- lapply(1:length(cellsix), function(i) sqrt(var_est[[i]]))
    LBa <- sapply(1:length(cellsix), function(i) exp(uniroot(f_LB,c(-3,3),
                                                             sesi=se_est[[i]][cellsix[i],],
                                                             thetaii=thetai[,cellsix[i]],
                                                             w=w,
                                                             alpha2=alpha2)$root))
    UBa <- sapply(1:length(cellsix), function(i) exp(uniroot(f_UB,c(-3,3),
                                                             sesi=se_est[[i]][cellsix[i],],
                                                             thetaii=thetai[,cellsix[i]],
                                                             w=w,
                                                             alpha2=alpha2)$root))
    return(list(mata.LB = LBa,
            clusterstack= thetaa,
                mata.UB = UBa))

}
 
 
