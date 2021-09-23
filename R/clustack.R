# clustack functions
#Load Libraries
library(clusso)
library(testthat)
library(Matrix)


############################################################


#' @title likeweights
#' @description Generate likelihood-based initial weights.
#' @param liki Likelihood for each potential cluster.
#' @return A vector of weights.
likweights <- function(liki){
    wi <- liki/sum(liki)
    if(any(is.na(liki))){
        message("There are NA likelihoods")
    }
    return(wi)
}


#' @title overdisp
#' @description Calculate the overdispersion parameter (\eqn{\phi}),
#' @param offset_reg Object of class glm or lm.
#' @param overdispfloor Boolean. When \code{TRUE}, overdispersion parameter (\eqn{\phi}) is limited to not be less than 1. If \code{FALSE}, underdispersion can be estimated. 
#' @return An estimate of \eqn{\phi} (overdispersion or underdispersion)
overdisp <- function(offset_reg, overdispfloor = TRUE) {
    stopifnot(inherits(offset_reg, c("glm", "lm")))
    phi <- max(unlist(deviance(offset_reg)/df.residual(offset_reg)))
    if(overdispfloor == TRUE & phi < 1){
        message(paste("Underdispersion detected (", phi,"). Setting phi to 1"))
        phi <- 1
    }
    if(overdispfloor == FALSE & phi < 1){
        message(paste("Underdispersion detect (", phi,").\n 'Floor' argument was set to FALSE to underdispersed model will be run"))
    }
    return(phi)
}

#'@title poisLik
#'@description Poisson-based LRT comparing cluster model to null model.
#'@param Ex Vector of expected counts.
#'@param Yx Vector of observed counts.
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@export
#'@return Likelihood for each potential cluster.
poisLik <-function(Ex, Yx, sparsemat){
    outExp <- sparsemat%*%Ex
    outObs <- sparsemat%*%Yx
    #calc Lambda
    lambdahat <- outObs/outExp
    Lambda <- as.vector(lambdahat)*sparsemat #big Lambda matrix
    Lambda_dense <- as.matrix(Lambda)
    Lambda_dense[Lambda_dense == 0] <- 1
    #Get scaled likelihood
    outlogLik <- outObs*log(outObs/outExp) + (outExp-outObs)
    outlogLik <- ifelse(is.na(outlogLik@x),0, outlogLik@x) #likelihood should be zero here
    outlogLik_scaled <- outlogLik-max(outlogLik)
    Lik <- exp(outlogLik_scaled)
    return(list(Lik=Lik,
                Lambda_dense=Lambda_dense))
}

#' #'@title binomLik
#' #'@description Binomial-based likelihood LRT comparing cluster model to null model.
#' #'@param Ex Vector of expected counts.
#' #'@param Yx Vector of observed counts.
#' #'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#' #'@export
#' #'@return Likelihood for each potential cluster.
#' binomLik <-function(nx, Yx, sparsemat){
#'     #outExp <- sparsemat%*%Ex
#'     outnx <- sparsemat%*%nx
#'     outObs <- sparsemat%*%Yx
#'     #calc Lambda
#'     lambdahat <- outObs/outExp
#'     Lambda <- as.vector(lambdahat)*sparsemat #big Lambda matrix
#'     Lambda_dense <- as.matrix(Lambda)
#'     Lambda_dense[Lambda_dense == 0] <- 1
#'     #Get scaled likelihood
#'     
#'     #Lik <- (((outObs/nx)/(sum(outObs)/sum(nx)))^(outobs))*(((sum(outObs)-outObs)/(sum(nx)-nx))/(sum(outObs)/sum(nx)))^(sum(outObs)-outObs)
#'     #((outObs/outExp)/(sum(outObs)/sum(outExp)))^outObs #TODO CHECK THIS
#'     outnt <- sum(outnx)
#'     outObst <- sum(outObs)
#'     Lik <- (((outObs/outnx)^outObs)*(1-(outObs/outnx))^(outnx-outObs))*(((outObst - outObs)/(outnt-outnx))^(outObst-outObs))*(1-((outObst-outObs)/(outnt - outnx)))^(outnt - outnx - outObst+outObs)
#'     
#'     outlogLik <- log(Lik)
#'     outlogLik_scaled <- outlogLik-max(outlogLik)
#'     Lik <- exp(outlogLik_scaled)
#'     return(list(Lik=Lik,
#'                 Lambda_dense=Lambda_dense))
#' }

#'@title bylocation
#'@description Cluster detection based on maximum location.
#'@param Lik Likelihood for each potential cluster.
#'#'@param Lambda_dense Initial estimates for relative risk for each potential cluster;\eqn{\Lambda}.
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@param maxclust Maximum number of clusters to be detected.
#'@return List. First element is a Weighted relative risks for identified locations. Second element is a large matrix of weights for each iteration.
#'
bylocation <- function(Lik, Lambda_dense,sparsemat, maxclust){
    wtMAT <- matrix(rep(NA, maxclust*nrow(sparsemat)), ncol=maxclust)
    wt0 <- likweights(Lik)
    ixall <- NULL
    ix <- NULL
    maxlocs <- rep(NA, maxclust)
    stop <- FALSE
    for(i in 1:maxclust){
        wtmp <-wt0
        wtmp[ixall] <-0
        #Find maxloc
        wi_loc <- t(wtmp)%*%sparsemat
        maxloc <- which.max(as.vector(wi_loc))
        message(paste0("Location identified: ",(maxloc)))
        maxlocs[i] <- maxloc 
        #find all potential clusters that overlap that location
        locmax <- rep(0,numCenters*Time); 
        locmax[maxloc] <-1;
        locmax <- matrix(locmax,ncol=1)
        pclocmax <- as.vector(t(locmax)%*%t(sparsemat))
        ix <- which(pclocmax!=0) #indices of all PCs that overlap max location
        #upweight Lik* to 1 in all Pcs that overlap max cluster
        #reweight everything inside cluster to sum to 1
        wtmp[ix] <- likweights(Lik[ix])
        #save wtmp vector into master weight matrix
        wtMAT[,i] <- wtmp
        if(all(ix %in% ixall)){
            stop <- TRUE
            maxi = i
            break
        }
        if(stop){break}
        #save everything that's already been considered
        ixall <- c(ixall, ix)
        maxi = i 
        #Lik[ixall] <-0
    }
    wLambda <- crossprod(wtMAT[,1:(maxi)], Lambda_dense)
    return(list(wLambda = wLambda,
                wtMAT = wtMAT,
                #wtMAT0 = wtMAT0,
                maxlocs = maxlocs))
}




#'@title bycluster
#'@description Cluster detection based on maximum potential cluster.
#'@param Lik Likelihood for each potential cluster.
#'#'@param Lambda_dense Initial estimates for relative risk for each potential cluster;\eqn{\Lambda}.
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@param maxclust Maximum number of clusters to be detected.
#'@return Weighted relative risks for identified locations.
bycluster <-  function(Lik, Lambda_dense, sparsemat,maxclust){
    wtMAT <- matrix(rep(NA, maxclust*nrow(sparsemat)), ncol=maxclust)
    wtMAT0 <- matrix(rep(NA, maxclust*nrow(sparsemat)), ncol=maxclust)
    maxpcs <- rep(NA, maxclust)
    stop = FALSE
    ixall <- NULL
        for (i in 1:maxclust){
            Lik[ixall] <-0
            wtmp <- likweights(Lik)
            wtmp[ixall] <-0
            #find potential cluster with largest weight
            maxpc <- which.max(as.vector(wtmp))
            message(paste0("Potential cluster identified: ",(maxpc)))
            maxpcs[i] <- maxpc 
            wtMAT0[,i] <- wtmp
            #find all potential clusters that overlap that PC
            pcmax <- rep(0,length(wtmp)); pcmax[maxpc] <-1; pcmax <- matrix(pcmax,ncol=1)
            pcpcmax <- t(t(sparsemat)%*%pcmax)%*%t(sparsemat); pcpcmax <- ifelse(pcpcmax!=0,1,0)
            ix <- which(pcpcmax!=0)
            wtmp[ix] <- likweights(Lik[ix]) 
            wtMAT[,i] <- wtmp
            if(all(ix %in% ixall)){
                stop = TRUE
                maxi = i
                break
            }
            if (stop){break}
            #save everything that's already been considered
            ixall <- c(ixall,ix)
            maxi = i
        }
    wLambda <- crossprod(wtMAT[,1:(maxi)], Lambda_dense)
    return(list(wLambda = wLambda,
                wtMAT = wtMAT,
                #wtMAT0 = wtMAT0,
                maxpcs = maxpcs))
}


#'@title detectclusters
#'@description Detect disease clusters either by location or by potential cluster.
#'@param sparsemat Large sparsematrix 
#'@param Ex Vector of expected counts.
#'@param Yx Vector of observed counts.
#'@param numCenters Number of centroids.
#'@param Time Number of time periods.
#'@param maxclust Maximum number of clusters allowed. TODO - allow this to be unknown.
#'@param byloc If clusters should be identified by maximum location (\code{TRUE}) or maximum potential cluster (\code{FALSE}). Default is \code{TRUE}.
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"}.
#'@param cv option for cross-validation instead of AIC/BIC. Default is set to FALSE
#'@param overdisp.est Overdispersion parameter estimated across all simulations (max).
#'@return Returns list for each iteration with weighted relative risks by location inside identified cluster.
#'@export
detectclusters <- function(sparsemat, Ex, Yx,numCenters,Time, maxclust,byloc=TRUE, model=c("poisson", "binomial"),overdisp.est) {
    if(is.null(overdisp.est)){
        quasi <- FALSE
        message(paste0("Model specified: ", model))
    } else{ quasi <- TRUE
    message(paste0("Model specified: ", "quasi-",model))
    }
    tsparsemat <- Matrix::t(sparsemat) 
    if (model=="poisson"){
        out <- poisLik(Ex, Yx, tsparsemat)  
    }
    else if (model=="binomial"){
        out <- binomLik(Ex, Yx, tsparsemat) #in dev
    }
    else {
        stop("Model not specified. Please indicate `poisson` or `binomial`.")
    }
    Lik <- out$Lik
    Lambda_dense <- out$Lambda_dense
    if(byloc==FALSE){
        message("Cluster detection by potential cluster")
        res <- bycluster(Lik, Lambda_dense, tsparsemat, maxclust)
        selection <- clusterselect(res[[1]], Yx, Ex, model,maxclust, numCenters, Time, quasi,overdisp.est)
        return(list(wLambda = res[[1]],
                    loglik = selection$loglik,
                    selection.bic = selection$select.bic,
                    selection.aic = selection$select.aic,
                    selection.aicc = selection$select.aicc,
                    selection.bic_forceid = ifelse(selection$select.bic==0,1, selection$select.bic),
                    selection.aic_forceid = ifelse(selection$select.aic==0,1, selection$select.aic),
                    selection.aicc_forceid = ifelse(selection$select.aicc==0,1,selection$select.aicc),
                    wtMAT = res[[2]],
                    maxid = res[[3]],
                    Lambda_dense = Lambda_dense))
    }
    else{
        #default
        message("Cluster detection by location")
        res <- bylocation(Lik, Lambda_dense, tsparsemat, maxclust)
        selection <- clusterselect(res[[1]], Yx, Ex, model,maxclust, numCenters, Time, quasi,overdisp.est)
        return(list(wLambda = res[[1]],
                    loglik = selection$loglik,
                    selection.bic = selection$select.bic,
                    selection.aic = selection$select.aic,
                    selection.aicc = selection$select.aicc,
                    selection.bic_forceid = ifelse(selection$select.bic==0,1,selection$select.bic),
                    selection.aic_forceid = ifelse(selection$select.aic==0,1, selection$select.aic),
                    selection.aicc_forceid = ifelse(selection$select.aic==0,1,selection$select.aicc),
                    wtMAT = res[[2]],
                    maxid = res[[3]],
                    Lambda_dense = Lambda_dense))
    }
}

#'@title clusterselect
#'@description Select the optimal number of (overlapping) clusters using either information criteria or cross-validation
#'@param wLambda Stacked estimates for relative risks in all locations.
#'@param maxclust Maximum number of clusters to consider.
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"}.
#'@param numCenters Number of polygon centroids (total number of locations).
#'@param Time Number of time periods.
#'@param quasi Boolean. \code{TRUE} indicates a quasi-Poisson or quasi-binomial model that accounts for overdispersion. \code{FALSE} indicates a Poisson or binomial model without adjustment for overdispersion or underdispersion.
#'@param overdisp.est Overdispersion (or underdispersion) parameter estimate.
#'@return Returns a list with the loglikelihood, selection by BIC (or QBIC), selection by AIC (or QAIC), and selection by AICc (or QAICc).
clusterselect <- function(wLambda,Yx, Ex, model,maxclust, numCenters, Time,quasi,overdisp.est){
    remtab <- which(sapply(1:nrow(wLambda), function(k) all(is.nan(wLambda[k,])))==FALSE)
    loglik <- sapply(1:length(remtab), function(i) dpoisson(Yx, wLambda[i,], Ex))
    if (model=="poisson"){
        loglik <- c(dpoisson(Yx, rep(1,length(Yx)), Ex), loglik)    
        if (quasi== TRUE){
            message(paste("Overdispersion estimate:", round(overdisp.est,4)))
            if(quasi == TRUE & is.null(overdisp.est)) warning("No overdispersion for quasi-Poisson model. Please check.")
        } else {
            message("no overdispersion being estimated")
        }
    }
    #find K clusters
    K <- seq(from=0, to=nrow(wLambda))
    #calc
    if (quasi== TRUE){
        PLL.bic <- (-2*loglik/overdisp.est) + (K+1)*log(sum(Yx)) 
        PLL.aic <- 2*(K+1) - 2*(loglik/overdisp.est)
        PLL.aicc <- 2*(K+1) - 2*(loglik/overdisp.est) +
            ((2*(K+1)*(K+1 + 1))/(sum(Yx) - K+1 - 1))
    } else {
        PLL.bic <- (-2*loglik) + K*log(sum(Yx)) 
        PLL.aic <- 2*K - 2*loglik
        PLL.aicc <- 2*(K) - 2*(loglik) +
            ((2*K*(K + 1))/(sum(Yx) - K - 1))
    }
    if(any(is.infinite(PLL.aicc))) {
        idx <- which(is.infinite(PLL.aicc))
        PLL.aicc[idx:length(PLL.aicc)] <- NA
    }
    #select
    select.bic <- which.min(PLL.bic)-1 #-1 because the first element corresponds to 0 clusters (null)
    select.aic <- which.min(PLL.aic)-1
    select.aicc <- which.min(PLL.aicc)-1 
    cat(paste0("\t","Num. selected clusters by BIC: ", (select.bic),"\n",
               "\t","Num. selected clusters by AIC: ",(select.aic),"\n",
               "\t","Num. selected clusters by AICc: ",(select.aicc)))
    return(list(loglik = loglik,
                select.bic = select.bic,
                select.aic = select.aic,
                select.aicc = select.aicc))
    
}


####################################################################
#CLUSTACKBOUNDS FUNCTIONS
####################################################################


#'@title bucklandbounds.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on Buckland et al. for specific cells (given by \code{cellsix_out}) 
#'@param thetaa Stacked relative risk estimate for the cluster(s)
#'@param thetai Individual relative risk estimates for each potential cluster.
#'@param res Resultant object from \code{detectclusters()}.
#'@param w_q Likelihood-weights for all Q<K potential clusters.
#'@param outExp_out Expected counts for each potential cluster.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC)
#'@param transform If a transformation is used. Default is \code{"none"}. Other transformation currently available is \code{"log"}.
#'@param tsparsemat Transpose of large sparsematrix where columns are potential clusters and rows are space-time locations.
#'@param overdisp.est Estimate of overdispersion (or underdispersion).
#'@param cellsix_out Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate for each of the cells in \code{cellsix_out}
#'@export
#'
bucklandbounds.cells <- function(thetaa, res, w_q, outExp_out ,IC,transform="none",tsparsemat,overdisp.est, cellsix_out, conf.level) {
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
            varthetai <- sapply(1:nrow(tsparsemat), function(k) overdisp.est*(1/(thetai[k,]*outExp_out )))
        } else {
            varthetai <- sapply(1:nrow(tsparsemat), function(k) 1/(thetai[k,]*outExp_out))
        }
        withintheta <- sapply(1:length(cellsix_out), function(j) (log(thetai[,cellsix_out][,j]) - log(thetaa[j])))^2
        wtnbtn <- sapply(1:nrow(tsparsemat), function(k) sqrt(varthetai[cellsix_out,k] + withintheta[k,]))
        varthetas_w <- matrix(w_q, nrow = 1)%*%t(wtnbtn)
        var_thetaa <- (as.vector(varthetas_w))^2
        UBa = exp(as.vector(log(thetaa)) + critval*sqrt(var_thetaa))
        LBa = exp(as.vector(log(thetaa)) - critval*sqrt(var_thetaa))
        
    } else {
      #browser()
        if(!is.null(overdisp.est)){
            varthetai <- sapply(1:nrow(tsparsemat), function(k) overdisp.est*thetai[k,]/outExp_out)
        } else {
            varthetai <- sapply(1:nrow(tsparsemat), function(k) thetai[k,]/outExp_out)
        }
        withintheta <- sapply(1:length(cellsix_out), function(j) (thetai[,cellsix_out][,j] - thetaa[j]))^2
        wtnbtn <- sapply(1:nrow(tsparsemat), function(k) sqrt(varthetai[cellsix_out,k] + withintheta[k,]))
        varthetas_w <- matrix(w_q, nrow = 1)%*%t(wtnbtn)
        var_thetaa <- (as.vector(varthetas_w))^2
        UBa = as.vector(thetaa) + critval*sqrt(var_thetaa)
        LBa = as.vector(thetaa) - critval*sqrt(var_thetaa)
    }
    return(list(buckland.LB = LBa,
                clusterMA = thetaa,
                buckland.UB = UBa))
}


#'@title maw2.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on Buckland et al. for specific cells (given by \code{cellsix_out}) 
#'@param thetaa Stacked relative risk estimate for the cluster(s)
#'@param thetai Individual relative risk estimates for each potential cluster.
#'@param res Resultant object from \code{detectclusters()}.
#'@param w_q Likelihood-weights for all Q<K potential clusters.
#'@param outExp_out Expected counts for each potential cluster.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC).
#'@param transform If a transformation is used. Default is \code{"none"}. Other transformation currently available is \code{"log"}.
#'@param tsparsemat Transpose of large sparsematrix where columns are potential clusters and rows are space-time locations.
#'@param overdisp.est Estimate of overdispersion (or underdispersion).
#'@param cellsix_out Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate for each of the cells in \code{cellsix_out}
#'@export
#'
maw2.cells <- function(thetaa,res, w_q, outExp, IC, transform,tsparsemat, overdisp.est, cellsix_out, conf.level=0.95){
    thetai <- res$Lambda_dense
    critval <- qnorm(1-(1-conf.level)/2)
    if(IC=="aic"){
        thetaa <- res$wLambda[res$selection.aic,][cellsix_out]
    } else {
        thetaa <- res$wLambda[res$selection.bic,][cellsix_out]
    }
        if(transform=="log"){
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(tsparsemat), function(k) overdisp.est*(1/(thetai[k,]*outExp)))
            } else {
                varthetai <- sapply(1:nrow(tsparsemat), function(k) 1/(thetai[k,]*outExp))
            }
            withintheta <- sapply(1:length(cellsix_out), function(j) (log(thetai[, cellsix_out][,j]) - log(thetaa[j])))^2
            wtnbtn <- sapply(1:nrow(tsparsemat), function(k) (varthetai[cellsix_out,k] + withintheta[k,]))
            #var_thetaa <- sum(w_q*(varthetai + withintheta))
            var_thetaa <- as.vector(matrix(w_q, nrow=1)%*%t(wtnbtn))
            UBa = exp(as.vector(log(thetaa)) + critval*sqrt(var_thetaa))
            LBa = exp(as.vector(log(thetaa)) - critval*sqrt(var_thetaa))
        } else{
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(tsparsemat), function(k) overdisp.est*thetai[k,]/outExp)
            } else {
                varthetai <- sapply(1:nrow(tsparsemat), function(k) thetai[k,]/outExp)
            }
            withintheta <- sapply(1:length(cellsix_out), function(j) (thetai[,cellsix_out][,j] - thetaa[j]))^2
            wtnbtn <- sapply(1:nrow(tsparsemat), function(k) sqrt(varthetai[cellsix_out,k] + withintheta[k,]))
            #var_thetaa <- sum(w_q*(varthetai + withintheta))
            var_thetaa <- as.vector(matrix(w_q, nrow=1)%*%t(wtnbtn))
            UBa = as.vector(thetaa) + critval*sqrt(var_thetaa)
            LBa = as.vector(thetaa) - critval*sqrt(var_thetaa)
        }
    return(list(maw2.LB = LBa,
                clusterMA = thetaa,
                maw2.UB = UBa))
}

#' @title mata_tailareazcore
#' @description Calculate MATA tail area weighted z-score based on Turek et al.
#' @param thetai Each individual cluster relative risk.
#' @param thetaa Stacked relative risk estimate.
#' @param sd.thetai Standard deviation of estimated variance for each individual cluster relative risk.
#' @param w_q Likelihood-weights for all Q<K potential clusters.
#' @param alpha \eqn{\alpha}, probability of rejecting the null hypothesis when it is true.
mata_tailareazscore <- function(thetai, thetaa, sd.thetai, w_q, alpha){
    thetai <- as.vector(thetai)
    zval <- (thetaa - thetai)/sd.thetai
    zpnorm <- pnorm(zval)
    w_zpnorm <- sum((w_q*zpnorm))-alpha
}


#'@title matabounds.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on MATA (model-averaged tail area) intervals based on Turek et al. for specific cells (given by \code{cellsix_out}).
#'@param thetaa Stacked relative risk estimate for the cluster(s).
#'@param res Resultant object from \code{detectclusters()}.
#'@param w_q Likelihood-weights for all Q<K potential clusters.
#'@param outExp Expected counts.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC).
#'@param transform If a transformation is used. Default is \code{"none"}. Other transformation currently available are \code{"log"} and \code{"sqrt"}.
#'@param tsparsemat Transpose of large sparsematrix where columns are potential clusters and rows are space-time locations.
#'@param overdisp.est Estimate of overdispersion (or underdispersion).
#'@param cellsix_out Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate
#'@export
#'
matabounds.cells <- function(thetaa,res, w_q, outExp, IC, transform=c("none","log", "sqrt"),tsparsemat,overdisp.est, cellsix_out, conf.level=0.95) {
    if(is.null(transform)){
        none = matabounds_none.cells(thetaa, res, w_q, outExp, IC, tsparsemat, overdisp.est, cellsix_out, conf.level)
    }
    switch(transform,
           none = matabounds_none.cells(thetaa, res, w_q, outExp, IC, tsparsemat, overdisp.est, cellsix_out, conf.level),
           log = matabounds_log.cells(thetaa, res, w_q, outExp, IC, tsparsemat, overdisp.est, cellsix_out, conf.level),
           sqrt = matabounds_sqrt.cells(thetaa, res, w_q, outExp, IC, tsparsemat, overdisp.est, cellsix_out, conf.level))
}



#'@title matabounds_none.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on MATA (model-averaged tail area) intervals based on Turek et al. for specific cells (given by \code{cellsix_out}) with no transformation(\code{transform="none"})..
#'@param thetaa Stacked relative risk estimate for the cluster(s).
#'@param res Resultant object from \code{detectclusters()}.
#'@param w_q Likelihood-weights for all Q<K potential clusters.
#'@param outExp Expected counts.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC).
#'@param tsparsemat Transpose of large sparsematrix where columns are potential clusters and rows are space-time locations.
#'@param overdisp.est Estimate of overdispersion (or underdispersion).
#'@param cellsix_out Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate
#'@export
#'
matabounds_none.cells <- function(thetaa, res, w_q,outExp, IC, tsparsemat, overdisp.est,cellsix_out, conf.level=0.95) {
    thetai <- res$Lambda_dense
    alpha <- (1-conf.level)/2
    if(IC=="aic"){
        thetaa <- res$wLambda[res$selection.aic,][cellsix_out]
    } else {
        thetaa <- res$wLambda[res$selection.bic,][cellsix_out]
    }
    if(!is.null(overdisp.est)){
        varthetai <- sapply(1:nrow(tsparsemat), function(k) overdisp.est*(1/(thetai[k,]*outExp)))
    } else {
        varthetai <- sapply(1:nrow(tsparsemat), function(k) thetai[k,]/outExp)
    }
    mataLB <- sapply(1:length(cellsix_out),
                     function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                         thetai= thetai[,cellsix_out[k]],
                                         sd.thetai=sqrt(varthetai[cellsix_out[k],]),
                                         w_q=w_q, alpha=alpha, tol=1e-8)$root)
    
    mataUB <- sapply(1:length(cellsix_out),
                     function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                         thetai= thetai[,cellsix_out[k]],
                                         sd.thetai=sqrt(varthetai[cellsix_out[k],]),
                                         w_q=w_q, alpha=1-alpha, tol=1e-8)$root)
    
    return(list(mata.LB = mataLB,
                clusterMA = thetaa,
                mata.UB = mataUB))
}



#'@title matabounds_sqrt.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on MATA (model-averaged tail area) intervals based on Turek et al. for specific cells (given by \code{cellsix_out}) with square root transformation(\code{transform="sqrt"}).
#'@param thetaa Stacked relative risk estimate for the cluster(s).
#'@param res Resultant object from \code{detectclusters()}.
#'@param w_q Likelihood-weights for all Q<K potential clusters.
#'@param outExp Expected counts.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC).
#'@param tsparsemat Transpose of large sparsematrix where columns are potential clusters and rows are space-time locations.
#'@param overdisp.est Estimate of overdispersion (or underdispersion).
#'@param cellsix_out Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate
#'@export
#'
matabounds_sqrt.cells <- function(thetaa, res, w_q,outExp, IC, tsparsemat, overdisp.est,cellsix_out, conf.level=0.95) {
    thetai <- res$Lambda_dens
    alpha <- (1-conf.level)/2
    if(IC=="aic"){
        thetaa <- res$wLambda[res$selection.aic,][cellsix_out]
    } else {
        thetaa <- res$wLambda[res$selection.bic,][cellsix_out]
    }
    if(!is.null(overdisp.est)){
        varthetai <- overdisp.est*(1/(4*outExp))
    } else {
        Tvarthetai <- 1/(4*outExp)
    }
        mataLB <- sapply(1:length(cellsix_out), 
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                             thetai= sqrt(thetai[,cellsix_out[k]]),
                                             sd.thetai=sqrt(Tvarthetai[k]),
                                             w_q=w_q, alpha=conf.level, tol=1e-8)$root)
        
        mataUB <- sapply(1:length(cellsix_out),
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                             thetai= sqrt(thetai[,cellsix_out[k]]),
                                             sd.thetai=sqrt(Tvarthetai[k]),
                                             w_q=w_q, alpha=1-conf.level, tol=1e-8)$root)
    return(list(matasqrt.LB = (mataLB)^2,
                clusterMA = thetaa,
                matasqrt.UB = (mataUB)^2))
}


#'@title matabounds_log.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates based on MATA (model-averaged tail area) intervals based on Turek et al. for specific cells (given by \code{cellsix_out}) with natural log transformation(\code{transform="log"}).
#'@param thetaa Stacked relative risk estimate for the cluster(s).
#'@param res Resultant object from \code{detectclusters()}.
#'@param w_q Likelihood-weights for all Q<K potential clusters.
#'@param outExp Expected counts.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC).
#'@param tsparsemat Transpose of large sparsematrix where columns are potential clusters and rows are space-time locations.
#'@param overdisp.est Estimate of overdispersion (or underdispersion).
#'@param cellsix_out Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate
#'@export
#'
matabounds_log.cells<- function(thetaa, res, w_q,outExp, IC, tsparsemat, overdisp.est,cellsix_out, conf.level=0.95) {
    thetai <- res$Lambda_dense
    alpha <- (1-conf.level)/2
    if(IC=="aic"){
        thetaa <- res$wLambda[res$selection.aic,][cellsix_out]
    } else {
        thetaa <- res$wLambda[res$selection.bic,][cellsix_out]
    }
    if(!is.null(overdisp.est)){
        logTvarthetai <- sapply(1:nrow(tsparsemat), function(k) overdisp.est*(1/(thetai[k,]*outExp)))
    } else{
        logTvarthetai <- sapply(1:nrow(tsparsemat), function(k) 1/(thetai[k,]*outExp))
    }
    

    mataLB <- sapply(1:length(cellsix_out),
                     function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                         thetai= log(thetai[,cellsix_out[k]]),
                                         sd.thetai=sqrt(logTvarthetai[cellsix_out[k],]),
                                         w_q=w_q, alpha=alpha, tol=1e-8)$root)
    mataUB <- sapply(1:length(cellsix_out),
                     function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                                         thetai= log(thetai[,cellsix_out[k]]),
                                         sd.thetai=sqrt(logTvarthetai[cellsix_out[k],]),
                                         w_q=w_q, alpha=1-alpha, tol=1e-8)$root)
    return(list(matalog.LB = exp(mataLB),
                clusterMA = thetaa,
                matalog.UB = exp(mataUB)))
}




#'@title nonma.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates assuming the cluster was known \textit{a priori}  for specific cells (given by \code{cellsix_out}).
#'@param thetaa Stacked relative risk estimate for the cluster(s)
#'@param res Resultant object from \code{detectclusters()}.
#'@param w Likelihood-weights for all potential clusters.
#'@param id_ic The number of clusters identified either by BIC (QBIC) or AIC (QAIC).
#'@param outExp_out Expected counts for each potential cluster.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC). This should match \code{id_ic} argument.
#'@param transform If a transformation is used. Default is \code{"none"}. Other transformation currently available is \code{"log"}.
#'@param byloc Boolean. Use \code{TRUE} when stacking by location. Use \code{FALSE} when stacking by potential cluster.
#'@param cellrisk_wt_out Stacked relative risk estimates for cells indicted by \code{cellsix_out}.
#'@param cellsix_out Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate
#'@export
#'
nonma.cells <- function(thetaa, thetai,res,w,id_ic, outExp_out,IC, transform, byloc=TRUE, cellrisk_wt_out=NULL, cellsix_out, conf.level=0.95){
    critval <- qnorm(1-(1-conf.level)/2)
    if(byloc==FALSE){
        if(transform=="log"){
            if(IC=="aic") {
                thetai <- res$Lambda_dense[res$maxid[res$selection.aic],][cellsix_out]
                se_thetai <- sqrt(1/(thetai*outExp_out))
                clusterMA <- res$wLambda[res$selection.aic,][cellsix_out]
                
            } else {
                thetai <- res$Lambda_dense[res$maxid[res$selection.bic],][cellsix_out]
                se_thetai <- sqrt(1/(thetai*outExp_out))
                clusterMA <- res$wLambda[res$selection.bic,][cellsix_out]
            }
            nonma.theta.time <- system.time(nonma.theta <- cbind(lb=exp(log(clusterMA)-critval*se_thetai), 
                                                                 clusterMA = clusterMA,
                                                                 ub=exp(log(clusterMA)+critval*se_thetai)))
        }else {
            if(IC=="aic") {
                thetai <- res$Lambda_dense[res$maxid[res$selection.aic],][cellsix_out]
                se_thetai <- sqrt(thetai/outExp_out)
                clusterMA <- res$wLambda[res$selection.aic,][cellsix_out]
            } else {
                thetai <- res$Lambda_dense[res$maxid[res$selection.bic],][cellsix_out]
                se_thetai <- sqrt(thetai/outExp_out)
                clusterMA <- res$wLambda[res$selection.bic,][cellsix_out]
            }
            nonma.theta.time <- system.time(nonma.theta <- cbind(lb=clusterMA-critval*se_thetai, 
                                                                 clusterMA = clusterMA,
                                                                 ub=clusterMA+critval*se_thetai))
        } 
        

    }else {
        if(transform=="log"){
            if(IC=="aic") {
                se_thetai <- sqrt(1/(cellrisk_wt_out*outExp_out))
            } else {
                se_thetai <- sqrt(1/(cellrisk_wt_out*outExp_out))
            }
            nonma.theta.time <- system.time(nonma.theta <- cbind(lb=exp(log(cellrisk_wt_out)-critval*se_thetai), 
                                                                 clusterMA = cellrisk_wt_out,
                                                                 ub=exp(log(cellrisk_wt_out) +critval*se_thetai)))
        } else {
            if(IC=="aic") {
                se_thetai <- sqrt(cellrisk_wt_out/outExp_out)
            } else {
                se_thetai <- sqrt(cellrisk_wt_out/outExp_out)
            }
            nonma.theta.time <- system.time(nonma.theta <- cbind(lb=cellrisk_wt_out-critval*se_thetai, 
                                                                 clusterMA = cellrisk_wt_out,
                                                                 ub=cellrisk_wt_out +critval*se_thetai))
        }
    }
    return(list(nonma.theta.time = nonma.theta.time[[3]],
                nonma.theta = nonma.theta))
}






#'@title nonma_asymp.cells
#'@description Calculate confidence bounds for stacked cluster relative risk estimates assuming the cluster was known \textit{a priori} and the sampling distribution is asymptotic for specific cells (given by \code{cellsix_out}).
#'@param thetaa Stacked relative risk estimate for the cluster(s)
#'@param thetai Individual relative risk estimates for each potential cluster.
#'@param res Resultant object from \code{detectclusters()}.
#'@param w Likelihood-weights for all potential clusters.
#'@param id_ic The number of clusters identified either by BIC (QBIC) or AIC (QAIC).
#'@param outExp_out Expected counts for each potential cluster.
#'@param IC Information criterion used. Currently available for BIC (QBIC) and AIC (QAIC). This should match \code{id_ic} argument.
#'@param transform If a transformation is used. Default is \code{"none"}. Other transformation currently available is \code{"log"}.
#'@param byloc Boolean. Use \code{TRUE} when stacking by location. Use \code{FALSE} when stacking by potential cluster.
#'@param cellrisk_wt_out Stacked relative risk estimates for cells indicted by \code{cellsix_out}.
#'@param cellsix_out Indices of the cells to calculate bounds for.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return List: lower bound estimate, stacked cluster relative risk estimate, upper bound estimate.
#'@export
#'
nonma_asymp.cells <-function(thetaa, thetai,res, w,id_ic, outObs_out,IC, transform, byloc=TRUE, cellrisk_wt_out, cellsix_out, conf.level=0.95){
    critval <- qnorm(1-(1-conf.level)/2)
    if(byloc==FALSE){
        #print("by PC")
        if(transform=="log"){
            if(IC=="aic") {
                thetai <- res$Lambda_dense[res$maxid[res$selection.aic],][cellsix_out]
                se_thetai <- sqrt(1/(thetai*outObs_out))
                clusterMA <- res$wLambda[res$selection.aic,][cellsix_out]
                
            } else {
                thetai <- res$Lambda_dense[res$maxid[res$selection.bic],][cellsix_out]
                se_thetai <- sqrt(1/(thetai*outObs_out))
                clusterMA <- res$wLambda[res$selection.bic,][cellsix_out]
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- cbind(lb=exp(log(clusterMA)-critval*se_thetai), 
                                                                 clusterMA = clusterMA,
                                                                 ub=exp(log(clusterMA)+critval*se_thetai)))
        }else {
            if(IC=="aic") {
                thetai <- res$Lambda_dense[res$maxid[res$selection.aic],][cellsix_out]
                se_thetai <- sqrt(thetai/outObs_out)
                clusterMA <- res$wLambda[res$selection.aic,][cellsix_out]
            } else {
                thetai <- res$Lambda_dense[res$maxid[res$selection.bic],][cellsix_out]
                se_thetai <- sqrt(thetai/outObs_out)
                clusterMA <- res$wLambda[res$selection.bic,][cellsix_out]
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- cbind(lb=clusterMA-critval*se_thetai, 
                                                                 clusterMA = clusterMA,
                                                                 ub=clusterMA+critval*se_thetai))
        } 
    }else {
       # print("by loc")
        #print(cellrisk_wt_out)
        if(transform=="log"){
            if(IC=="aic") {
                se_thetai <- sqrt(1/(cellrisk_wt_out*outObs_out))
            } else {
                se_thetai <- sqrt(1/(cellrisk_wt_out*outObs_out))
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- cbind(lb=exp(log(cellrisk_wt_out)-critval*se_thetai), 
                                                                 clusterMA = cellrisk_wt_out,
                                                                 ub=exp(log(cellrisk_wt_out) +critval*se_thetai)))
        } else {
            if(IC=="aic") {
                se_thetai <- sqrt(cellrisk_wt_out/outObs_out)
            } else {
                se_thetai <- sqrt(cellrisk_wt_out/outObs_out)
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- cbind(lb=cellrisk_wt_out-critval*se_thetai, 
                                                                 clusterMA = cellrisk_wt_out,
                                                                 ub=cellrisk_wt_out +critval*se_thetai))
        }
    }
    return(list(nonma_asymp.theta.time = nonma_asymp.theta.time[[3]],
                nonma_asymp.theta = nonma_asymp.theta))
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
#'@param Obs Observed counts for each cell.
#'@param cellsix_out Indices of the cells to calculate bounds for. 
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@details The confidence bounds methods returned are:
#'\itemize{
#'\item \code{nonma}: Non-stacked (assumes the cluster was known \textit{a priori}).
#'\item \code{nonmaTlog}: Non-stacked  (assumes the cluster was known \textit{a priori}); natural log-transformed.
#'\item \code{nonma_asymp}: Non-stacked (assumes the cluster was known \textit{a priori}) assuming asymptotic sampling distribution.
#'\item \code{nonma_asympTlog}: Non-stacked (assumes the cluster was known \textit{a priori}) assuming asymptotic sampling distribution; natural log-transformed.
#'\item \code{buck}: Stacked confidence bounds based on Buckland et al.
#'\item \code{buckTlog}: Stacked confidence bounds based on Buckland et al.; natural log-transformed.
#'\item \code{maw}: Model-average weighted confidence bounds based on Burnham and Anderson
#'\item \code{mawTlog}: Model-average weighted confidence bounds based on Burnham and Anderson.
#'\item \code{mata}: Model-average tail area (MATA) intervals based on Turek et al.
#'\item \code{mataTlog}: Model-average tail area (MATA) intervals based on Turek et al.; natural log-transformed.
#'
#'}
#calcbounds <- function(id_ic, IC, res, byloc, Ex, Obs,target=c("cluster", "cells"), cellsix=NULL, sparsemat, conf.level=0.95){
calcbounds <- function(id_ic, IC, res, byloc, Ex, Obs, cellsix=NULL, sparsemat, conf.level=0.95){
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
      out[[i]] <- calcbounds.cells(id_ic, IC, res, byloc, Ex, Obs,as.vector(w[,i]), thetaa,thetai, sparsemat, cellsix, conf.level)
      
    }
  } else {
    out <- vector(mode = "list", length = res$selection.bic)
    w <- matrix(res$wtMAT[,1:res$selection.bic], ncol=res$selection.bic)
    for(i in 1:res$selection.bic){
      thetaa <- sum(thetai*w[,i])
      out[[i]] <-calcbounds.cells(id_ic, IC, res, byloc, Ex, Obs,as.vector(w[,i]), thetaa,thetai, sparsemat, cellsix, conf.level)
      
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
#'@param Obs Observed counts for each cell.
#'@param w Likelihood-weights for all potential clusters.
#'@param thetaa Stacked relative risk estimate for the cluster(s)
#'@param thetai Individual relative risk estimates for each potential cluster.
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@param cellsix Indices of the cells to calculate bounds for. 
#'@param conf.level Confidence level for the interval. Default is 0.95. 
#'@return Returns large list of confidence bounds and stacked estimates in addition to timings for each of the confidence bounds methods.    
calcbounds.cells <- function(id_ic, IC, res, byloc, Ex, Obs,w, thetaa,thetai, sparsemat, cellsix, conf.level=0.95){
    if(id_ic==0){
        emptynonma <- vector(mode = "list", length = 2)
        emptynonma[[1]] <- 0
        emptynonma[[2]] <- matrix(rep(0, length(cellsix)*3), ncol=3)
        outnonma <- outnonmaTlog <- outnonma_asymp <- outnonma_asympTlog <- emptynonma
        names(outnonma) <- c("nonma.theta.time", "nonma.theta")
        names(outnonmaTlog) <- c("nonma.theta.time", "nonma.theta")
        names(outnonma_asymp) <- c("nonma_asymp.theta.time", "nonma_asymp.theta")
        names(outnonma_asympTlog) <- c("nonma_asymp.theta.time", "nonma_asymp.theta")
        return(list(
            outnonma = outnonma,
            outnonmaTlog = outnonmaTlog,
            outnonma_asymp = outnonma_asymp,
            outnonma_asympTlog = outnonma_asympTlog,
            outbuck.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            outbuckTlog.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            outmaw2.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            outmaw2Tlog.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            outmata.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            outmataTlog.theta = list(as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix))), as.vector(rep(0, length(cellsix)))),
            
            outnonma.time = 000,
            outnonmaTlog.time = 000,
            outnonma_asymp.time = 000,
            outnonma_asympTlog.time = 000,
            outbuck.theta.time = 000,
            outbuckTlog.theta.time = 000,
            outmaw2.theta.time = 000,
            outmaw2Tlog.theta.time = 000,
            outmata.theta.time = 000,
            outmataTlog.theta.time = 000
        ))
        
        
    } else {
        cellrisk_wt_out <- rep(NA, length(cellsix))
        cellsix_out <- cellsix
        if(byloc==TRUE){
            #print(paste0("byloc", byloc))
            outExp <- t(sparsemat)%*%Ex
            outObs <- t(sparsemat)%*%Obs
            outExp_out <- Ex[cellsix]
            outObs_out <- Obs[cellsix]
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
            outExp <- t(sparsemat)%*%Ex
            outObs <- t(sparsemat)%*%Obs
            #print(paste0("cellsix: ", cellsix))
            outExp_out <- outExp@x[cellsix]
            outObs_out <- outObs@x[cellsix]
            cellrisk_wt_out <- NULL
        }    
        
        outnonma.time <- system.time(outnonma <- nonma.cells(thetaa, thetai,res, w, id_ic,
                                                             outExp_out, IC=IC, transform="none", byloc, 
                                                             cellrisk_wt_out, cellsix_out,
                                                             conf.level))
        outnonmaTlog.time <- system.time(outnonmaTlog <- nonma.cells(thetaa, thetai,res, w, id_ic,
                                                                     outExp_out, IC=IC, transform="log", byloc, 
                                                                     cellrisk_wt_out, cellsix_out,
                                                                     conf.level))
        
        outnonma_asymp.time <- system.time(outnonma_asymp <- nonma_asymp.cells(thetaa, thetai,res, w, id_ic,
                                                                               outObs_out, IC=IC, transform="none", byloc, cellrisk_wt_out, cellsix_out))
        outnonma_asympTlog.time <- system.time(outnonma_asympTlog <- nonma_asymp.cells(thetaa, thetai, res, w, id_ic,
                                                                                       outObs_out, IC=IC, transform="log", byloc, 
                                                                                       cellrisk_wt_out, cellsix_out, conf.level))
        message("Non-model averaged bounds finished")
        
        outbuck.theta.time <- system.time(outbuck.theta <- bucklandbounds.cells(thetaa,
                                                                                res,
                                                                                w_q=w,
                                                                                Ex,
                                                                                IC=IC,
                                                                                transform="none",
                                                                                tsparsemat=t(sparsemat),
                                                                                overdisp.est, cellsix_out,
                                                                                conf.level))
        outbuckTlog.theta.time <- system.time(outbuckTlog.theta <- bucklandbounds.cells(thetaa,
                                                                                        res,
                                                                                        w_q=w,
                                                                                        Ex,
                                                                                        IC=IC,
                                                                                        transform="log",
                                                                                        tsparsemat=t(sparsemat),
                                                                                        overdisp.est, cellsix_out,
                                                                                        conf.level))
        message("Buckland bounds finished")
        outmaw2.theta.time <- system.time(outmaw2.theta <- maw2.cells(thetaa,
                                                                      res,
                                                                      w_q=w,
                                                                      Ex,
                                                                      IC,
                                                                      transform="none",
                                                                      tsparsemat=t(sparsemat),
                                                                      overdisp.est,
                                                                      cellsix_out,
                                                                      conf.level))
        outmaw2Tlog.theta.time <- system.time(outmaw2Tlog.theta <- maw2.cells(thetaa,
                                                                              res,
                                                                              w_q=w,
                                                                              Ex,
                                                                              IC,
                                                                              transform="log",
                                                                              tsparsemat=t(sparsemat),
                                                                              overdisp.est,
                                                                              cellsix_out,
                                                                              conf.level))
        
        message("Burnham & Anderson bounds finished")
        outmata.theta.time <- system.time(outmata.theta <- matabounds.cells(thetaa,
                                                                            res,
                                                                            w_q=w,
                                                                            Ex,
                                                                            IC,
                                                                            transform="none",
                                                                            tsparsemat=t(sparsemat),
                                                                            overdisp.est,
                                                                            cellsix_out,
                                                                            conf.level))
        outmataTlog.theta.time <- system.time(outmataTlog.theta <- matabounds.cells(thetaa,
                                                                                    res,
                                                                                    w_q=w,
                                                                                    Ex,
                                                                                    IC,
                                                                                    transform="log",
                                                                                    tsparsemat=t(sparsemat),
                                                                                    overdisp.est,
                                                                                    cellsix_out,
                                                                                    conf.level))
        message("MATA bounds finished")
        
        
        return(list(
            outnonma = outnonma,
            outnonmaTlog = outnonmaTlog,
            outnonma_asymp = outnonma_asymp,
            outnonma_asympTlog = outnonma_asympTlog,
            outbuck.theta = outbuck.theta,
            outbuckTlog.theta = outbuckTlog.theta,
            outmaw2.theta = outmaw2.theta,
            outmaw2Tlog.theta = outmaw2Tlog.theta,
            outmata.theta = outmata.theta,
            outmataTlog.theta = outmataTlog.theta,
            
            outnonma.time = outnonma.time[[3]],
            outnonmaTlog.time = outnonmaTlog.time[[3]],
            outnonma_asymp.time = outnonma_asymp.time[[3]],
            outnonma_asympTlog.time = outnonma_asympTlog.time[[3]],
            outbuck.theta.time = outbuck.theta.time[[3]],
            outbuckTlog.theta.time = outbuckTlog.theta.time[[3]],
            outmaw2.theta.time = outmaw2.theta.time[[3]],
            outmaw2Tlog.theta.time = outmaw2Tlog.theta.time[[3]],
            outmata.theta.time = outmata.theta.time[[3]],
            outmataTlog.theta.time = outmataTlog.theta.time[[3]]
        ))
    }
    
}

