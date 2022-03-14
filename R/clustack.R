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
#'@description Poisson-based LRT statistic comparing cluster model to null model.
#'@param Ex Vector of expected counts.
#'@param Yx Vector of observed counts.
#'@param sparsemat Large sparse matrix where rows are potential clusters and columns are space-time locations.
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



#'@title bylocation
#'@description Cluster detection based on maximum location.
#'@param Lik Likelihood for each potential cluster.
#'@param Lambda_dense Initial estimates for relative risk for each potential cluster;\eqn{\Lambda}.
#'@param sparsemat Large sparse matrix where rows are potential clusters and columns are space-time locations.
#'@param maxclust Maximum number of clusters to be detected.
#'@return List. First element is a Weighted relative risks for identified locations; second element is a large matrix of weights for each iteration; third element are the maximum location IDs.
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
        pclocmax <- as.vector(t(locmax)%*%Matrix::t(sparsemat))
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
                maxlocs = maxlocs))
}




#'@title bycluster
#'@description Cluster detection based on maximum potential cluster.
#'@param Lik Likelihood for each potential cluster.
#'#'@param Lambda_dense Initial estimates for relative risk for each potential cluster;\eqn{\Lambda}.
#'@param sparsemat Large sparse matrix where rows are potential clusters and columns are space-time locations.
#'@param maxclust Maximum number of clusters to be detected.
#'@return List. First element is a Weighted relative risks for identified locations; second element is a large matrix of weights for each iteration; third element are the maximum potential cluster IDs.
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
            pcpcmax <- Matrix::t(Matrix::t(sparsemat)%*%pcmax)%*%Matrix::t(sparsemat); pcpcmax <- ifelse(pcpcmax!=0,1,0)
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
                maxpcs = maxpcs))
}


#'@title detectclusters
#'@description Detect disease clusters either by location or by potential cluster.
#'@param sparsemat Large sparse matrix of potential clusters.
#'@param Ex Vector of expected counts.
#'@param Yx Vector of observed counts.
#'@param numCenters Number of centroids.
#'@param Time Number of time periods.
#'@param maxclust Maximum number of clusters allowed. TODO - allow this to be unknown.
#'@param byloc If clusters should be identified by maximum location (\code{TRUE}) or maximum potential cluster (\code{FALSE}). Default is \code{FALSE} (detection by potential cluster).
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"}.
#'@param overdisp.est Overdispersion estimate.
#'@return Returns a large list.
#'@details The elements of the returned list: 1) the weighted relative risks for each cell; 2) the LRT statistic; 3) the selected number of clusters by (Q)BIC; 4) the selected number of clusters by (Q)AIC; 5) the selected number of clusters by (Q)AICc; 6) the selected number of clusters by (Q)BIC when forced to detect a cluster; 7) the selected number of clusters by (Q)AIC when forced to detect a cluster; 8) the selected number of clusters by (Q)AICc when forced to detect a cluster; 9) matrix of weights; 10) maximum potential cluster or location IDs identified; 11) large matrix of single-cluster relative risk estimates.
#'@export
detectclusters <- function(sparsemat, Ex, Yx,numCenters,Time, maxclust,byloc=FALSE, model=c("poisson", "binomial"),overdisp.est) {
    #browser()
    if(is.null(overdisp.est)){
        isoverdisp <- FALSE
        message(paste0("Model specified: ", model))
    } else{ isoverdisp <- TRUE
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
        selection <- clusterselect(res[[1]], Yx, Ex, model,maxclust, numCenters, Time, isoverdisp,overdisp.est)
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
        message("Cluster detection by location")
        res <- bylocation(Lik, Lambda_dense, tsparsemat, maxclust)
        selection <- clusterselect(res[[1]], Yx, Ex, model,maxclust, numCenters, Time, isoverdisp,overdisp.est)
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
#'@description Select the optimal number of (overlapping) clusters using either information criteria.
#'@param wLambda Stacked estimates for relative risks in all locations.
#'@param Yx Vector of observed counts.
#'@param Ex Vector of expected counts.
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"}.
#'@param maxclust Maximum number of clusters to consider.
#'@param numCenters Number of polygon centroids (total number of locations).
#'@param Time Number of time periods.
#'@param isoverdisp Boolean. \code{TRUE} indicates a quasi-Poisson or quasi-binomial model that accounts for overdispersion. \code{FALSE} indicates a Poisson or binomial model without adjustment for overdispersion or underdispersion.
#'@param overdisp.est Overdispersion (or underdispersion) parameter estimate.
#'@return Returns a list with the loglikelihood, selection by (Q)BIC, selection by (Q)AIC, and selection by (Q)AICc.
clusterselect <- function(wLambda,Yx, Ex, model,maxclust, numCenters, Time,isoverdisp,overdisp.est){
    remtab <- which(sapply(1:nrow(wLambda), function(k) all(is.nan(wLambda[k,])))==FALSE)
    loglik <- sapply(1:length(remtab), function(i) dpoisson(Yx, wLambda[i,], Ex))
    if (model=="poisson"){
        loglik <- c(dpoisson(Yx, rep(1,length(Yx)), Ex), loglik)    
        if (isoverdisp== TRUE){
            message(paste("Overdispersion estimate:", round(overdisp.est,4)))
            if(isoverdisp == TRUE & is.null(overdisp.est)) warning("No overdispersion for quasi-Poisson model. Please check.")
        } else {
            message("no overdispersion being estimated")
        }
    }
    #find K clusters
    K <- seq(from=0, to=nrow(wLambda))
    #calc
    if (isoverdisp== TRUE){
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


