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
        print("There are NA likelihoods")
    }
    return(wi)
}

####NOTE: this overdisp function is part of clusso. Omit after testing


#' @title overdisp
#' @description Calculate the overdispersion parameter (\eqn{\phi}),
#' @param offset_reg Object of class glm or lm.
#' @param overdispfloor Boolean. When \code{TRUE}, overdispersion parameter (\eqn{\phi}) is limited to not be less than 1. If \code{FALSE}, underdispersion can be estimated. 
#' @return An estimate of \eqn{\phi} (overdispersion or underdispersion)
overdisp <- function(offset_reg, overdispfloor = TRUE) {
    # if(sim==TRUE){
    #     stopifnot(inherits(offset_reg[[1]], c("glm", "lm")))
    #     phi <- max(unlist(lapply(1:length(offset_reg), function(i) deviance(offset_reg[[i]])/df.residual(offset_reg[[i]]))))
    # }
    # else{
    stopifnot(inherits(offset_reg, c("glm", "lm")))
    phi <- max(unlist(deviance(offset_reg)/df.residual(offset_reg)))
    #}
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

#'@title binomLik
#'@description Binomial-based likelihood LRT comparing cluster model to null model.
#'@param Ex Vector of expected counts.
#'@param Yx Vector of observed counts.
#'@param sparsemat Large sparsematrix where rows are potential clusters and columns are space-time locations.
#'@export
#'@return Likelihood for each potential cluster.
binomLik <-function(nx, Yx, sparsemat){
    #outExp <- sparsemat%*%Ex
    outnx <- sparsemat%*%nx
    outObs <- sparsemat%*%Yx
    #calc Lambda
    lambdahat <- outObs/outExp
    Lambda <- as.vector(lambdahat)*sparsemat #big Lambda matrix
    Lambda_dense <- as.matrix(Lambda)
    Lambda_dense[Lambda_dense == 0] <- 1
    #Get scaled likelihood
    
    #Lik <- (((outObs/nx)/(sum(outObs)/sum(nx)))^(outobs))*(((sum(outObs)-outObs)/(sum(nx)-nx))/(sum(outObs)/sum(nx)))^(sum(outObs)-outObs)
    #((outObs/outExp)/(sum(outObs)/sum(outExp)))^outObs #TODO CHECK THIS
    outnt <- sum(outnx)
    outObst <- sum(outObs)
    Lik <- (((outObs/outnx)^outObs)*(1-(outObs/outnx))^(outnx-outObs))*(((outObst - outObs)/(outnt-outnx))^(outObst-outObs))*(1-((outObst-outObs)/(outnt - outnx)))^(outnt - outnx - outObst+outObs)
    
    outlogLik <- log(Lik)
    outlogLik_scaled <- outlogLik-max(outlogLik)
    Lik <- exp(outlogLik_scaled)
    return(list(Lik=Lik,
                Lambda_dense=Lambda_dense))
}

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
#'@param bylocation If clusters should be identified by maximum location (\code{TRUE}) or maximum potential cluster (\code{FALSE}). Default is \code{TRUE}.
#'@param model A string specifying which model to use, Poisson or binomial. For Poisson, specify \code{"poisson"} and both the Poisson and quasi-Poisson model results are returned. For binomial, specify \code{"binomial"}.
#'@param cv option for cross-validation instead of AIC/BIC. Default is set to FALSE
#'@param overdisp.est Overdispersion parameter estimated across all simulations (max).
#'@return Returns list for each iteration with weighted relative risks by location inside identified cluster.
#'@export
detectclusters <- function(sparsemat, Ex, Yx,numCenters,Time, maxclust,bylocation=TRUE, model=c("poisson", "binomial"),
                           cv=FALSE,  overdisp.est) {
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
    if(bylocation==FALSE){
        print("Cluster detection by potential cluster")
        res <- bycluster(Lik, Lambda_dense, tsparsemat, maxclust)
        selection <- clusterselect(res[[1]], Yx, Ex, model,maxclust, numCenters, Time, quasi,cv=FALSE,overdisp.est)
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
        print("Cluster detection by location")
        res <- bylocation(Lik, Lambda_dense, tsparsemat, maxclust)
        selection <- clusterselect(res[[1]], Yx, Ex, model,maxclust, numCenters, Time, quasi,cv=FALSE,overdisp.est)
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
            print("no overdispersion being estimated")
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


bucklandbounds <- function(thetai,thetaa, w_q,sparsematrix, outExp,overdisp.est, transform=NULL, cellrates=FALSE) {
    if(cellrates == TRUE){
        if(!is.null(transform)){
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*(1/(thetai[k,]*outExp)))
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) 1/thetai[k,]*outExp)
            }
            withintheta <- (log(thetai) - log(thetaa))^2
            wtnbtn <- sapply(1:nrow(sparsematrix), function(k) sqrt(varthetai[,k] + withintheta[k,]))
            varthetas_w <- matrix(w_q, nrow = 1)%*%t(wtnbtn)
            var_thetaa <- (as.vector(varthetas_w))^2
            UBa = exp(as.vector(log(thetaa)) + 1.96*sqrt(var_thetaa))
            LBa = exp(as.vector(log(thetaa)) - 1.96*sqrt(var_thetaa))
            
        } else {
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k,]/outExp)
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k,]/outExp)
            }
            withintheta <- (thetai - thetaa)^2
            wtnbtn <- sapply(1:nrow(sparsematrix), function(k) sqrt(varthetai[,k] + withintheta[k,]))
            varthetas_w <- matrix(w_q, nrow = 1)%*%t(wtnbtn)
            var_thetaa <- (as.vector(varthetas_w))^2
            UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
            LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
        }
    } else {
        if(!is.null(transform)){
            #print("log-scale")
            #log transform
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*(1/(thetai[k]*outExp[k])))
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(thetai[k]*outExp[k]))
            }
            withintheta <- (log(thetai) - log(thetaa))^2
            varthetas_w <- sum(w_q*sqrt(varthetai + withintheta))^2
            var_thetaa <- (as.vector(varthetas_w))^2
            UBa = exp(as.vector(log(thetaa)) + 1.96*sqrt(var_thetaa))
            LBa = exp(as.vector(log(thetaa)) - 1.96*sqrt(var_thetaa))
            
        } else {
            #print("No transform")
            #no transform
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k]/outExp[k])
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k]/outExp[k])
            }
            withintheta <- (thetai - thetaa)^2
            varthetas_w <- sum(w_q*sqrt(varthetai + withintheta))
            var_thetaa <- (as.vector(varthetas_w))^2
            UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
            LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
        }
    }
    
    return(list(buckland.LB = LBa,
                clusterMA = thetaa,
                buckland.UB = UBa))
}


maw2 <- function(thetai,thetaa, w_q,sparsematrix, outExp, overdisp.est, transform=NULL, cellrates=FALSE) {
    if(cellrates==TRUE){
        if(!is.null(transform)){
            #print("log-scale")
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*(1/(thetai[k,]*outExp)))
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(thetai[k,]*outExp))
            }
            withintheta <- (log(thetai) - log(thetaa))^2
            wtnbtn <- sapply(1:nrow(sparsematrix), function(k) (varthetai[,k] + withintheta[k,]))
            #var_thetaa <- sum(w_q*(varthetai + withintheta))
            var_thetaa <- as.vector(matrix(w_q, nrow=1)%*%t(wtnbtn))
            UBa = exp(as.vector(log(thetaa)) + 1.96*sqrt(var_thetaa))
            LBa = exp(as.vector(log(thetaa)) - 1.96*sqrt(var_thetaa))
        } else{
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k,]/outExp)
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k,]/outExp)
            }
            withintheta <- (thetai - thetaa)^2
            wtnbtn <- sapply(1:nrow(sparsematrix), function(k) (varthetai[,k] + withintheta[k,]))
            #var_thetaa <- sum(w_q*(varthetai + withintheta))
            var_thetaa <- as.vector(matrix(w_q, nrow=1)%*%t(wtnbtn))
            UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
            LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
        }
    } else {
        if(!is.null(transform)){
            #print("log-scale")
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*(1/(thetai[k]*outExp[k])))
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(thetai[k]*outExp[k]))
            }
            withintheta <- (log(thetai) - log(thetaa))^2
            var_thetaa <- sum(w_q*(varthetai + withintheta))
            var_thetaa <- as.vector(var_thetaa)
            UBa = exp(as.vector(log(thetaa)) + 1.96*sqrt(var_thetaa))
            LBa = exp(as.vector(log(thetaa)) - 1.96*sqrt(var_thetaa))
        } else{
            if(!is.null(overdisp.est)){
                varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k]/outExp[k])
            } else {
                varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k]/outExp[k])
            }
            withintheta <- (thetai - thetaa)^2
            var_thetaa <- sum(w_q*(varthetai + withintheta))
            var_thetaa <- as.vector(var_thetaa)
            UBa = as.vector(thetaa) + 1.96*sqrt(var_thetaa)
            LBa = as.vector(thetaa) - 1.96*sqrt(var_thetaa)
        }
    }
    return(list(maw2.LB = LBa,
                clusterMA = thetaa,
                maw2.UB = UBa))
}


mata_tailareazscore <- function(thetaii, thetaaa, se.thetaii, w_q, alpha){
    thetaii <- as.vector(thetaii)
    zval <- (thetaaa - thetaii)/se.thetaii
    zpnorm <- pnorm(zval)
   # print(str(zpnorm))
#    print(str(w_q))
    w_zpnorm <- sum((w_q*zpnorm))-alpha
    
}

matabounds <- function(thetai,thetaa, w_q,sparsematrix, outExp, overdisp.est,transform=c("none","log", "sqrt"), cellrates=FALSE) {
    #print(cellrates)
    #NT <- rowSums(sparsematrix)
    if(is.null(transform)){
        none = matabounds_none(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp, cellrates)
    }
    switch(transform,
           none = matabounds_none(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp, cellrates),
           log = matabounds_log(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp,cellrates),
           sqrt = matabounds_sqrt(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp,cellrates))
}

matabounds_none <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp, cellrates) {
    if(cellrates==TRUE){
        #browser()
        if(!is.null(overdisp.est)){
            varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k,]/outExp)
        } else {
            varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k,]/outExp)
        }
        mataLB <- sapply(1:ncol(thetai), 
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= thetai[,k],
                          se.thetaii=sqrt(varthetai[k,]),
                          w_q=w_q, alpha=0.025, tol=1e-8)$root)    
        
        mataUB <- sapply(1:ncol(thetai),
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= thetai[,k],
                          se.thetaii=sqrt(varthetai[k,]),
                          w_q=w_q, alpha=1-0.025, tol=1e-8)$root) 
    } else{
        if(!is.null(overdisp.est)){
            varthetai <- sapply(1:nrow(sparsematrix), function(k) overdisp.est*thetai[k]/outExp[k])
        } else {
            varthetai <- sapply(1:nrow(sparsematrix), function(k) thetai[k]/outExp[k])
        }
        mataLB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= thetai,
                          se.thetaii=sqrt(varthetai),
                          w_q=w_q, alpha=0.025, tol=1e-8)$root    
        
        mataUB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= thetai,
                          se.thetaii=sqrt(varthetai),
                          w_q=w_q, alpha=1-0.025, tol=1e-8)$root 
    }
    
    return(list(mata.LB = mataLB,
                clusterMA = thetaa,
                mata.UB = mataUB))
}

matabounds_sqrt <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp, cellrates) {
   # print("sqrt version")
    if(cellrates==TRUE){
    #    print("cellrates true")
        #print(str(outExp))
        #browser()
        Tvarthetai <-  1/(4*outExp)
        # test <- sapply(1:5, 
        #                function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
        #                                    thetaii= sqrt(thetai[,k]),
        #                                    se.thetaii=sqrt(Tvarthetai[k]),
        #                                    w_q=w_q, alpha=0.025, tol=1e-8)$root)
        
        #TODO: fix herel it should be thetai[,k] where k corresponds to the kth cell of cellsix, not 1:5
        mataLB <- sapply(1:ncol(thetai), 
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= sqrt(thetai[,k]),
                          se.thetaii=sqrt(Tvarthetai[k]),
                          w_q=w_q, alpha=0.025, tol=1e-8)$root)
        
        mataUB <- sapply(1:ncol(thetai),
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= sqrt(thetai[,k]),
                          se.thetaii=sqrt(Tvarthetai[k]),
                          w_q=w_q, alpha=1-0.025, tol=1e-8)$root)
    } else{
        Tvarthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(4*outExp[k]))
        mataLB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= sqrt(thetai),
                          se.thetaii=sqrt(Tvarthetai),
                          w_q=w_q, alpha=0.025, tol=1e-8)$root
        
        mataUB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= sqrt(thetai),
                          se.thetaii=sqrt(Tvarthetai),
                          w_q=w_q, alpha=1-0.025, tol=1e-8)$root
    }
    return(list(matasqrt.LB = (mataLB)^2,
                clusterMA = thetaa,
                matasqrt.UB = (mataUB)^2))
}

matabounds_log <- function(thetai,thetaa, w_q,sparsematrix, overdisp.est, outExp, cellrates) {
    if(cellrates==TRUE){
        logTvarthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(thetai[k,]*outExp))
        mataLB <- sapply(1:ncol(thetai),
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= log(thetai[,k]),
                          se.thetaii=sqrt(logTvarthetai[k,]),
                          w_q=w_q, alpha=0.025, tol=1e-8)$root)
        
        mataUB <- sapply(1:ncol(thetai),
                         function(k) uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= log(thetai[,k]),
                          se.thetaii=sqrt(logTvarthetai[k,]),
                          w_q=w_q, alpha=1-0.025, tol=1e-8)$root)
    } else {
        logTvarthetai <- sapply(1:nrow(sparsematrix), function(k) 1/(thetai[k]*outExp[k]))   
        mataLB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= log(thetai),
                          se.thetaii=sqrt(logTvarthetai),
                          w_q=w_q, alpha=0.025, tol=1e-8)$root
        
        mataUB <- uniroot(f=mata_tailareazscore, interval=c(-10, 10),
                          thetaii= log(thetai),
                          se.thetaii=sqrt(logTvarthetai),
                          w_q=w_q, alpha=1-0.025, tol=1e-8)$root
    }
    return(list(matalog.LB = exp(mataLB),
                clusterMA = thetaa,
                matalog.UB = exp(mataUB)))
}

selectuniqRR <- function(uniqRRs){
    clusterRR_i <- rep(NA, nrow(uniqRRs))
    flag1 <- sapply(1:nrow(uniqRRs), function(k) length(unique(uniqRRs[k,])))
    non1s <- uniqRRs[which(flag1==2),]
    clusterRR_i[which(flag1==2)] <- sapply(1:nrow(non1s), function(k) non1s[k,which(non1s[k,]!=1)])
    clusterRR_i[which(flag1==1)] <- 1    
    return(unlist(clusterRR_i))
    
}

##############################################################################
#SIMULATION HELPER FUNCTIONS
##############################################################################

#nonma
nonma <- function(cluster_thetaa,res, clusterRR_ilarge,wslarge,idix, IC, transform, bylocation=TRUE){
    if(bylocation==FALSE){
        print("by pc")
        # browser()
        if(transform=="log"){
            if(IC=="aic") {
                #clusterRRlarge <- lapply(1:length(idix), function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]])[2])
                clusterRRlarge <- lapply(1:length(idix), function(j) unique(res[[idix[[j]]]]$Lambda_dense[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic],])[2])
                se_clusterRRlarge <- lapply(1:length(idix), function(j)sqrt(1/(clusterRRlarge[[j]]*outExp[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]])))
                
            } else {
                # clusterRRlarge <- lapply(1:length(idix), function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]])[2])
                # se_clusterRRlarge <- lapply(1:length(idix), function(j)sqrt(1/(clusterRRlarge[[j]]*outExp[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]])))
                clusterRRlarge <- lapply(1:length(idix), function(j) unique(res[[idix[[j]]]]$Lambda_dense[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic],])[2])
                se_clusterRRlarge <- lapply(1:length(idix), function(j)sqrt(clusterRRlarge[[j]]/outExp[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]]))
                
            }
            nonma.theta.time <- system.time(nonma.theta <- lapply(1:length(idix), function(i) cbind(lb=exp(log(cluster_thetaa[[i]])-1.96*se_clusterRRlarge[[i]]), 
                                                                                                    clusterMA = cluster_thetaa[[i]],
                                                                                                    ub=exp(log(cluster_thetaa[[i]])+1.96*se_clusterRRlarge[[i]]))))
        } else {
            if(IC=="aic") {
                #clusterRRlarge <- lapply(1:length(idix), function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]])[2])
                clusterRRlarge <- lapply(1:length(idix), function(j) unique(res[[idix[[j]]]]$Lambda_dense[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic],])[2])
                se_clusterRRlarge <- lapply(1:length(idix), function(j)sqrt(clusterRRlarge[[j]]/outExp[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]]))
                
            } else {
                #browser()
                #clusterRRlarge <- lapply(1:length(idix), function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]])[2])
                clusterRRlarge <- lapply(1:length(idix), function(j) unique(res[[idix[[j]]]]$Lambda_dense[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic],])[2])
                se_clusterRRlarge <- lapply(1:length(idix), function(j)sqrt(clusterRRlarge[[j]]/outExp[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]]))
                
                
            }
            nonma.theta.time <- system.time(nonma.theta <- lapply(1:length(idix), function(i) cbind(lb=cluster_thetaa[[i]]-1.96*se_clusterRRlarge[[i]], 
                                                                                                    clusterMA = cluster_thetaa[[i]],
                                                                                                    ub=cluster_thetaa[[i]]+1.96*se_clusterRRlarge[[i]])))
            
        }
    } else {
        print("by loc")
        # browser()
        if(transform=="log"){
            if(IC=="aic") {
                clusterRRlarge <- lapply(1:length(idix), function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]])[2])
                se_clusterRRlarge <- lapply(1:length(idix), function(j)sqrt(1/(clusterRRlarge[[j]]*outExp[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]])))
                
            } else {
                clusterRRlarge <- lapply(1:length(idix), function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]])[2])
                se_clusterRRlarge <- lapply(1:length(idix), function(j)sqrt(1/(clusterRRlarge[[j]]*outExp[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]])))
                
            }
            nonma.theta.time <- system.time(nonma.theta <- lapply(1:length(idix), function(i) cbind(lb=exp(log(cluster_thetaa[[i]])-1.96*se_clusterRRlarge[[i]]), 
                                                                                                    clusterMA = cluster_thetaa[[i]],
                                                                                                    ub=exp(log(cluster_thetaa[[i]])+1.96*se_clusterRRlarge[[i]]))))
        } else {
            if(IC=="aic") {
                clusterRRlarge <- lapply(1:length(idix), function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]])[2])
                se_clusterRRlarge <- lapply(1:length(idix), function(j)sqrt(clusterRRlarge[[j]]/outExp[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]]))
                
            } else {
                #browser()
                clusterRRlarge <- lapply(1:length(idix), function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]])[2])
                se_clusterRRlarge <- lapply(1:length(idix), function(j)sqrt(clusterRRlarge[[j]]/outExp[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]]))
                
            }
            nonma.theta.time <- system.time(nonma.theta <- lapply(1:length(idix), function(i) cbind(lb=cluster_thetaa[[i]]-1.96*se_clusterRRlarge[[i]], 
                                                                                                    clusterMA = cluster_thetaa[[i]],
                                                                                                    ub=cluster_thetaa[[i]]+1.96*se_clusterRRlarge[[i]])))
            
        }
    }
   
    return(list(nonma.theta.time = nonma.theta.time[[3]],
                nonma.theta = nonma.theta))
}


#nonma_asymp
nonma_asymp <- function(cluster_thetaa,res,clusterRR_ilarge,wslarge,idix, IC, transform, bylocation=TRUE){
    if(bylocation==FALSE){
        if(transform=="log"){
            if(IC=="aic") {
                # clusterRRlarge <- lapply(1:length(idix), 
                #                          function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]])[2])
                clusterRRlarge <- lapply(1:length(idix), 
                                         function(j) unique(res[[idix[[j]]]]$Lambda_dense[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic],])[2])
                se_clusterRRlarge_asymp <- lapply(1:length(idix), function(j) sqrt(1/(clusterRRlarge[[j]]*outObs[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]])))
                
            } else {
                # clusterRRlarge <- lapply(1:length(idix), 
                #                          function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]])[2])
                clusterRRlarge <- lapply(1:length(idix), 
                                         function(j) unique(res[[idix[[j]]]]$Lambda_dense[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic],])[2])
                se_clusterRRlarge_asymp <- lapply(1:length(idix), function(j) sqrt(1/(clusterRRlarge[[j]]*outObs[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]])))
                
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- lapply(1:length(idix), function(i) cbind(lbasymp=exp(log(cluster_thetaa[[i]])-1.96*se_clusterRRlarge_asymp[[i]]), 
                                                                                                                clusterMA = cluster_thetaa[[i]],
                                                                                                                ubasymp=exp(log(cluster_thetaa[[i]])+1.96*se_clusterRRlarge_asymp[[i]]))))
        } else{
            if(IC=="aic") {
                # clusterRRlarge <- lapply(1:length(idix), 
                #                          function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]])[2])
                clusterRRlarge <- lapply(1:length(idix), 
                                         function(j) unique(res[[idix[[j]]]]$Lambda_dense[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic],])[2])
                se_clusterRRlarge_asymp <- lapply(1:length(idix), function(j) sqrt(clusterRRlarge[[j]]/outObs[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]]))
                
            } else {
                # clusterRRlarge <- lapply(1:length(idix), 
                #                          function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]])[2])
                clusterRRlarge <- lapply(1:length(idix), 
                                         function(j) unique(res[[idix[[j]]]]$Lambda_dense[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic],])[2])
                se_clusterRRlarge_asymp <- lapply(1:length(idix), function(j) sqrt(clusterRRlarge[[j]]/outObs[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]]))
                
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- lapply(1:length(idix), function(i) cbind(lbasymp=cluster_thetaa[[i]]-1.96*se_clusterRRlarge_asymp[[i]], 
                                                                                                                clusterMA = cluster_thetaa[[i]],
                                                                                                                ubasymp=cluster_thetaa[[i]]+1.96*se_clusterRRlarge_asymp[[i]])))
        }
        
    } else {
        if(transform=="log"){
            if(IC=="aic") {
                clusterRRlarge <- lapply(1:length(idix), 
                                         function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]])[2])
                se_clusterRRlarge_asymp <- lapply(1:length(idix), function(j) sqrt(1/(clusterRRlarge[[j]]*outObs[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]])))
                
            } else {
                clusterRRlarge <- lapply(1:length(idix), 
                                         function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]])[2])
                se_clusterRRlarge_asymp <- lapply(1:length(idix), function(j) sqrt(1/(clusterRRlarge[[j]]*outObs[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]])))
                
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- lapply(1:length(idix), function(i) cbind(lbasymp=exp(log(cluster_thetaa[[i]])-1.96*se_clusterRRlarge_asymp[[i]]), 
                                                                                                                clusterMA = cluster_thetaa[[i]],
                                                                                                                ubasymp=exp(log(cluster_thetaa[[i]])+1.96*se_clusterRRlarge_asymp[[i]]))))
        } else{
            if(IC=="aic") {
                clusterRRlarge <- lapply(1:length(idix), 
                                         function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]])[2])
                se_clusterRRlarge_asymp <- lapply(1:length(idix), function(j) sqrt(clusterRRlarge[[j]]/outObs[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.aic]]))
                
            } else {
                clusterRRlarge <- lapply(1:length(idix), 
                                         function(j) unique(res[[idix[[j]]]]$Lambda_dense[,res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]])[2])
                se_clusterRRlarge_asymp <- lapply(1:length(idix), function(j) sqrt(clusterRRlarge[[j]]/outObs[[j]]@x[res[[idix[[j]]]]$maxid[res[[idix[[j]]]]$selection.bic]]))
                
            }
            nonma_asymp.theta.time <- system.time(nonma_asymp.theta <- lapply(1:length(idix), function(i) cbind(lbasymp=cluster_thetaa[[i]]-1.96*se_clusterRRlarge_asymp[[i]], 
                                                                                                                clusterMA = cluster_thetaa[[i]],
                                                                                                                ubasymp=cluster_thetaa[[i]]+1.96*se_clusterRRlarge_asymp[[i]])))
        }
    }
    
    return(list(nonma_asymp.theta.time = nonma_asymp.theta.time[[3]],
                nonma_asymp.theta = nonma_asymp.theta))
}


calcbounds <- function(id, IC, res, bylocation){
    #do all diagnostics
    idix <- which(id!=0)
    #print(idix)
    print(paste0("calculating bounds for ", IC))
    #prep
    if(IC=="aic"){
        wslarge <- lapply(1:length(idix), function(j) res[[idix[[j]]]]$wtMAT[,res[[idix[[j]]]]$selection.aic])
    } else {
        wslarge <- lapply(1:length(idix), function(j) res[[idix[[j]]]]$wtMAT[,res[[idix[[j]]]]$selection.bic])
    }
    clusterRR_uniqlarge <- lapply(1:length(idix), function(j) sapply(1:nrow(res[[idix[[j]]]]$Lambda_dense), 
                                                                     function(k) unique(res[[idix[[j]]]]$Lambda_dense[k,]))) 
    
    clusterRR_ilarge <- lapply(1:length(idix), function(i) rep(NA, 66870))
    clusterRR_uniq_ilarge <- lapply(1:length(idix), function(j) as.matrix(do.call(rbind, clusterRR_uniqlarge[[j]]), ncol=2))
    clusterRR_ilarge <- lapply(1:length(idix), function(j) selectuniqRR(clusterRR_uniq_ilarge[[j]]))
    cluster_thetaa <- lapply(1:length(idix), function(j) sum(clusterRR_ilarge[[j]]*wslarge[[j]]))
    
    
    
    #Perform
    outnonma.time <- system.time(outnonma <- nonma(cluster_thetaa, res, clusterRR_ilarge, wslarge, idix, IC=IC, transform="none", bylocation))
    outnonmaTlog.time <- system.time(outnonmaTlog <- nonma(cluster_thetaa, res, clusterRR_ilarge, wslarge, idix, IC=IC, transform="log", bylocation))
    outnonma_asymp.time <- system.time(outnonma_asymp <- nonma_asymp(cluster_thetaa ,res, clusterRR_ilarge, wslarge, idix, IC=IC, transform="none", bylocation))
    outnonma_asympTlog.time <- system.time(outnonma_asympTlog <- nonma_asymp(cluster_thetaa ,res, clusterRR_ilarge, wslarge, idix, IC=IC, transform="log", bylocation))
    print("nonma finished")
    outbuck.theta.time <- system.time(outbuck.theta <- lapply(1:length(idix), function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]], 
                                                                                                 thetaa = cluster_thetaa[[i]], 
                                                                                                 w_q=wslarge[[i]], 
                                                                                                 sparsematrix=t(sparsematrix), 
                                                                                                 outExp[[i]],overdisp.est = NULL)))
    outbuckTlog.theta.time <- system.time(outbuckTlog.theta <- lapply(1:length(idix), function(i) bucklandbounds(thetai=clusterRR_ilarge[[i]],
                                                                                                         thetaa =cluster_thetaa[[i]], 
                                                                                                         w_q=wslarge[[i]], 
                                                                                                         sparsematrix=t(sparsematrix),
                                                                                                         outExp[[i]],
                                                                                                         overdisp.est = NULL, 
                                                                                                         transform=TRUE)))
    print("buckland finished")
    outmaw2.theta.time <- system.time(outmaw2.theta <- lapply(1:length(idix), function(i) maw2(thetai=clusterRR_ilarge[[i]], 
                                                                                       thetaa = cluster_thetaa[[i]], 
                                                                                       w_q=wslarge[[i]],
                                                                                       sparsematrix=t(sparsematrix ), 
                                                                                       outExp[[i]],overdisp.est = NULL)))
    
    outmaw2Tlog.theta.time <- system.time(outmaw2Tlog.theta  <- lapply(1:length(idix), function(i) maw2(thetai=clusterRR_ilarge[[i]], 
                                                                                                thetaa = cluster_thetaa[[i]], 
                                                                                                w_q=wslarge[[i]], 
                                                                                                sparsematrix=t(sparsematrix), 
                                                                                                outExp[[i]], 
                                                                                                overdisp.est = NULL,
                                                                                                transform=TRUE)))
    print("maw2 finished")
    outmata.theta.time <- system.time(outmata.theta <- lapply(1:length(idix), function(i) matabounds(thetai=clusterRR_ilarge[[i]], 
                                                                                             thetaa = cluster_thetaa[[i]], 
                                                                                             w_q=wslarge[[i]], 
                                                                                             sparsematrix=t(sparsematrix ), 
                                                                                             outExp = outExp[[i]],
                                                                                             overdisp.est = NULL,
                                                                                             transform="none")))
    outmataTlog.theta.time <- system.time(outmataTlog.theta <- lapply(1:length(idix), function(i) matabounds(thetai=clusterRR_ilarge[[i]], 
                                                                                                     thetaa = cluster_thetaa[[i]], 
                                                                                                     w_q=wslarge[[i]], 
                                                                                                     sparsematrix=t(sparsematrix ),
                                                                                                     outExp = outExp[[i]], 
                                                                                                     overdisp.est = NULL,
                                                                                                     transform="log")))
    print("mata finished")
    #return
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
        outmataTlog.theta.time = outmataTlog.theta.time[[3]] ))
    
}
