set.seed(20190326)
library(clusso)
source("R/clustack.R")


#0) Setup
load("data/japanbreastcancer.RData")
cases <- japanbreastcancer$death
expected <- japanbreastcancer$expdeath
# centroids <- japanbreastcancer
# periods <- japanbreastcancer
x <- utmJapan$utmx/1000
y <- utmJapan$utmy/1000
japan.poly2 <- dframe.poly2[,2:3]
japan.prefect2 <- dframe.prefect2[,2:5]


#set global
rMax <- 20 
Time <- 5
#create set of potential clusters based on distances
potentialclusters <- clusso::clusters2df(x,y,rMax, utm = TRUE, length(x))
n_uniq <- length(unique(potentialclusters$center))
numCenters <- n_uniq
#create giant sparse design matrix (single potential clusters)
sparsematrix <- spacetimeMat(potentialclusters, numCenters, Time) 


#################################
#TESTER SIM
#################################
#quick function to recode
reval <- function(probs, ix){
    probs[ix] <-1
    return(probs)
}
#nsims
nsim <- 10
#theta - overdispersion
theta <- 500
maxclust <- 11


#scenarios
centers <- c(150,35)
radii <- c(9,18)
timeperiod <- c(3:5)
#timeperiods <- list(c(3:5), c(1:2))
risk.ratios <- c(1, 1.1, 1.5, 2)
table.detection.loc <- NULL
table.detection.pc <- NULL
table.detection.clusso <- NULL
eps <- 3

#test
cent <- 150
rad <- 11
risk <- 2
tim <- c(3:5)

# Start the clock!
ptm <- proc.time()
for (cent in centers){
    for (rad in radii){
        for (risk in risk.ratios){
            print(paste0("Params:", "center: ",cent,"; radius: ",
                         rad, "; Timeperiods: ", as.numeric(paste(tim, collapse = "")), 
                         "; RR: ", risk, "; Theta :", theta))
            simid <- 1:nsim
            #put the cluster in
            clusters <- clusso::clusters2df(x,y,rMax, utm = TRUE, length(x))
            n <- length(x)
            init <- clusso::setVectors(japanbreastcancer$period, 
                                       japanbreastcancer$expdeath, japanbreastcancer$death,
                                       covars=NULL, Time=Time)
            E1 <- init$E0
            #Ex <- clusso::scale(init, Time)
            Yx <- init$Y.vec
            vectors <- list(Period = init$Year, Ex=init$E0, E0_0=init$E0, Y.vec=init$Y.vec, covars = NULL)    
            n_uniq <- length(unique(clusters$center))
            numCenters <- n_uniq
            tmp <- clusters[clusters$center==cent,]
            cluster <- tmp[(tmp$r <= rad),]
            rr = matrix(1, nrow=n, ncol=Time)
            rr[cluster$last, tim[1]:tail(tim, n=1)] <- risk
            E1 <- as.vector(rr)*init$E0
            message(paste("Running model for periods",tim[1],"through", tail(tim, n=1)))
            #simulate data here
            YSIM <- lapply(1:nsim, function(i) MASS::rnegbin(E1, theta = theta))
            Ex <- clusso::scale_sim(YSIM, init, nsim, Time)
            #####################################################################################
            #####################################################################################
            #####################################################################################
            #CLUSTER DETECTION BY LOCATION
            #####################################################################################
            #####################################################################################
            #####################################################################################
            #run superclust by location for each sim (apply)
            sim_superclust_loc <- lapply(1:nsim, function(i) detectclusters(sparsematrix, Ex[[i]], YSIM[[i]],
                                                                            numCenters, Time, maxclust,
                                                                            bylocation = TRUE, model="poisson"))
            print(sim.i <- paste0("sim","_", "center",cent,"_" ,"radius", rad, "_",
                                  "risk", risk, "_", "theta", as.character(theta),
                                  as.numeric(paste(tim, collapse = "")), "_moderateoverdisperson", "_superclustLOC"))
            filename <- paste0(sim.i,".RData")
            #save .RData
            save(sim_superclust_loc, file=filename)
            ##################################
            #DIAGNOSTICS: #calc power and FB rate
            ##################################
            #Which PCs overlap true cluster?
            rrbin_cluster <- matrix(as.vector(ifelse(rr!=1,1,0)),nrow=1)
            clusteroverlap<- rrbin_cluster%*%sparsematrix
            rrbin_inside <- ifelse(sparsematrix%*%t(clusteroverlap)!=0,1,0)
            #what was identified in each sim by IC
            ident.bic <- lapply(1:nsim, function(i) 
                round(sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.bic,],eps))
            ident.aic <- lapply(1:nsim, function(i)
                round(sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.aic,],eps))
            ident.aicc <- lapply(1:nsim, function(i) 
                round(sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.aicc,],eps))
            #1) Did it find anything INSIDE the cluster?
            incluster.bic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.bic[[i]]==1,0,1),ncol=1))
            incluster.aic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.aic[[i]]==1,0,1),ncol=1))
            incluster.aicc <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.aicc[[i]]==1,0,1),ncol=1))
            #calc power
            pow.bic <- sum(ifelse(unlist(incluster.bic)!=0,1,0))/nsim
            pow.aic <- sum(ifelse(unlist(incluster.aic)!=0,1,0))/nsim
            pow.aicc <-sum(ifelse(unlist(incluster.aicc)!=0,1,0))/nsim
            outpow.bic <- paste0(pow.bic*100, "%")
            outpow.aic <- paste0(pow.aic*100, "%")
            outpow.aicc <- paste0(pow.aicc*100, "%")
            #2) Did it find anything OUTSIDE the cluster?
            rrbin_outside <- ifelse(sparsematrix%*%t(clusteroverlap)==0,1,0) 
            #this should be everything that doesn't touch the cluster
            outcluster.bic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.bic[[i]]==1,0,1),ncol=1))
            outcluster.aic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aic[[i]]==1,0,1),ncol=1))
            outcluster.aicc <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aicc[[i]]==1,0,1),ncol=1))
            #calc FP rate
            fp.bic <- sum(ifelse(unlist(outcluster.bic)!=0,1,0))/nsim
            fp.aic <- sum(ifelse(unlist(outcluster.aic)!=0,1,0))/nsim
            fp.aicc <- sum(ifelse(unlist(outcluster.aicc)!=0,1,0))/nsim
            outfp.bic <- paste0(fp.bic*100, "%")
            outfp.aic <- paste0(fp.aic*100, "%")
            outfp.aicc <- paste0(fp.aicc*100, "%")
            ##################################
            #plot probability maps
            ##################################
            #create empties
            vec <- rep(0, 208 * Time)
            position.bic <- list(vec)[rep(1, nsim)]
            position.aic <- list(vec)[rep(1, nsim)]
            position.aicc <- list(vec)[rep(1, nsim)]
            #recode identified cells as 1's, all other zeros
            ix.bic <- lapply(1:nsim, function(i) which(ifelse(ident.bic[[i]]==1,0,1) ==1))
            ix.aic <- lapply(1:nsim, function(i) which(ifelse(ident.aic[[i]]==1,0,1) ==1))
            ix.aicc <- lapply(1:nsim, function(i) which(ifelse(ident.aicc[[i]]==1,0,1) ==1))
            #creatematrix by location (rows) and sim (cols) with 1's indicating selection by superlearner
            simindicator.bic <- mapply(reval, position.bic, ix.bic)
            simindicator.aic <- mapply(reval, position.aic, ix.aic)
            simindicator.aicc <- mapply(reval, position.aicc, ix.aicc)
            #find probability of detection for each location in time
            probs.bic <- Matrix::rowSums(simindicator.bic)/nsim
            probs.aic <- Matrix::rowSums(simindicator.aic)/nsim
            probs.aicc <- Matrix::rowSums(simindicator.aicc)/nsim
            #map probability detections by IC to grey scale
            colprob <- colormapping(list(probs.bic,
                                         probs.aic,
                                         probs.aicc,
                                         as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
            #plot map with probability detection by each IC
            probplotmapAllIC(colprob,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_loc.pdf"))
            ##################################
            #RR maps
            ##################################
            #plot each RR for each sim            
            lapply(1:nsim, function(i) 
                plotmapAllIC(res.bic = sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.bic,],
                             res.aic = sim_superclust_loc[[i]]$wLambda[sim_superclust_loc[[i]]$selection.aic,],
                                                    oracle = rr,
                                                    genpdf = TRUE,
                                                    pdfname = paste0(sim.i,"_loc_simID",simid[i],".pdf")))
            ##################################
            #Add sim results to table
            ##################################
            tabn.loc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                                    time=as.numeric(paste(tim, collapse="")),
                                    mod="ST", pow=outpow.bic, fp = outfp.bic),
                              cbind(IC="AIC",rad, risk, cent, theta,
                                    time=as.numeric(paste(tim, collapse="")),
                                    mod="ST", pow=outpow.aic, fp = outfp.aic),
                              cbind(IC="AICc",rad, risk, cent, theta,
                                    time=as.numeric(paste(tim, collapse="")),
                                    mod="ST", pow=outpow.aicc, fp = outfp.aicc))
            table.detection.loc <- rbind(table.detection.loc, tabn.loc)
            
            
            #####################################################################################
            #####################################################################################
            #####################################################################################
            #CLUSTER DETECTION BY Potential Cluster
            #####################################################################################
            #####################################################################################
            #####################################################################################
            #run superclust by PC for each sim (apply)
            sim_superclust_pc<- lapply(1:nsim, function(i) detectclusters(sparsematrix, Ex[[i]], YSIM[[i]],
                                                                          numCenters, Time, maxclust,
                                                                          bylocation = FALSE, model="poisson"))
            print(sim.i <- paste0("sim","_", "center",cent,"_" ,"radius", rad, "_",
                                  "risk", risk, "_", "theta", as.character(theta),
                                  as.numeric(paste(tim, collapse = "")), "_moderateoverdisperson", "_pc"))
            filename <- paste0(sim.i,".RData")
            #save .RData
            save(sim_superclust_pc, file=filename)
            #################################
            DIAGNOSTICS: #calc power and FB rate
                #################################
            # #Which PCs overlap true cluster?
            #what was identified in each sim by IC
            ident.bic <- lapply(1:nsim, function(i) 
                round(sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.bic,],eps))
            ident.aic <- lapply(1:nsim, function(i) 
                round(sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.aic,],eps))
            ident.aicc <- lapply(1:nsim, function(i)
                round(sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.aicc,],eps))
            #1) Did it find anything INSIDE the cluster?
            incluster.bic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.bic[[i]]==1,0,1),ncol=1))
            incluster.aic <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.aic[[i]]==1,0,1),ncol=1))
            incluster.aicc <- lapply(1:nsim, function(i) rrbin_inside%*%matrix(ifelse(ident.aicc[[i]]==1,0,1),ncol=1))
            #calc power
            pow.bic <- sum(ifelse(unlist(incluster.bic)!=0,1,0))/nsim
            pow.aic <- sum(ifelse(unlist(incluster.aic)!=0,1,0))/nsim
            pow.aicc <-sum(ifelse(unlist(incluster.aicc)!=0,1,0))/nsim
            outpow.bic <- paste0(pow.bic*100, "%")
            outpow.aic <- paste0(pow.aic*100, "%")
            outpow.aicc <- paste0(pow.aicc*100, "%")
            #2) Did it find anything OUTSIDE the cluster?
            #rrbin_outside <- ifelse(sparsematrix%*%t(clusteroverlap)==0,1,0)
            #this should be everything that doesn't touch the cluster
            outcluster.bic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.bic[[i]]==1,0,1),ncol=1))
            outcluster.aic <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aic[[i]]==1,0,1),ncol=1))
            outcluster.aicc <- lapply(1:nsim, function(i) rrbin_outside%*%matrix(ifelse(ident.aicc[[i]]==1,0,1),ncol=1))
            #calc FP rate
            fp.bic <- sum(ifelse(unlist(outcluster.bic)!=0,1,0))/nsim
            fp.aic <- sum(ifelse(unlist(outcluster.aic)!=0,1,0))/nsim
            fp.aicc <- sum(ifelse(unlist(outcluster.aicc)!=0,1,0))/nsim
            outfp.bic <- paste0(fp.bic*100, "%")
            outfp.aic <- paste0(fp.aic*100, "%")
            outfp.aicc <- paste0(fp.aicc*100, "%")
            ##################################
            #plot probability maps
            ##################################
            #create empties
            vec <- rep(0, 208 * Time)
            position.bic <- list(vec)[rep(1, nsim)]
            position.aic <- list(vec)[rep(1, nsim)]
            position.aicc <- list(vec)[rep(1, nsim)]
            #recode identified cells as 1's, all other zeros
            ix.bic <- lapply(1:nsim, function(i) which(ifelse(ident.bic[[i]]==1,0,1) ==1))
            ix.aic <- lapply(1:nsim, function(i) which(ifelse(ident.aic[[i]]==1,0,1) ==1))
            ix.aicc <- lapply(1:nsim, function(i) which(ifelse(ident.aicc[[i]]==1,0,1) ==1))
            #creatematrix by location (rows) and sim (cols) with 1's indicating selection by superlearner
            simindicator.bic <- mapply(reval, position.bic, ix.bic)
            simindicator.aic <- mapply(reval, position.aic, ix.aic)
            simindicator.aicc <- mapply(reval, position.aicc, ix.aicc)
            #find probability of detection for each location in time
            probs.bic <- Matrix::rowSums(simindicator.bic)/nsim
            probs.aic <- Matrix::rowSums(simindicator.aic)/nsim
            probs.aicc <- Matrix::rowSums(simindicator.aicc)/nsim
            #map probability detections by IC to grey scale
            colprob <- colormapping(list(probs.bic,
                                         probs.aic,
                                         probs.aicc,
                                         as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
            #plot map with probability detection by each IC
            probplotmapAllIC(colprob,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_pc.pdf"))
            ##################################
            #RR maps
            ##################################
            #plot each RR for each sim
            lapply(1:nsim, function(i) plotmapAllIC(res.bic = sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.bic,],
                                                    res.aic = sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.aic,],
                                                    oracle = rr,
                                                    genpdf = TRUE,
                                                    pdfname = paste0(sim.i,"_pc_simID",simid[i],".pdf")))
            ##################################
            #Add sim results to table
            ##################################
            tabn.pc <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                                   time=as.numeric(paste(tim, collapse="")),
                                   mod="ST", pow=outpow.bic, fp = outfp.bic),
                             cbind(IC="AIC",rad, risk, cent, theta,
                                   time=as.numeric(paste(tim, collapse="")),
                                   mod="ST", pow=outpow.aic, fp = outfp.aic),
                             cbind(IC="AICc",rad, risk, cent, theta,
                                   time=as.numeric(paste(tim, collapse="")),
                                   mod="ST", pow=outpow.aicc, fp = outfp.aicc))
            table.detection.pc <- rbind(table.detection.pc, tabn.pc)
            
            #####################################################################################
            #####################################################################################
            #####################################################################################
            #CLUSTER DETECTION BY CLUSSO
            #####################################################################################
            #####################################################################################
            #####################################################################################
            #create list of dataframes
            jbcSIM <- lapply(1:nsim, function(i) cbind.data.frame(expected = Ex[[i]],
                                                                  observed = YSIM[[i]],
                                                                  period = japanbreastcancer$period))
            #run clusso
            sim_clusso <- lapply(1:nsim, function(i) clusso(df=jbcSIM[[i]],
                                 expected = expected,
                                 observed = observed,
                                 timeperiod = period,
                                 covars=FALSE,
                                 id = NULL,
                                 x = x,
                                 y = y,
                                 rMax = rMax,
                                 utm=TRUE, 
                                 analysis = "both",
                                 model = "poisson",
                                 maxclust = maxclust))
            print(sim.i <- paste0("sim","_", "center",cent,"_" ,"radius", rad, "_",
                                  "risk", risk, "_", "theta", as.character(theta),
                                  as.numeric(paste(tim, collapse = "")), "_moderateoverdisperson", "_clusso"))
            filename <- paste0(sim.i,".RData")
            #save .RData
            save(sim_superclust_pc, file=filename)
            ##################################
            #DIAGNOSTICS: #calc power and FB rate
            ##################################
            #detect 
            ##quasi-Poisson
            ident.bic.qp <- lapply(1:nsim, function(i) {
                out <- sparsematrix %*% ifelse(abs(sim_clusso[[i]]$lassoresult.qp.st$coefs.lasso.all[1:66870, sim_clusso[[i]]$lassoresult.qp.st$selections$select.qbic]) >= 1e-05,1,0);
                ifelse(out!=0,1,0)})
            ident.aic.qp <- lapply(1:nsim, function(i) {
                out <- sparsematrix %*% ifelse(abs(sim_clusso[[i]]$lassoresult.qp.st$coefs.lasso.all[1:66870, sim_clusso[[i]]$lassoresult.qp.st$selections$select.qaic]) >= 1e-05,1,0);
                ifelse(out!=0,1,0)})
            ident.aicc.qp <- lapply(1:nsim, function(i) {
                out <- sparsematrix %*% ifelse(abs(sim_clusso[[i]]$lassoresult.qp.st$coefs.lasso.all[1:66870, sim_clusso[[i]]$lassoresult.qp.st$selections$select.qaicc]) >= 1e-05,1,0);
                ifelse(out!=0,1,0)})
            ##Poisson
            ident.bic.p <- lapply(1:nsim, function(i) {
                out <- sparsematrix %*% ifelse(abs(sim_clusso[[i]]$lassoresult.p.st$coefs.lasso.all[1:66870, sim_clusso[[i]]$lassoresult.p.st$selections$select.qbic]) >= 1e-05,1,0);
                ifelse(out!=0,1,0)})
            ident.aic.p <- lapply(1:nsim, function(i) {
                out <- sparsematrix %*% ifelse(abs(sim_clusso[[i]]$lassoresult.p.st$coefs.lasso.all[1:66870, sim_clusso[[i]]$lassoresult.p.st$selections$select.qaic]) >= 1e-05,1,0);
                ifelse(out!=0,1,0)})
            ident.aicc.p <- lapply(1:nsim, function(i) {
                out <- sparsematrix %*% ifelse(abs(sim_clusso[[i]]$lassoresult.p.st$coefs.lasso.all[1:66870, sim_clusso[[i]]$lassoresult.p.st$selections$select.qaicc]) >= 1e-05,1,0);
                ifelse(out!=0,1,0)})
            #calculate
            ##quasi-Poisson
            incluster.bic.qp <- lapply(1:nsim, function(i)
                rrbin_inside%*% matrix(ident.bic.qp[[i]], ncol=1))
            incluster.aic.qp <- lapply(1:nsim, function(i)
                rrbin_inside%*% matrix(ident.aic.qp[[i]], ncol=1))
            incluster.aicc.qp <- lapply(1:nsim, function(i)
                rrbin_inside%*% matrix(ident.aicc.qp[[i]], ncol=1))
            ##Poisson
            incluster.bic.p <- lapply(1:nsim, function(i)
                rrbin_inside%*% matrix(ident.bic.p[[i]], ncol=1))
            incluster.aic.p <- lapply(1:nsim, function(i)
                rrbin_inside%*% matrix(ident.aic.p[[i]], ncol=1))
            incluster.aicc.p <- lapply(1:nsim, function(i)
                rrbin_inside%*% matrix(ident.aicc.p[[i]], ncol=1))
            # #Which PCs overlap true cluster?
            # #1) Did it find anything INSIDE the cluster?
            #calc power
            ##quasi-Poisson
            pow.bic.qp <- sum(ifelse(unlist(incluster.bic.qp)!=0,1,0))/nsim
            pow.aic.qp <- sum(ifelse(unlist(incluster.aic.qp)!=0,1,0))/nsim
            pow.aicc.qp <-sum(ifelse(unlist(incluster.aicc.qp)!=0,1,0))/nsim
            outpow.bic.qp <- paste0(pow.bic.qp*100, "%")
            outpow.aic.qp <- paste0(pow.aic.qp*100, "%")
            outpow.aicc.qp <- paste0(pow.aicc.qp*100, "%")
            ##Poisson
            pow.bic.p <- sum(ifelse(unlist(incluster.bic.p)!=0,1,0))/nsim
            pow.aic.p <- sum(ifelse(unlist(incluster.aic.p)!=0,1,0))/nsim
            pow.aicc.p <-sum(ifelse(unlist(incluster.aicc.p)!=0,1,0))/nsim
            outpow.bic.p <- paste0(pow.bic.p*100, "%")
            outpow.aic.p <- paste0(pow.aic.p*100, "%")
            outpow.aicc.p <- paste0(pow.aicc.p*100, "%")
            # #2) Did it find anything OUTSIDE the cluster?
            ##quasi-Poisson
            outcluster.bic.qp <- lapply(1:nsim, function(i)
                rrbin_outside%*% matrix(ident.bic.qp[[i]], ncol=1))
            outcluster.aic.qp <- lapply(1:nsim, function(i)
                rrbin_outside%*% matrix(ident.aic.qp[[i]], ncol=1))
            outcluster.aicc.qp <- lapply(1:nsim, function(i)
                rrbin_outside%*% matrix(ident.aicc.qp[[i]], ncol=1))
            ##Poisson
            outcluster.bic.p <- lapply(1:nsim, function(i)
                rrbin_outside%*% matrix(ident.bic.p[[i]], ncol=1))
            outcluster.aic.p <- lapply(1:nsim, function(i)
                rrbin_outside%*% matrix(ident.aic.p[[i]], ncol=1))
            outcluster.aicc.p <- lapply(1:nsim, function(i)
                rrbin_outside%*% matrix(ident.aicc.p[[i]], ncol=1))
            #calc FP rate
            ##quasi-Poisson
            fp.bic.qp <- sum(ifelse(unlist(outcluster.bic.qp)!=0,1,0))/nsim
            fp.aic.qp <- sum(ifelse(unlist(outcluster.aic.qp)!=0,1,0))/nsim
            fp.aicc.qp <- sum(ifelse(unlist(outcluster.aicc.qp)!=0,1,0))/nsim
            outfp.bic.qp <- paste0(fp.bic.qp*100, "%")
            outfp.aic.qp <- paste0(fp.aic.qp*100, "%")
            outfp.aicc.qp <- paste0(fp.aicc.qp*100, "%")
            ##Poisson
            fp.bic.p <- sum(ifelse(unlist(outcluster.bic.p)!=0,1,0))/nsim
            fp.aic.p <- sum(ifelse(unlist(outcluster.aic.p)!=0,1,0))/nsim
            fp.aicc.p <- sum(ifelse(unlist(outcluster.aicc.p)!=0,1,0))/nsim
            outfp.bic.p <- paste0(fp.bic.p*100, "%")
            outfp.aic.p <- paste0(fp.aic.p*100, "%")
            outfp.aicc.p <- paste0(fp.aicc.p*100, "%")

            # ##################################
            # #plot probability maps
            # ##################################
            # #create empties
            # vec <- rep(0, 208 * Time)
            # position.bic <- list(vec)[rep(1, nsim)]
            # position.aic <- list(vec)[rep(1, nsim)]
            # position.aicc <- list(vec)[rep(1, nsim)]
            # #recode identified cells as 1's, all other zeros
            # ix.bic <- lapply(1:nsim, function(i) which(ifelse(ident.bic[[i]]==1,0,1) ==1))
            # ix.aic <- lapply(1:nsim, function(i) which(ifelse(ident.aic[[i]]==1,0,1) ==1))
            # ix.aicc <- lapply(1:nsim, function(i) which(ifelse(ident.aicc[[i]]==1,0,1) ==1))
            # #creatematrix by location (rows) and sim (cols) with 1's indicating selection by superlearner
            # simindicator.bic <- mapply(reval, position.bic, ix.bic)
            # simindicator.aic <- mapply(reval, position.aic, ix.aic)
            # simindicator.aicc <- mapply(reval, position.aicc, ix.aicc)
            # #find probability of detection for each location in time
            # probs.bic <- Matrix::rowSums(simindicator.bic)/nsim
            # probs.aic <- Matrix::rowSums(simindicator.aic)/nsim
            # probs.aicc <- Matrix::rowSums(simindicator.aicc)/nsim
            # #map probability detections by IC to grey scale
            # colprob <- colormapping(list(probs.bic,
            #                              probs.aic,
            #                              probs.aicc,
            #                              as.vector(rrbin_cluster)), Time, cv = NULL,prob=TRUE)
            # #plot map with probability detection by each IC
            # probplotmapAllIC(colprob,genpdf = TRUE, pdfname=paste0(sim.i, "_prob_pc.pdf"))
            # ##################################
            # #RR maps
            # ##################################
            # #plot each RR for each sim            
            # lapply(1:nsim, function(i) plotmapAllIC(res.bic = sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.bic,],
            #                                         res.aic = sim_superclust_pc[[i]]$wLambda[sim_superclust_pc[[i]]$selection.aic,],
            #                                         oracle = rr,
            #                                         genpdf = TRUE,
            #                                         pdfname = paste0(sim.i,"_pc_simID",simid[i],".pdf")))
            ##################################
            #Add sim results to table
            ##################################
            tabn.clusso <- rbind(cbind(IC="BIC",rad, risk, cent, theta,
                                       time=as.numeric(paste(tim, collapse="")),
                                       mod="ST", pow=outpow.bic, fp = outfp.bic),
                                 cbind(IC="AIC",rad, risk, cent, theta,
                                       time=as.numeric(paste(tim, collapse="")),
                                       mod="ST", pow=outpow.aic, fp = outfp.aic),
                                 cbind(IC="AICc",rad, risk, cent, theta,
                                       time=as.numeric(paste(tim, collapse="")),
                                       mod="ST", pow=outpow.aicc, fp = outfp.aicc))
            table.detection.clusso <- rbind(table.detection.clusso, tabn.clusso)
            
            
            
            
        }   
    }
}
}
# Stop the clock
proc.time() - ptm
#####################################################################################
#####################################################################################
#####################################################################################
##WRITE TO CSV
#####################################################################################
#####################################################################################
#####################################################################################
#superclust by loc
print(table.detection.loc)
write.csv(table.detection.loc, file="test_loc.csv", row.names=TRUE)

#superclust by loc
print(table.detection.pc)
write.csv(table.detection.pc, file="test_pc.csv", row.names=TRUE)

#clusso
print(table.detection.clusso)
write.csv(table.detection.clusso, file="test_clusso.csv", row.names=TRUE)
