rm(list=ls())
set.seed(20190326)
library(clusso)
library(MASS)
#source("clustack_20191222.R")
source("../clustack.R")



##################################################
#SINGLE CLUSTER SIMS
##################################################

#0) Setup
dframe1 <- read.csv("../../../../../clusso_KamenetskyLeeZhuGangnon/data/JBC/jap.breast.F.9.10.11.csv")
dframe2 <- read.csv("../../../../../clusso_KamenetskyLeeZhuGangnon/data/JBC/utmJapan.csv")
dframe3 <- aggregate(dframe1, by=list(as.factor(rep(1:(nrow(dframe1)/4),each=4))), FUN="sum")
dframe=data.frame(id=as.factor(dframe3$id/4),period=as.factor(dframe3$year),death=dframe3$death,expdeath=dframe3$expdeath)
levels(dframe$period) <- c("1","2","3","4","5")

dframe.poly2 <- read.csv("../../../../../clusso_KamenetskyLeeZhuGangnon/data/JBC/japan_poly2.csv")
japan.poly2 <- dframe.poly2[,2:3]
dframe.prefect2 <- read.csv("../../../../../clusso_KamenetskyLeeZhuGangnon/data/JBC/japan_prefect2.csv")
japan.prefect2 <- dframe.prefect2[,2:5]
japanbreastcancer <- dframe3
japanbreastcancer$period <- japanbreastcancer$year
japanbreastcancer$year <- NULL

cases <- japanbreastcancer$death
expected <- japanbreastcancer$expdeath

x <- utmJapan$utmx/1000
y <- utmJapan$utmy/1000
japan.poly2 <- dframe.poly2[,2:3]
japan.prefect2 <- dframe.prefect2[,2:5]

#arguments passed
nsim <- 100
nsimstep <- 1000
maxclust <- 15
rMax <- 20
Time <- 5
thetas <- c(Inf, 60)
tim <- c(3:5)
centers <- c(150,35)
radii <- c(9,11,18)
timeperiod <- c(3:5)
#risk.ratios <- c(1.1, 1.5, 2, 0.5)
risk.ratios <- c(1,1.1, 1.5, 2, 0.5)
models <- c("spacetime", "space")
path.figures <- "./"
path.tables <- "./"

eps <- 3

#################################################################################################
#################################################################################################
#Set-Up for Condor
#################################################################################################
#################################################################################################
#nsettings <-96
nsettings <- 120
#2085  = 1040 (Yx) + 1040 (Ex) + cent + rad +risk +theta +mod (each row)
setsims <- matrix(rep(NA, 2085*nsettings*nsim), nrow=nsettings*nsim)
rrsims <- matrix(rep(NA, 1040*nsettings*nsim), nrow=nsettings*nsim)
start <- 0
for (theta in thetas){
    for (cent in centers){
        for (rad in radii){
            for (risk in risk.ratios){
                for (mod in models){
                # print(paste0("Params:", "center: ",cent,"; radius: ",
                #              rad, "; Timeperiods: ", as.numeric(paste(tim, collapse = "")),
                #              "; RR: ", risk, "; Theta :", theta, "; Model :", mod))
                simid <- 1:nsim
                #put the cluster in
                clusters <- clusters2df(x,y,rMax, utm = TRUE, length(x))
                n <- length(x)
                init <- setVectors(japanbreastcancer$period,
                                   japanbreastcancer$expdeath, japanbreastcancer$death,
                                   covars=NULL, Time=Time)
                E1 <- init$E0
                #Ex <- scale(init, Time)
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
                #if (theta=="Inf"){
                if (is.infinite(theta)){
                    #print("Inf theta: pure poisson")
                    YSIM <- lapply(1:nsim, function(i) rpois(length(E1), lambda = E1))
                } else {
                    #print("OverD: NB")
                    YSIM <- lapply(1:nsim, function(i) MASS::rnegbin(E1, theta = theta))
                }
               Ex <- scale_sim(YSIM, init, nsim, Time)
               for (i in 1:nsim){
                   #print(start)
                   #print(i)
                   setsims[i+start,] <- c(YSIM[[i]], Ex[[i]],cent, rad, risk, theta, mod)
                   rrsims[i+start,] <- as.vector(rr)
               }
               start <- start+nsim

                }

            }
        }
    }

}

write.csv(setsims, file="setsims.csv", row.names = FALSE)
write.csv(rrsims, file="rrsims.csv", row.names = FALSE)

# setsims <- read.csv("setsims.csv")
# rrsims <- read.csv("rrsims.csv")