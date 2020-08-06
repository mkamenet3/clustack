#Analyze clustackbounds simulation
#M.Kamenetsky
#2020-07-18
#2020-07-27

#master <- read.csv(master, file="clustackbounds_sim_bypc.csv")
master2 <- read.csv("../../../results/CLUSTACKBOUNDS/clustackbounds_sim_bypc_theta60.csv")
nsim <-100

##########################################################################################
#POST ANALYSIS
##########################################################################################


#need to figure out st locs that are in each cluster I created (should be according to 3 radii (9, 11, 18))
library(MASS)
library(clusso)
library(tidyverse)
source("clustack.R")
load("../data/japanbreastcancer.RData")

#set global
nsim <-100
rMax <- 20 
Time <- 5
x <- utmJapan$utmx/1000
y <- utmJapan$utmy/1000
#create set of potential clusters based on distances
# potentialclusters <- clusters2df(x,y,rMax, utm = TRUE, length(x))
# n_uniq <- length(unique(potentialclusters$center))
# numCenters <- n_uniq

cent <- 150
tim <- c(1:5)

#Sim-Over Params
clusters <- clusters2df(x,y,rMax, utm = TRUE, length(x))
n <- length(x)
tmp <- clusters[clusters$center==cent,]

#9km
cluster9km <- tmp[(tmp$r <= 9),]
rr9km = matrix(1, nrow=n, ncol=Time)
rr9km[cluster9km$last, tim[1]:tail(tim, n=1)] <- NA
clusterix9km = which(is.na(rr9km))
#11km
cluster11km <- tmp[(tmp$r <= 11),]
rr11km = matrix(1, nrow=n, ncol=Time)
rr11km[cluster11km$last, tim[1]:tail(tim, n=1)] <- NA
clusterix11km = which(is.na(rr11km))
#18km
cluster18km <- tmp[(tmp$r <= 18),]
rr18km = matrix(1, nrow=n, ncol=Time)
rr18km[cluster18km$last, tim[1]:tail(tim, n=1)] <- NA
clusterix18km = which(is.na(rr18km))


##############################################
master2$locid <- rep(1:1040, times=nsim)
master2$locidF <- as.factor(rep(1:1040, times=nsim))

master2 <- master2 %>%
    mutate(clusterix = if_else(radius==9 & (locid %in% clusterix9km),1,
                               if_else(radius==11 & (locid %in% clusterix11km),1,
                                       if_else(radius==18 & (locid %in% clusterix18km),1,0)))) 

# masterselects <- master2 %>%
#     dplyr::filter(select.aic!=0 & select.bic!=0) %>%
#     mutate(coverage_buckland.bic = ifelse((risk < adjusted.bic.log.UB & risk > adjusted.bic.log.LB & clusterix==1), 1, 0),
#            coverage_mata.bic = ifelse((risk < mata.bic.UB & risk > mata.bic.LB & clusterix==1), 1, 0),
#            coverage_buckland.aic = ifelse((risk < adjusted.aic.log.UB & risk > adjusted.aic.log.LB & clusterix==1), 1, 0),
#            coverage_mata.aic = ifelse((risk < mata.aic.UB & risk > mata.aic.LB & clusterix==1), 1, 0))

# masterselects <- master2 %>%
#     #dplyr::filter(select.aic!=0 & select.bic!=0) %>%
#     mutate(coverage_buckland.bic = ifelse((risk < adjusted.bic.log.UB & risk > adjusted.bic.log.LB), 1, 0),
#            coverage_mata.bic = ifelse((risk < mata.bic.UB & risk > mata.bic.LB), 1, 0),
#            coverage_buckland.aic = ifelse((risk < adjusted.aic.log.UB & risk > adjusted.aic.log.LB), 1, 0),
#            coverage_mata.aic = ifelse((risk < mata.aic.UB & risk > mata.aic.LB), 1, 0)) %>%
#     #dplyr::select(risk, radius, coverage_mata.bic, coverage_mata.aic, clusterix, simID) %>%
#     group_by(risk, radius) %>%
#     summarize(clustercoverage_mata.bic = mean(ifelse(coverage_mata.bic==1 & clusterix==1,1,0), na.rm=TRUE),
#            clustercoverage_mata.aic = mean(ifelse(coverage_mata.aic==1 & clusterix==1,1,0), na.rm=TRUE))
# masterselects <- master2 %>%
#     mutate(coverage_buckland.bic = ifelse((risk < adjusted.bic.log.UB & risk > adjusted.bic.log.LB), 1, 0),
#            coverage_mata.bic = ifelse((risk < mata.bic.UB & risk > mata.bic.LB), 1, 0),
#            coverage_buckland.aic = ifelse((risk < adjusted.aic.log.UB & risk > adjusted.aic.log.LB), 1, 0),
#            coverage_mata.aic = ifelse((risk < mata.aic.UB & risk > mata.aic.LB), 1, 0)) %>%
#     mutate(cases.bic = case_when(coverage_mata.bic==1 & clusterix==1 ~ 1,
#                                  coverage_mata.bic==1 & clusterix==0 ~ 0,
#                                                 coverage_mata.bic==0 & clusterix==0 ~ 0,
#                                                 coverage_mata.bic==0 & clusterix==1 ~ 0,
#                                                 is.na(coverage_mata.bic) & clusterix==0 ~ 0,
#                                                 is.na(coverage_mata.bic) & clusterix==1 ~ 0)) %>%
#     group_by(risk, radius, simID) %>%
#     mutate(test = ifelse(any(cases.bic==1),1,0))


master<- master2 %>%
    mutate(coverage_buckland.bic = ifelse((risk < adjusted.bic.log.UB & risk > adjusted.bic.log.LB), 1, 0),
           coverage_mata.bic = ifelse((risk < mata.bic.UB & risk > mata.bic.LB), 1, 0),
           coverage_buckland.aic = ifelse((risk < adjusted.aic.log.UB & risk > adjusted.aic.log.LB), 1, 0),
           coverage_mata.aic = ifelse((risk < mata.aic.UB & risk > mata.aic.LB), 1, 0))
    
mastercoverage <- master
    group_by(risk, radius, simID) %>%
    mutate(Indclustercoverage_mata.bic = ifelse(any(coverage_mata.bic==1 & clusterix==1),1,0),
           Indclustercoverage_mata.aic = ifelse(any(coverage_mata.aic==1 & clusterix==1),1,0)) %>%
    replace_na(list(Indclustercoverage_mata.bic = 0, Indclustercoverage_mata.aic = 0)) %>%
    group_by(risk, radius) %>%
    summarise(clustercoverage_mata.bic = mean(Indclustercoverage_mata.bic),
              clustercoverage_mata.aic = mean(Indclustercoverage_mata.aic)) %>%
    pivot_longer(cols=clustercoverage_mata.bic:clustercoverage_mata.aic, names_to="IC", values_to="prob") %>%
    mutate(IC = ifelse(IC=="clustercoverage_mata.aic", "QAIC", "QBIC"),
           method = "MATA")

mastercoverage_select <- master2 %>%
    mutate(coverage_buckland.bic = ifelse((risk < adjusted.bic.log.UB & risk > adjusted.bic.log.LB), 1, 0),
           coverage_mata.bic = ifelse((risk < mata.bic.UB & risk > mata.bic.LB), 1, 0),
           coverage_buckland.aic = ifelse((risk < adjusted.aic.log.UB & risk > adjusted.aic.log.LB), 1, 0),
           coverage_mata.aic = ifelse((risk < mata.aic.UB & risk > mata.aic.LB), 1, 0)) %>%
    dplyr::filter(select.aic!=0 & select.bic!=0) %>%
    group_by(risk, radius, simID) %>%
    mutate(Indclustercoverage_mata.bic = ifelse(any(coverage_mata.bic==1 & clusterix==1),1,0),
           Indclustercoverage_mata.aic = ifelse(any(coverage_mata.aic==1 & clusterix==1),1,0)) %>%
    replace_na(list(Indclustercoverage_mata.bic = 0, Indclustercoverage_mata.aic = 0)) %>%
    group_by(risk, radius) %>%
    summarise(clustercoverage_mata.bic = mean(Indclustercoverage_mata.bic),
              clustercoverage_mata.aic = mean(Indclustercoverage_mata.aic)) %>%
    pivot_longer(cols=clustercoverage_mata.bic:clustercoverage_mata.aic, names_to="IC", values_to="prob") %>%
    mutate(IC = ifelse(IC=="clustercoverage_mata.aic", "QAIC", "QBIC"),
           method = "MATA")

mastercoverage_selectK <- master2 %>%
    mutate(coverage_buckland.bic = ifelse((risk < adjusted.bic.log.UB & risk > adjusted.bic.log.LB), 1, 0),
           coverage_mata.bic = ifelse((risk < mata.bic.UB & risk > mata.bic.LB), 1, 0),
           coverage_buckland.aic = ifelse((risk < adjusted.aic.log.UB & risk > adjusted.aic.log.LB), 1, 0),
           coverage_mata.aic = ifelse((risk < mata.aic.UB & risk > mata.aic.LB), 1, 0)) %>%
    dplyr::filter(select.aic!=0 & select.bic!=0) %>%
    group_by(risk, radius, simID) %>%
    mutate(Indclustercoverage_mata.bic = ifelse(any(coverage_mata.bic==1 & clusterix==1),1,0),
           Indclustercoverage_mata.aic = ifelse(any(coverage_mata.aic==1 & clusterix==1),1,0)) %>%
    replace_na(list(Indclustercoverage_mata.bic = 0, Indclustercoverage_mata.aic = 0)) %>%
    group_by(risk, radius, select.bic, select.aic) %>%
    summarise(clustercoverage_mata.bic = mean(Indclustercoverage_mata.bic),
              clustercoverage_mata.aic = mean(Indclustercoverage_mata.aic))%>%
    pivot_longer(cols=clustercoverage_mata.bic:clustercoverage_mata.aic, names_to="IC", values_to="prob") %>%
    mutate(IC = ifelse(IC=="clustercoverage_mata.aic", "QAIC", "QBIC"),
           method = "MATA") %>%
    pivot_longer(select.bic:select.aic, names_to="ICselect", values_to="K") %>%
    dplyr::select(-IC)



################################
#Coverage probability: by simulation
#if at least 1 cell was within the bounds 
################################
ggplot(data=mastercoverage_select, aes(x=risk, y=prob)) +
    geom_line(col="blue", size=2) +
    facet_grid(radius~IC) +
    xlab("Relative risk (RR)") +
    ylab("MATA-Interval Coverage Probability") +
    theme_bw()+
    ggtitle("MATA-Interval Coverage") +
    scale_x_continuous(breaks=c(1.1, 1.5,2)) +
    ylim(0,1)+
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
    labs(caption="Coverage probability: at least 1 ST cell was within the bounds. \nAll simulations considered (including those with no detection)")
    
ggplot(data=mastercoverage_select, aes(x=risk, y=prob)) +
    geom_line(col="blue", size=2) +
    facet_grid(radius~IC) +
    xlab("Relative risk (RR)") +
    ylab("MATA-Interval Coverage Probability") +
    theme_bw()+
    ggtitle("MATA-Interval Coverage") +
    scale_x_continuous(breaks=c(1.1, 1.5,2)) +
    ylim(0,1)+
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
    labs(caption="Coverage probability: at least 1 ST cell was within the bounds. \nExcluding sims with no detection.")

mastercoverage_selectK %>%
    dplyr::filter(ICselect=="select.bic") %>%
    ggplot(data=, aes(x=risk, y=prob)) +
    geom_point(col="blue", size=2) +
    facet_grid(radius~K)


    xlab("Relative risk (RR)") +
    ylab("MATA-Interval Coverage Probability") +
    theme_bw()+
    ggtitle("MATA-Interval Coverage") +
    scale_x_continuous(breaks=c(1.1, 1.5,2)) +
    ylim(0,1)+
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
    labs(caption="Coverage probability: at least 1 ST cell was within the bounds. \nExcluding sims with no detection.")

################################
#Coverage probability: by location 
################################
masterlocs <- master %>%
    dplyr::filter(clusterix==1) %>%
    dplyr::filter(select.aic!=0 & select.bic!=0) %>%
    arrange(thetaa.bic) %>%
    droplevels() %>%
    group_by(locid, risk, radius) %>%
    summarize(thetaa.bic.mean = mean(thetaa.bic),
              thetaa.aic.mean = mean(thetaa.aic)) %>%
    pivot_longer(cols=thetaa.bic.mean:thetaa.aic.mean, names_to="IC", values_to="thetaa")%>%
    mutate(IC = ifelse(IC=="thetaa.aic.mean", "QAIC", "QBIC"))

# masterlocs_all <- master %>%
#         #dplyr::filter(clusterix==1) %>%
#         dplyr::filter(select.aic!=0 & select.bic!=0) %>%
#         arrange(thetaa.bic) %>%
#         droplevels() %>%
#         group_by(locid, risk, radius, clusterix) %>%
#         summarize(thetaa.bic.mean = mean(thetaa.bic),
#                   thetaa.aic.mean = mean(thetaa.aic), .groups="keep") %>%
#         pivot_longer(cols=thetaa.bic.mean:thetaa.aic.mean, names_to="IC", values_to="thetaa")%>%
#         mutate(IC = ifelse(IC=="thetaa.aic.mean", "QAIC", "QBIC"))        

    masterlocs_all <- master %>%
        #dplyr::filter(clusterix==1) %>%
        dplyr::filter(select.aic!=0 & select.bic!=0) %>%
        arrange(thetaa.bic) %>%
        droplevels() %>%
        group_by(locid, risk, radius, clusterix) %>%
        summarize(thetaa.bic.mean = mean(thetaa.bic),
                  thetaa.aic.mean = mean(thetaa.aic),
                  mata.bic.UB.mean = mean(mata.bic.UB, na.rm=TRUE),
                  mata.bic.LB.mean = mean(mata.bic.LB, na.rm=TRUE),
                  mata.aic.UB.mean = mean(mata.aic.UB, na.rm=TRUE),
                  mata.aic.LB.mean = mean(mata.aic.LB, na.rm=TRUE),.groups="keep") %>%
        pivot_longer(cols=thetaa.bic.mean:thetaa.aic.mean, names_to="IC", values_to="thetaa")%>%
        #pivot_longer(cols=thetaa.bic.mean:mata.bic.UB.mean, names_to="IC", values_to="thetaa")%>%
        mutate(IC = ifelse(IC=="thetaa.aic.mean", "QAIC", "QBIC")) %>%
        pivot_longer(cols=c(mata.bic.UB.mean,mata.aic.UB.mean), names_to="mataUB", values_to="thetamata.UB") %>%
        pivot_longer(cols=c(mata.bic.LB.mean,mata.aic.LB.mean), names_to="mataLB", values_to="thetamata.LB") %>%
        mutate(mataUB = ifelse(mataUB=="mata.aic.UB.mean", "QAIC", "QBIC"),
               mataLB = ifelse(mataLB=="mata.aic.LB.mean", "QAIC", "QBIC")) %>%
        dplyr::filter(IC == mataUB & IC==mataLB)
    
clusters_uniq <- clusters %>%
    dplyr::select(center, x,y) %>%
    dplyr::distinct()
#masterlocs_dist <- merge(masterlocs, clusters_uniq, by.x="locid", by.y="center", all.x=TRUE)
distances <- as.matrix(dist(cbind(clusters_uniq$x, clusters_uniq$y)))[,cent]
distancesfromcenter <- cbind.data.frame(distfromcent=rep(distances, times=Time),
                                        center = rep(as.numeric(names(distances)), times=Time),
                                        Period = rep(c("Period 1","Period 2","Period 3","Period 4","Period 5"),each=208),
                                        centerTime = 1:1040)
masterlocs_dist <- merge(masterlocs, distancesfromcenter, by.x="locid", by.y="centerTime", all.x=TRUE)
masterlocs_all_dist <- merge(masterlocs_all, distancesfromcenter, by.x="locid", by.y="centerTime", all.x=TRUE)
masterlocs_all_dist <- masterlocs_all_dist %>%
    dplyr::mutate(radius2=ifelse(clusterix==1,radius, "Outside of Cluster"))

    
    
masterlocsK <- master %>%
        dplyr::filter(clusterix==1) %>%
        dplyr::filter(select.aic!=0 & select.bic!=0) %>%
        arrange(thetaa.bic) %>%
        droplevels() %>%
        group_by(locidF, risk, radius,select.bic, select.aic) %>%
        summarize(thetaa.bic.mean = mean(thetaa.bic),
                  thetaa.aic.mean = mean(thetaa.aic)) %>%
        pivot_longer(cols=thetaa.bic.mean:thetaa.aic.mean, names_to="IC", values_to="thetaa")%>%
        mutate(IC = ifelse(IC=="thetaa.aic.mean", "QAIC", "QBIC")) %>%
        pivot_longer(select.bic:select.aic, names_to="ICselect", values_to="K")%>%
        dplyr::select(-IC) %>%
        distinct()
    
#this plot is good
ggplot() +
    geom_line(data=masterlocs, aes(x=fct_inorder(locidF), y=thetaa, group=radius, color=factor(radius)), size=1.5, alpha=0.9)+
    facet_grid(IC ~ risk) +
    theme_bw() +
    scale_x_discrete(breaks=0) +
    ylab("Mean Thetaa") +
    xlab("Cell i at Time t in cluster C") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
    labs(title="Mean Model-Averaged thetaa Estimates by Cell i, Time t",
         caption="Only ST cells inside each respective cluster.\n Mean thetaa calculted across 100 simulations for radius-risk pair. \nIt looks like as cluster increases, \nthe mean thetaa shrinks due to more cells being in the cluster and the relative risk weights dispersed more.\nConditional on a cluster being selected")

#####################################################
#let x by spatial distance form center
ggplot() +
    geom_line(data=subset(masterlocs_all_dist, clusterix==0), aes(x=distfromcent, 
                                                                  y=thetaa),
              color="black", fill="black",
              size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1), aes(x=distfromcent, 
                                            y=thetaa, group=radius, color=factor(radius)), size=1, alpha=0.8) +
    facet_grid(IC ~ risk) +
    theme_bw() +
    ylab(bquote(bar(theta)[a])) +
    xlab("Distance from cluster center (km)") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
    labs(title=bquote(bar(theta)[a] ~  "Estimates by Distance from Cluster Center"),
         caption="Only ST cells inside each respective cluster.\n Mean thetaa calculted across 100 simulations for radius-risk pair. \nIt looks like as cluster increases, \nthe mean thetaa shrinks due to more cells being in the cluster and the relative risk weights dispersed more.\nConditional on a cluster being selected",
         color = "Cluster radius (km)")

# ggplot() +
#     geom_line(data=subset(masterlocs_all_dist, clusterix==0 & IC=="QBIC"), aes(x=distfromcent, 
#                                                                   y=thetaa),
#               color="black",
#               size=1, alpha=0.8) +
#     geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC"), aes(x=distfromcent, 
#                                                                   y=thetaa, 
#                                                                   group=radius, color=factor(radius)), size=1, alpha=0.8) +
#     facet_grid(factor(risk) ~ Period) +
#     theme_minimal() +
#     ylab(bquote(bar(theta)[a])) +
#     xlab("Distance from cluster center (km)") +
#     theme(plot.title = element_text(hjust = 0.5),
#           plot.caption = element_text(face="italic")) +
#     labs(title=bquote(bar(theta)[a] ~  "Estimates by Distance from Cluster Center - QBIC"),
#          caption="Only ST cells inside each respective cluster.\n Mean thetaa calculted across 100 simulations for radius-risk pair. \nIt looks like as cluster increases, \nthe mean thetaa shrinks due to more cells being in the cluster and the relative risk weights dispersed more.\nConditional on a cluster being selected",
#          color = "Cluster radius (km)") +
#     scale_y_continuous(
#         # Add a second axis and specify its features
#         sec.axis = sec_axis(~., name="Cluster Relative Risk", breaks=0) 
#    ) 
cols <- c("0" = "red", "9" = "blue", "11" = "darkgreen", "18" = "orange")
masterlocs_all_dist_maxthetaa <- masterlocs_all_dist %>%
    group_by(radius, risk, IC, Period) %>%
    dplyr::filter(thetaa== max(thetaa))
masterlocs_all_dist_maxthetaa_noclust <- masterlocs_all_dist %>%
    dplyr::filter(clusterix==0) %>%
    group_by(risk, IC, Period) %>%
    dplyr::filter(thetaa== max(thetaa))
    #summarize(maxthetaa = max(thetaa, na.rm=TRUE),.groups="keep")
pa <- ggplot() +
    geom_line(data=subset(masterlocs_all_dist, clusterix==0 & IC=="QBIC"), 
              aes(x=distfromcent,  y=thetaa, color="Outside of Cluster"),
             # color="black",
              size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==9), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="9 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==11), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="11 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==18), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="18 km"), size=1, alpha=0.8) +
    geom_point(data=subset(masterlocs_all_dist_maxthetaa_noclust, IC=="QBIC"), 
               aes(x=distfromcent, y=thetaa, color="Outside of Cluster"),shape=8) + 
    geom_point(data=subset(masterlocs_all_dist_maxthetaa, IC=="QBIC" & radius ==9), 
               aes(x=distfromcent, y=thetaa, color="9 km"),shape=8) + 
    geom_point(data=subset(masterlocs_all_dist_maxthetaa, IC=="QBIC" & radius ==11), 
               aes(x=distfromcent, y=thetaa, color="11 km"),shape=8) + 
    geom_point(data=subset(masterlocs_all_dist_maxthetaa, IC=="QBIC" & radius ==18), 
               aes(x=distfromcent, y=thetaa, color="18 km"),shape=8) + 
    facet_grid(factor(risk) ~ Period) +
    theme_minimal() +
    ylab(bquote(bar(theta)[a])) +
    xlab("Distance from cluster center (km)") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
    labs(title=bquote(bar(theta)[a] ~  "Estimates by Distance from Cluster Center - QBIC"),
         color = "Cluster radius (km)") +
    scale_y_continuous(sec.axis = sec_axis(~., name="Cluster Relative Risk", breaks=0)) +
    scale_color_manual(name="Cluster radius", values=c("Outside of Cluster"="grey20",
                                         "9 km" = "red",
                                         "11 km" = "orange",
                                         "18 km" = "blue"),
                       breaks=c("9 km", "11 km", "18 km", "Outside of Cluster")) +
    geom_hline(yintercept = 1, linetype="dashed") 
pa
ggsave(file="../../../figures/SimulationResults/qbic_meanthetaa_sim.pdf", height=8, width=12, units="in")

pb <- ggplot() +
    geom_line(data=subset(masterlocs_all_dist, clusterix==0 & IC=="QAIC"), 
              aes(x=distfromcent,  y=thetaa, color="Outside of Cluster"),
              # color="black",
              size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==9), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="9 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==11), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="11 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==18), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="18 km"), size=1, alpha=0.8) +
    geom_point(data=subset(masterlocs_all_dist_maxthetaa_noclust, IC=="QAIC"), 
               aes(x=distfromcent, y=thetaa, color="Outside of Cluster"),shape=8) + 
    geom_point(data=subset(masterlocs_all_dist_maxthetaa, IC=="QAIC" & radius ==9), 
               aes(x=distfromcent, y=thetaa, color="9 km"),shape=8) + 
    geom_point(data=subset(masterlocs_all_dist_maxthetaa, IC=="QAIC" & radius ==11), 
               aes(x=distfromcent, y=thetaa, color="11 km"),shape=8) + 
    geom_point(data=subset(masterlocs_all_dist_maxthetaa, IC=="QAIC" & radius ==18), 
               aes(x=distfromcent, y=thetaa, color="18 km"),shape=8) + 
    facet_grid(factor(risk) ~ Period) +
    theme_minimal() +
    ylab(bquote(bar(theta)[a])) +
    xlab("Distance from cluster center (km)") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
    labs(title=bquote(bar(theta)[a] ~  "Estimates by Distance from Cluster Center - QAIC"),
         color = "Cluster radius (km)") +
    scale_y_continuous(sec.axis = sec_axis(~., name="Cluster Relative Risk", breaks=0)) +
    scale_color_manual(name="Cluster radius", values=c("Outside of Cluster"="grey20",
                                                       "9 km" = "red",
                                                       "11 km" = "orange",
                                                       "18 km" = "blue"),
                       breaks=c("9 km", "11 km", "18 km", "Outside of Cluster")) +
    geom_hline(yintercept = 1, linetype="dashed") 
pb
ggsave(file="../../../figures/SimulationResults/qaic_meanthetaa_sim.pdf", height=8, width=12, units="in")

#combine
paa <- ggplot() +
    geom_line(data=subset(masterlocs_all_dist, clusterix==0 & IC=="QBIC"), 
              aes(x=distfromcent,  y=thetaa, color="Outside of Cluster"),
              # color="black",
              size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==9), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="9 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==11), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="11 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==18), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="18 km"), size=1, alpha=0.8) +
    geom_point(data=subset(masterlocs_all_dist_maxthetaa_noclust, IC=="QBIC"), 
               aes(x=distfromcent, y=thetaa, color="Outside of Cluster"),shape=8) + 
    geom_point(data=subset(masterlocs_all_dist_maxthetaa, IC=="QBIC" & radius ==9), 
               aes(x=distfromcent, y=thetaa, color="9 km"),shape=8) + 
    geom_point(data=subset(masterlocs_all_dist_maxthetaa, IC=="QBIC" & radius ==11), 
               aes(x=distfromcent, y=thetaa, color="11 km"),shape=8) + 
    geom_point(data=subset(masterlocs_all_dist_maxthetaa, IC=="QBIC" & radius ==18), 
               aes(x=distfromcent, y=thetaa, color="18 km"),shape=8) + 
    facet_grid(factor(risk) ~ Period) +
    theme_minimal() +
    ylab(bquote(bar(theta)[a])) +
    xlab("") +
    #xlab("Distance from cluster center (km)") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
    # labs(title=bquote(bar(theta)[a] ~  "Estimates by Distance from Cluster Center - QBIC"),
    #      color = "Cluster radius (km)") +
    scale_y_continuous(sec.axis = sec_axis(~., name="Cluster Relative Risk", breaks=0)) +
    scale_color_manual(name="Cluster radius", values=c("Outside of Cluster"="grey20",
                                                       "9 km" = "red",
                                                       "11 km" = "orange",
                                                       "18 km" = "blue"),
                       breaks=c("9 km", "11 km", "18 km", "Outside of Cluster")) +
    geom_hline(yintercept = 1, linetype="dashed") 

pbb <- ggplot() +
    geom_line(data=subset(masterlocs_all_dist, clusterix==0 & IC=="QAIC"), 
              aes(x=distfromcent,  y=thetaa, color="Outside of Cluster"),
              # color="black",
              size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==9), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="9 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==11), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="11 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==18), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="18 km"), size=1, alpha=0.8) +
    geom_point(data=subset(masterlocs_all_dist_maxthetaa_noclust, IC=="QAIC"), 
               aes(x=distfromcent, y=thetaa, color="Outside of Cluster"),shape=8) + 
    geom_point(data=subset(masterlocs_all_dist_maxthetaa, IC=="QAIC" & radius ==9), 
               aes(x=distfromcent, y=thetaa, color="9 km"),shape=8) + 
    geom_point(data=subset(masterlocs_all_dist_maxthetaa, IC=="QAIC" & radius ==11), 
               aes(x=distfromcent, y=thetaa, color="11 km"),shape=8) + 
    geom_point(data=subset(masterlocs_all_dist_maxthetaa, IC=="QAIC" & radius ==18), 
               aes(x=distfromcent, y=thetaa, color="18 km"),shape=8) + 
    facet_grid(factor(risk) ~ Period) +
    theme_minimal() +
    ylab(bquote(bar(theta)[a])) +
    xlab("Distance from cluster center (km)") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
    # labs(title=bquote(bar(theta)[a] ~  "Estimates by Distance from Cluster Center - QAIC"),
    #      color = "Cluster radius (km)") +
    scale_y_continuous(sec.axis = sec_axis(~., name="Cluster Relative Risk", breaks=0)) +
    scale_color_manual(name="Cluster radius", values=c("Outside of Cluster"="grey20",
                                                       "9 km" = "red",
                                                       "11 km" = "orange",
                                                       "18 km" = "blue"),
                       breaks=c("9 km", "11 km", "18 km", "Outside of Cluster")) +
    geom_hline(yintercept = 1, linetype="dashed")

plotrow <- cowplot::plot_grid(paa+ theme(legend.position="none"),
                              pbb+ theme(legend.position="none"),
                              nrow=2, labels=c("QBIC", "QAIC"))
title <- ggdraw() + 
    draw_label(label=bquote(bar(theta)[a] ~  "Estimates by Distance from Cluster Center"),
        hjust = 0.5) 
# extract the legend from one of the plots
legend <- get_legend(
    # create some space to the left of the legend
    paa + theme(legend.box.margin = margin(0, 0, 0, 12))
)
out1 <- plot_grid(
    title, plotrow,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
)
plot_grid(out1, legend, rel_widths = c(3, .4))
ggsave(file="../../../figures/SimulationResults/qbicqaic_meanthetaa_sim.pdf", height=14, width=14, units="in")



#################################
#Plotting MATA BOUNDS
ggplot() +
    geom_line(data=subset(masterlocs_all_dist, clusterix==0 & IC=="QAIC"), 
              aes(x=distfromcent,  y=thetaa, color="Outside of Cluster"),
              # color="black",
              size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==9), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="9 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==11), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="11 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==18), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="18 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==0 & IC=="QAIC"), 
              aes(x=distfromcent,  y=thetamata.UB, color="Outside of Cluster"),
              # color="black",
              size=0.6, linetype="dotted") +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==9), 
              aes(x=distfromcent, y=thetamata.UB, 
                  group=radius, color="9 km"), size=0.6,  linetype="dotted") +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==11), 
              aes(x=distfromcent, y=thetamata.UB, 
                  group=radius, color="11 km"), size=0.6, linetype="dotted") +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==18), 
              aes(x=distfromcent, y=thetamata.UB, 
                  group=radius, color="18 km"), size=0.6, linetype="dotted") +
    #lower
    geom_line(data=subset(masterlocs_all_dist, clusterix==0 & IC=="QAIC"), 
              aes(x=distfromcent,  y=thetamata.LB, color="Outside of Cluster"),
              # color="black",
              size=0.6, linetype="dotted") +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==9), 
              aes(x=distfromcent, y=thetamata.LB, 
                  group=radius, color="9 km"), size=0.6,  linetype="dotted") +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==11), 
              aes(x=distfromcent, y=thetamata.LB, 
                  group=radius, color="11 km"), size=0.6, linetype="dotted") +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QAIC" & radius==18), 
              aes(x=distfromcent, y=thetamata.LB, 
                  group=radius, color="18 km"), size=0.6, linetype="dotted") +
    facet_grid(factor(risk) ~ Period) +
    theme_minimal() +
    ylab(bquote(bar(theta)[a])) +
    xlab("Distance from cluster center (km)") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
     labs(title=bquote(bar(theta)[a] ~  "Estimates by Distance from Cluster Center - QAIC"),
          color = "Cluster radius (km)") +
    scale_y_continuous(sec.axis = sec_axis(~., name="Cluster Relative Risk", breaks=0)) +
    scale_color_manual(name="Cluster radius", values=c("Outside of Cluster"="grey20",
                                                       "9 km" = "red",
                                                       "11 km" = "orange",
                                                       "18 km" = "blue"),
                       breaks=c("9 km", "11 km", "18 km", "Outside of Cluster")) +
    geom_hline(yintercept = 1, linetype="dashed")
ggsave(file="../../../figures/SimulationResults/qaic_meanthetaa_bounds_sim.pdf", height=8, width=12, units="in")


ggplot() +
    geom_line(data=subset(masterlocs_all_dist, clusterix==0 & IC=="QBIC"), 
              aes(x=distfromcent,  y=thetaa, color="Outside of Cluster"),
              # color="black",
              size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==9), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="9 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==11), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="11 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==18), 
              aes(x=distfromcent, y=thetaa, 
                  group=radius, color="18 km"), size=1, alpha=0.8) +
    geom_line(data=subset(masterlocs_all_dist, clusterix==0 & IC=="QBIC"), 
              aes(x=distfromcent,  y=thetamata.UB, color="Outside of Cluster"),
              # color="black",
              size=0.6, linetype="dotted") +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==9), 
              aes(x=distfromcent, y=thetamata.UB, 
                  group=radius, color="9 km"), size=0.6,  linetype="dotted") +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==11), 
              aes(x=distfromcent, y=thetamata.UB, 
                  group=radius, color="11 km"), size=0.6, linetype="dotted") +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==18), 
              aes(x=distfromcent, y=thetamata.UB, 
                  group=radius, color="18 km"), size=0.6, linetype="dotted") +
    #lower
    geom_line(data=subset(masterlocs_all_dist, clusterix==0 & IC=="QBIC"), 
              aes(x=distfromcent,  y=thetamata.LB, color="Outside of Cluster"),
              # color="black",
              size=0.6, linetype="dotted") +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==9), 
              aes(x=distfromcent, y=thetamata.LB, 
                  group=radius, color="9 km"), size=0.6,  linetype="dotted") +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==11), 
              aes(x=distfromcent, y=thetamata.LB, 
                  group=radius, color="11 km"), size=0.6, linetype="dotted") +
    geom_line(data=subset(masterlocs_all_dist, clusterix==1 & IC=="QBIC" & radius==18), 
              aes(x=distfromcent, y=thetamata.LB, 
                  group=radius, color="18 km"), size=0.6, linetype="dotted") +
    facet_grid(factor(risk) ~ Period) +
    theme_minimal() +
    ylab(bquote(bar(theta)[a])) +
    xlab("Distance from cluster center (km)") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
    labs(title=bquote(bar(theta)[a] ~  "Estimates by Distance from Cluster Center - QBIC"),
         color = "Cluster radius (km)") +
    scale_y_continuous(sec.axis = sec_axis(~., name="Cluster Relative Risk", breaks=0)) +
    scale_color_manual(name="Cluster radius", values=c("Outside of Cluster"="grey20",
                                                       "9 km" = "red",
                                                       "11 km" = "orange",
                                                       "18 km" = "blue"),
                       breaks=c("9 km", "11 km", "18 km", "Outside of Cluster")) +
    geom_hline(yintercept = 1, linetype="dashed")
ggsave(file="../../../figures/SimulationResults/qbic_meanthetaa_bounds_sim.pdf", height=8, width=12, units="in")



########################################################
#Plotting all sims by bounds

ggplot() +
    geom_line(data=subset(masterlocs, IC=="QAIC"), aes(x=fct_inorder(locidF), y=thetaa, group=radius),color="black", size=1.5)+
    facet_grid(radius ~ risk) +
    theme_bw() +
    scale_x_discrete(breaks=0) +
    ylab("Mean Thetaa by QAIC") +
    xlab("Cell i at Time t in cluster C") +
    geom_line(data=masterlocsMATA, aes(x=fct_inorder(locidF), y=mata.aic.LB, group=simID, color=factor(radius)), size=0.2, alpha=0.3) +
    geom_line(data=masterlocsMATA, aes(x=fct_inorder(locidF), y=mata.aic.UB, group=simID, color=factor(radius)), size=0.2, alpha=0.3) +
    geom_hline(data=subset(masterlocs, IC=="QAIC"),aes(yintercept = risk), size=1, linetype="dashed") 

########################################################
pb1 <- masterlocsK %>%
    dplyr::filter(ICselect=="select.bic") %>%
ggplot() +
    geom_line(aes(x=fct_inorder(locidF), y=thetaa, group=radius, color=factor(radius)), size=1.5, alpha=0.6)+
    facet_grid(K ~ risk) +
    theme_bw() +
    scale_x_discrete(breaks=0) +
    ylab("Mean Thetaa by BIC") +
    xlab("Cell i at Time t in cluster C") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
    labs(title="Mean Model-Averaged thetaa Estimates by Cell i, Time t - QBIC",
         caption="Only ST cells inside each respective cluster.\n Mean thetaa calculted across 100 simulations for radius-risk pair. \nIt looks like as cluster increases, \nthe mean thetaa shrinks due to more cells being in the cluster and the relative risk weights dispersed more.\nConditional on a cluster being selected") +
    scale_y_continuous(
        # Features of the first axis
        name = "Mean Thetaa",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~., name="K Clusters")
    ) +
    geom_hline(aes(yintercept = risk), size=1, linetype="dashed") 

pa1 <- masterlocsK %>%
    dplyr::filter(ICselect=="select.aic") %>%
    ggplot() +
    geom_line(aes(x=fct_inorder(locidF), y=thetaa, group=radius, color=factor(radius)), size=1.5, alpha=0.6)+
    facet_grid(K ~ risk) +
    theme_bw() +
    scale_x_discrete(breaks=0) +
    ylab("Mean Thetaa by BIC") +
    xlab("Cell i at Time t in cluster C") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
    labs(title="Mean Model-Averaged thetaa Estimates by Cell i, Time t - QAIC",
         caption="Only ST cells inside each respective cluster.\n Mean thetaa calculted across 100 simulations for radius-risk pair. \nIt looks like as cluster increases, \nthe mean thetaa shrinks due to more cells being in the cluster and the relative risk weights dispersed more.\nConditional on a cluster being selected") +
    scale_y_continuous(
        # Features of the first axis
        name = "Mean Thetaa",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~., name="K Clusters")
    ) +
    geom_hline(aes(yintercept = risk), size=1, linetype="dashed") 

# masterlocsK %>%
#     dplyr::filter(ICselect=="select.aic") %>%
    ggplot(data=masterlocsK) +
    geom_point( aes(x=thetaa, y=K, group=radius, color=factor(radius)), size=1.5, alpha=0.6)+
    geom_jitter(aes(x=thetaa, y=K, group=radius, color=factor(radius)), width=0) +
    geom_vline(aes(xintercept = risk), size=1, linetype="dashed") +
    facet_grid(ICselect ~ risk) +
     theme_bw() +
    ylab("K identified clusters") +
    xlab("Model-Averaged Thetaa") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(face="italic")) +
    ggtitle("Thetaa vs. Number of Clusters, RRR, and selection crit")
    


masterlocsMATA <- master %>%
    dplyr::filter(clusterix==1) %>%
    dplyr::filter(select.aic!=0 & select.bic!=0) %>%
    arrange(thetaa.bic) %>%
    droplevels() 
    
p1 <- ggplot() +
    geom_line(data=subset(masterlocs, IC=="QAIC"), aes(x=fct_inorder(locidF), y=thetaa, group=radius),color="black", size=1.5)+
    facet_grid(radius ~ risk) +
    theme_bw() +
    scale_x_discrete(breaks=0) +
    ylab("Mean Thetaa by QAIC") +
    xlab("Cell i at Time t in cluster C") +
    geom_line(data=masterlocsMATA, aes(x=fct_inorder(locidF), y=mata.aic.LB, group=simID, color=factor(radius)), size=0.2, alpha=0.3) +
    geom_line(data=masterlocsMATA, aes(x=fct_inorder(locidF), y=mata.aic.UB, group=simID, color=factor(radius)), size=0.2, alpha=0.3) +
    geom_hline(data=subset(masterlocs, IC=="QAIC"),aes(yintercept = risk), size=1, linetype="dashed") 
        
p1

p2 <- ggplot() +
    geom_line(data=subset(masterlocs, IC=="QBIC"), aes(x=fct_inorder(locidF), y=thetaa, group=radius),color="black", size=1.5)+
    facet_grid(radius ~ risk) +
    theme_bw() +
    scale_x_discrete(breaks=0) +
    ylab("Mean Thetaa by QBIC") +
    xlab("Cell i at Time t in cluster C") +
    geom_line(data=masterlocsMATA, aes(x=fct_inorder(locidF), y=mata.bic.LB, group=simID, color=factor(radius)), size=0.2, alpha=0.3) +
    geom_line(data=masterlocsMATA, aes(x=fct_inorder(locidF), y=mata.bic.UB, group=simID, color=factor(radius)), size=0.2, alpha=0.3) +
    geom_hline(data=subset(masterlocs, IC=="QAIC"),aes(yintercept = risk), size=1, linetype="dashed") 

gridExtra::grid.arrange(p1, p2, ncol=2, left="QAIC", right="QBIC", bottom="MATA-based Intervals")


gridExtra::grid.arrange(pa1, pb1, ncol=2, left="QAIC", right="QBIC", bottom="MATA-based Intervals")

##################################################################################################################
##################################################################################################################
##################################################################################################################

################################
#Coverage probability: given that a cluster was selected (aic and bic selects !=0),
#what is the coverage probability?
################################
masterselects %>%
    dplyr::filter(risk==2 & radius==18) %>%
    arrange(thetaa.aic) %>%
    ggplot(aes(x=fct_inorder(locidF), y=mata.aic.UB)) +
    geom_line()
    #geom_ribbon(aes(ymin=mata.aic.LB, ymax=mata.aic.UB), alpha=0.1) 

masterselects %>%
    dplyr::filter(risk==2 & radius==18) %>%
    dplyr::filter(clusterix==1) %>%
    arrange(thetaa.bic) %>%
    droplevels() %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.bic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = mata.bic.LB, ymax = mata.bic.UB), alpha=0.5, size=1) 
    # geom_hline(yintercept = 2, color="red") +
    # ggtitle("estimates & bounds inside RR=2 cluster (MATA, BIC)")

coverage <- masterselects %>%
    dplyr::filter(clusterix==1) %>%
    dplyr::select(locid,risk, radius, coverage_mata.bic, coverage_mata.aic) %>%
    #View()
    group_by(locid, radius, risk) %>%
    summarise(coverage.mata.aic = mean(coverage_mata.aic), 
              coverage.mata.bic = mean(coverage_mata.bic),
              .groups="keep") 

coverage %>%
    arrange(locid) %>%
    ggplot(aes(y=coverage.mata.aic, x=risk)) +
        geom_line() +
    facet_grid(.~radius)



#####################################################################################################################
#####################################################################################################################
#####################################################################################################################





ggplot(data=master2,aes(x=locid, y=thetaa.bic, group=factor(simID), color=factor(simID))) +
    theme_bw() +
    geom_ribbon(aes(ymin = adjusted.bic.LB, ymax = adjusted.bic.UB, fill=factor(simID)), alpha=0.1) +
    geom_line()

master2 %>%
    arrange(thetaa.bic) %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.bic, group=factor(simID), color=factor(simID))) +
    theme_bw() +
    geom_ribbon(aes(ymin = adjusted.bic.LB, ymax = adjusted.bic.UB, fill=factor(simID)), alpha=0.1) +
    geom_line()

master2 %>%
    dplyr::filter(coverage.bic==1) %>%
    arrange(thetaa.bic) %>%
    droplevels() %>%
    ggplot(aes(x=locidF, y=thetaa.bic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = adjusted.bic.LB, ymax = adjusted.bic.UB,color=factor(simID)), alpha=0.9, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("nsim=10, estimates & bounds where RR=2 within bounds (adequate coverage)")

master2 %>%
    dplyr::filter(clusterix==1) %>%
    arrange(thetaa.bic) %>%
    droplevels() %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.bic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = adjusted.bic.log.LB, ymax = adjusted.bic.log.UB,color=factor(simID)), alpha=0.5, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("estimates & bounds inside RR=2 cluster (Buckland, BIC)")
master2 %>%
    dplyr::filter(clusterix==1) %>%
    arrange(thetaa.aic) %>%
    droplevels() %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.aic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = adjusted.aic.log.LB, ymax = adjusted.aic.log.UB,color=factor(simID)), alpha=0.5, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("estimates & bounds inside RR=2 cluster (Buckland, AIC)")


master2 %>%
    dplyr::filter(clusterix==1) %>%
    arrange(thetaa.bic) %>%
    droplevels() %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.bic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = mata.bic.LB, ymax = mata.bic.UB,color=factor(simID)), alpha=0.5, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("estimates & bounds inside RR=2 cluster (MATA, BIC)")
master2 %>%
    dplyr::filter(clusterix==1) %>%
    arrange(thetaa.aic) %>%
    droplevels() %>%
    ggplot(aes(x=fct_inorder(locidF), y=thetaa.aic))+
    geom_point() +
    theme_bw() +
    geom_pointrange(aes(ymin = mata.aic.LB, ymax = mata.aic.UB,color=factor(simID)), alpha=0.5, size=1) +
    geom_hline(yintercept = 2, color="red") +
    ggtitle("estimates & bounds inside RR=2 cluster (MATA, AIC)")

##########################################
#coverage probability
##########################################
#make a quick and dirty plot
dat <- cbind.data.frame(thetaa=as.vector(modelthetas.bic$thetaa), UBa = as.vector(bounds.bic$ma_adjusted[[2]]), 
                        LBa=as.vector(bounds.bic$ma_adjusted[[1]]))

dat <- cbind.data.frame(UBa= UBa, LBa=LBa, thetaa=as.vector(thetaa))

dat <- dat %>%arrange(thetaa) %>%mutate(id = 1:nrow(.))
ggplot(data=dat) +
    geom_line(aes(x=id, y=thetaa), size=2) +
    theme_bw() +
    geom_line(aes(x=id, y=LBa), col="blue", size=1.5, alpha=0.5) +
    geom_line(aes(x=id, y=UBa), col="blue", size=1.5, alpha=0.5) +
    geom_hline(yintercept = 1, col="red", linetype="dashed", size=1) +
    geom_vline(xintercept = 823, col="red", linetype="dashed", size=1)


# #https://github.com/jlaake/RMark/blob/master/RMark/R/mata.wald.r
# thetai <- log(modelthetas.bic$thetai)
# thetaa <- log(modelthetas.bic$thetaa)
# se.thetai <- sqrt(apply(thetai, 2, function(x) 1/outEx@x[clusterlocs.bic$param_ix]))
# w_q <- clusterlocs.bic$w_q
# thetaaa <- log(as.vector(thetaa)[827])
# thetaii <- log(thetai[,827])
# se.thetaii <- se.thetai[,827]
# 
# ########################################################
# #MATA-Intervals
# ########################################################
# thetai <- as.matrix(thetai)
# diff1 = apply(thetai,1, function(x) (as.vector(thetaa)-x))
# zquant <- sapply(1:ncol(diff1), function(x) diff1[,x]/se.thetai[x,])
# test1 <- apply(zquant, 2, function(x) pnorm(zquant))

######################################################
mataL <- uniroot(f=mata_tailareazscore, interval=c(-3, 3),
                 thetaii= thetaii,
                 se.thetaii=se.thetaii,
                 w_q=w_q, alpha=0.025, tol=1e-8)$root

mataU <- uniroot(f=mata_tailareazscore, interval=c(-3, 3),
                 thetaii= thetaii,
                 se.thetaii=se.thetaii,
                 w_q=w_q, alpha=1-0.025, tol=1e-8)$root
exp(cbind(mataL, mataU))



mata_solvebounds(thetaii = thetaii, se.thetaii = se.thetaii, w_q, alpha = 0.025, tol=1e-8)
mata_solvebounds(thetaii = thetaii, se.thetaii = se.thetaii, w_q, alpha = 0.025, tol=1e-8)
#do this for all ST locations

mata_st <- sapply(1:1040, function(x) mata_solvebounds(thetai[,x],
                                                       se.thetaii = se.thetai[,x],
                                                       w_q = w_q,
                                                       alpha = 0.025,
                                                       tol=1e-8))
mata_test <- sapply(1:10, function(x) mata_solvebounds(thetai[,x],
                                                       se.thetaii = se.thetai[,x],
                                                       w_q = w_q,
                                                       alpha = 0.025,
                                                       tol=1e-8))
mata_stt <- t(mata_st)
mata_df <- as.data.frame(mata_stt)
names(mata_df) <- c("mataLB", "mataUB")
mata_df$locid <- 1:1040

mata_df <- mata_df %>%
    mutate(clusterix = ifelse(locid %in% clusterix,1,0)) 
mata_df <- mata_df %>%
    mutate(coverage.bic = ifelse((risk < mataUB & risk > mataLB & clusterix==1), 1, 0))





