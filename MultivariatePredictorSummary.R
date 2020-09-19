##############################################################
# Bronwyn Rayfield, ApexRMS                                  #
# R version 3.5.1, ade version 1.7.11, vegan version 2.5.2,  #
# plyr version 1.8.4, ggplot2 version 3.0.0,                 #  
# gclus version 1.3.1, ape version 5.1, FD version 1.0.12    #
#                                                            #
# Script to produce multivariate summary of predictors.      #
#                                                            #
# Scripts from Rayfield et al. (in press)                    #
# Influence of habitat availability and fire disturbance on  #
# the northern range boundary of eastern white cedar (Thuja  #
# occidentalis L.). Journal of Biogeography.                 #
#                                                            #
# Figure 2                                                   #
#                                                            #
# Data available from Dryad:                                 # 
# https://doi.org/10.5061/dryad.sbcc2fr4j                    #
##############################################################

library(plyr)
library(ggplot2)
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(FD)

# Functions
# coldiss()
# Color plots of a dissimilarity matrix, without and with ordering
#
# License: GPL-2 
# Author: Francois Gillet, August 2009
#
"coldiss" <- function(D, nc = 4, byrank = TRUE, diag = FALSE){
  require(gclus)
  
  if (max(D)>1) D <- D/max(D)
  
  if (byrank) {
    spe.color = dmat.color(1-D, cm.colors(nc))
  }
  else {
    spe.color = dmat.color(1-D, byrank=FALSE, cm.colors(nc))
  }
  
  spe.o = order.single(1-D)
  speo.color = spe.color[spe.o,spe.o]
  
  op = par(mfrow=c(1,2), pty="s")
  
  if (diag) {
    plotcolors(spe.color, rlabels=attributes(D)$Labels, 
               main="Dissimilarity Matrix", 
               dlabels=attributes(D)$Labels)
    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
               main="Ordered Dissimilarity Matrix", 
               dlabels=attributes(D)$Labels[spe.o])
  }
  else {
    plotcolors(spe.color, rlabels=attributes(D)$Labels, 
               main="Dissimilarity Matrix")
    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
               main="Ordered Dissimilarity Matrix")
  }
  
  par(op)
}

#Read in full dataset
dataName <- "2020RayfieldForestSampleSiteData.csv"
alldata <- data.frame(read.csv(file=dataName, header=TRUE))
#edaphic and climate
alldata<-alldata[,c("DOM_BIO", "TOUT_THO", "Lat", "Long", "CL_OM", "CL_DRAI", "HUMUS2", "EXPO2", "SIT_PENTE", "CL_PENT", "CL_textB", "CL_dep_sur", "dd5_Proj", "mmin_tenth", "gsp_Proj")]

#Re-order the factor levels
alldata$CL_OM=as.factor(alldata$CL_OM)
levels(alldata$CL_OM)<-c("0-10","10-20","20-30",">=30")
alldata$CL_DRAI=as.factor(alldata$CL_DRAI)
alldata$CL_DRAI<-factor(alldata$CL_DRAI, levels=c("Poor","Moderate","Good"))
alldata$HUMUS2=as.factor(alldata$HUMUS2)
alldata$HUMUS2<-factor(alldata$HUMUS2, levels=c("Org", "MD", "MR"))
levels(alldata$HUMUS2)<-c("Organic", "Moder", "Mor")
alldata$EXPO2=as.factor(alldata$EXPO2)
alldata$EXPO2<-factor(alldata$EXPO2, levels=c("West", "North", "South", "East", "Total"))
alldata$SIT_PENTE=as.factor(alldata$SIT_PENTE)
alldata$SIT_PENTE<-factor(alldata$SIT_PENTE, levels=c("Bottom", "Middle", "Top"))
alldata$CL_PENT=as.factor(alldata$CL_PENT)
levels(alldata$CL_PENT)<-c("0-3", "4-8", "9-15", "16-30", ">=30")
alldata$CL_dep_sur=as.factor(alldata$CL_dep_sur)
alldata$CL_dep_sur<-factor(alldata$CL_dep_sur, levels=c("Organic", "Fine", "Medium", "Coarse", "Rock"))
alldata$CL_textB=as.factor(alldata$CL_textB)
alldata$CL_textB<-factor(alldata$CL_textB, levels=c("Organic", "Fine", "Medium", "Coarse"))
alldata$TOUT_THO=as.numeric(alldata$TOUT_THO)

# Keep only habitat predictors (soil and topographic predictors)
mydata<-alldata[,1:12]

mydata.d<-gowdis(mydata[,5:12], ord="podani")
coldiss(mydata.d, diag=TRUE)
 
n <- length(10)
dim <- data.frame(k=1:10, s=0)
mydata.NMDS <- metaMDS(mydata.d, k=1)
dim$s[1] <- mydata.NMDS$stress
mydata.NMDS <- metaMDS(mydata.d, k=2)
dim$s[2] <- mydata.NMDS$stress
mydata.NMDS <- metaMDS(mydata.d, k=3)
dim$s[3] <- mydata.NMDS$stress
mydata.NMDS <- metaMDS(mydata.d, k=4)
dim$s[4] <- mydata.NMDS$stress
mydata.NMDS <- metaMDS(mydata.d, k=5)
dim$s[5] <- mydata.NMDS$stress
mydata.NMDS <- metaMDS(mydata.d, k=6)
dim$s[6] <- mydata.NMDS$stress
mydata.NMDS <- metaMDS(mydata.d, k=7)
dim$s[7] <- mydata.NMDS$stress
mydata.NMDS <- metaMDS(mydata.d, k=8)
dim$s[8] <- mydata.NMDS$stress
mydata.NMDS <- metaMDS(mydata.d, k=9)
dim$s[9] <- mydata.NMDS$stress
mydata.NMDS <- metaMDS(mydata.d, k=10)
dim$s[10] <- mydata.NMDS$stress

write.table(dim,"ScreePlotInput.txt",sep=",", row.names=FALSE)
tiff("NMDS_screeplot_soil.tiff", width=5, height=5, units='in', res=300)
  plot(dim$k, dim$s, type="b", main="Stress vs dimensions", xlab="Dimensions (k)", ylab="Stress")
dev.off()


# Final NMDS result	- two dimensions
mydata.NMDS <- metaMDS(mydata.d, trace=T, plot=T)

#Join NMDS results to dataset
mydata<-data.frame(mydata, scores(mydata.NMDS))

save.image("nmds_edaphic.RData")

#plot 4 classes of domaine biologique
mydata$DOM_BIO <- factor(mydata$DOM_BIO)
levels(mydata$DOM_BIO) <- rev(levels(mydata$DOM_BIO))
mydata$TOUT_THO <- factor(mydata$TOUT_THO)

vec.soil<-envfit(mydata.NMDS~mydata$CL_DRAI + mydata$HUMUS2 + mydata$EXPO2, data=mydata, perm=999)
centroids<-vec.soil$factors$centroids
centroids<-cbind(centroids, Label=c("Poor drainage", "Moderate drainage", "Good drainage", "Organic", "Moder", "Mor", "Western aspect", "Northern aspect", "Southern aspect", "Eastern aspect", "Total aspect"))

# Plot with centroids
tiff("NMDS_2biodom_edaphic_centroids.tif", width=5, height=5, units='in', res=300)
ggplot(mydata, aes(x=MDS1, y=MDS2, color=DOM_BIO, shape=TOUT_THO)) + theme(legend.position=c(1,0.3),legend.justification=c(1,1)) + geom_point() + scale_colour_manual(values=c("black","black","grey69","grey69"), name ="Bioclimatic\ndomains") + scale_shape_manual(values = c(1, 19), name="Presence of\nwhite cedar") + theme_bw()+ annotate(geom="label", x=as.numeric(centroids[1:3,"NMDS1"]), y=as.numeric(centroids[1:3,"NMDS2"]), label=centroids[1:3, "Label"], color="black")
dev.off()

# No legend
tiff("NMDS_2biodom_edaphic_centroids_nolegend.tif", width=5, height=5, units='in', res=300)
ggplot(mydata, aes(x=MDS1, y=MDS2, color=DOM_BIO, shape=TOUT_THO)) + geom_point() + scale_colour_manual(values=c("black","black","grey69","grey69")) + scale_shape_manual(values = c(1, 19)) + theme_bw()+ annotate(geom="label", x=as.numeric(centroids[1:3,"NMDS1"]), y=as.numeric(centroids[1:3,"NMDS2"]), label=centroids[1:3, "Label"], color="black") + theme(legend.position="none")
dev.off()





