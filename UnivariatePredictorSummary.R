##############################################################
# Bronwyn Rayfield, ApexRMS                                  #
# R version 3.5.1, plyr version 1.8.4, ggplot2 version 3.0.0 #
# gridExtra version 2.3, agricolae version 1.2.8             #
#                                                            #
# Script to produce univariate summaries of predictors.      #
#                                                            #
# Scripts from Rayfield et al. (in press)                    #
# Influence of habitat availability and fire disturbance on  #
# the northern range boundary of eastern white cedar (Thuja  #
# occidentalis L.). Journal of Biogeography.                 #
#                                                            #
# Table 1 and Figure S1, S2                                  #
#                                                            #
# Data available from Dryad:                                 # 
# https://doi.org/10.5061/dryad.sbcc2fr4j                    #
##############################################################

library(plyr)
library(ggplot2)
library(gridExtra)
library(agricolae)

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

attach(alldata)

col3<-rgb(255,255,255,maxColorValue=255)
col4<-rgb(225,225,225,maxColorValue=255)
col5<-rgb(178,178,178,maxColorValue=255)
col6<-rgb(130,130,130,maxColorValue=255)

##########################
# Categorical predictors #
##########################

# Function to compute Mode/Mean of predictors
OrdinalMode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

MeanMode<-ddply(alldata, "DOM_BIO", summarize, 
                CL_OM = OrdinalMode(CL_OM),
                CL_DRAI = OrdinalMode(CL_DRAI),
                HUMUS2 = OrdinalMode(HUMUS2),
                EXPO2 = OrdinalMode(EXPO2),
                SIT_PENTE = OrdinalMode(SIT_PENTE),
                CL_PENT = OrdinalMode(CL_PENT),
                CL_textB = OrdinalMode(CL_textB),
                CL_dep_sur = OrdinalMode(CL_dep_sur),
                dd5_Proj = mean(dd5_Proj),
                mmin_tenth = mean(mmin_tenth),
                gsp_Proj = mean(gsp_Proj)
)

alldata.CL_OM<-ddply(alldata, c("DOM_BIO","CL_OM"), summarise, N=length(CL_OM))
alldata.CL_OM<-ddply(alldata.CL_OM, c("DOM_BIO"), mutate, Pct=N/sum(N)*100)
alldata.CL_DRAI<-ddply(alldata, c("DOM_BIO","CL_DRAI"), summarise, N=length(CL_DRAI))
alldata.CL_DRAI<-ddply(alldata.CL_DRAI, c("DOM_BIO"), mutate, Pct=N/sum(N)*100)
alldata.HUMUS2<-ddply(alldata, c("DOM_BIO","HUMUS2"), summarise, N=length(HUMUS2))
alldata.HUMUS2<-ddply(alldata.HUMUS2, c("DOM_BIO"), mutate, Pct=N/sum(N)*100)
alldata.EXPO2<-ddply(alldata, c("DOM_BIO","EXPO2"), summarise, N=length(EXPO2))
alldata.EXPO2<-ddply(alldata.EXPO2, c("DOM_BIO"), mutate, Pct=N/sum(N)*100)
alldata.SIT_PENTE<-ddply(alldata, c("DOM_BIO","SIT_PENTE"), summarise, N=length(SIT_PENTE))
alldata.SIT_PENTE<-ddply(alldata.SIT_PENTE, c("DOM_BIO"), mutate, Pct=N/sum(N)*100)
alldata.CL_PENT<-ddply(alldata, c("DOM_BIO","CL_PENT"), summarise, N=length(CL_PENT))
alldata.CL_PENT<-ddply(alldata.CL_PENT, c("DOM_BIO"), mutate, Pct=N/sum(N)*100)
alldata.CL_textB<-ddply(alldata, c("DOM_BIO","CL_textB"), summarise, N=length(CL_textB))
alldata.CL_textB<-ddply(alldata.CL_textB, c("DOM_BIO"), mutate, Pct=N/sum(N)*100)
alldata.CL_dep_sur<-ddply(alldata, c("DOM_BIO","CL_dep_sur"), summarise, N=length(CL_dep_sur))
alldata.CL_dep_sur<-ddply(alldata.CL_dep_sur, c("DOM_BIO"), mutate, Pct=N/sum(N)*100)

plot1<-ggplot(alldata.CL_OM, aes(x=CL_OM, y=Pct, fill=factor(DOM_BIO))) + geom_bar(stat="identity", colour="black", position="dodge") + scale_x_discrete(name="Organic layer thickness (cm)") + scale_y_continuous(name="Percent", limits=c(0,75)) + scale_fill_manual(values=c(col3,col4,col5,col6), guide=FALSE) + theme_bw() + theme(axis.text=element_text(size=8), axis.title=element_text(size=10), plot.margin = unit(c(1,0,0,0), "lines"))
plot2<-ggplot(alldata.CL_DRAI, aes(x=CL_DRAI, y=Pct, fill=factor(DOM_BIO))) + geom_bar(stat="identity", colour="black", position="dodge") + scale_x_discrete(name="Drainage") + scale_y_continuous(name="", limits=c(0,75)) + scale_fill_manual(values=c(col3,col4,col5,col6), guide=FALSE) + theme_bw()  + theme(axis.text=element_text(size=8), axis.title=element_text(size=10), plot.margin = unit(c(1,0,0,0), "lines"))
plot3<-ggplot(alldata.HUMUS2, aes(x=HUMUS2, y=Pct, fill=factor(DOM_BIO))) + geom_bar(stat="identity", colour="black", position="dodge") + scale_x_discrete(name="Humus type") + scale_y_continuous(name="Percent", limits=c(0,75)) + scale_fill_manual(values=c(col3,col4,col5,col6), guide=FALSE) + theme_bw() + theme(axis.text=element_text(size=8), axis.title=element_text(size=10), plot.margin = unit(c(1,0,0,0), "lines"))
plot4<-ggplot(alldata.EXPO2, aes(x=EXPO2, y=Pct, fill=factor(DOM_BIO))) + geom_bar(stat="identity", colour="black", position="dodge") + scale_x_discrete(name="Aspect") + scale_y_continuous(name="", limits=c(0,75)) + scale_fill_manual(values=c(col3,col4,col5,col6), guide=FALSE) + theme_bw() + theme(axis.text=element_text(size=8), axis.title=element_text(size=10), plot.margin = unit(c(1,0,0,0), "lines"))
plot5<-ggplot(alldata.SIT_PENTE, aes(x=SIT_PENTE, y=Pct, fill=factor(DOM_BIO))) + geom_bar(stat="identity", colour="black", position="dodge") + scale_x_discrete(name="Slope position") + scale_y_continuous(name="Percent", limits=c(0,75)) + scale_fill_manual(values=c(col3,col4,col5,col6), guide=FALSE) + theme_bw() + theme(axis.text=element_text(size=8), axis.title=element_text(size=10), plot.margin = unit(c(1,0,0,0), "lines"))
plot6<-ggplot(alldata.CL_PENT, aes(x=CL_PENT, y=Pct, fill=factor(DOM_BIO))) + geom_bar(stat="identity", colour="black", position="dodge") + scale_x_discrete(name="Slope angle (%)") + scale_y_continuous(name="", limits=c(0,75)) + scale_fill_manual(values=c(col3,col4,col5,col6), guide=FALSE) + theme_bw() + theme(axis.text=element_text(size=8), axis.title=element_text(size=10), plot.margin = unit(c(1,0,0,0), "lines"))
plot7<-ggplot(alldata.CL_textB, aes(x=CL_textB, y=Pct, fill=factor(DOM_BIO))) + geom_bar(stat="identity", colour="black", position="dodge") + scale_x_discrete(name="Texture of B horizon") + scale_y_continuous(name="Percent", limits=c(0,75)) + scale_fill_manual(values=c(col3,col4,col5,col6), guide=FALSE) + theme_bw() + theme(axis.text=element_text(size=8), axis.title=element_text(size=10), plot.margin = unit(c(1,0,0,0), "lines"))
plot8<-ggplot(alldata.CL_dep_sur, aes(x=CL_dep_sur, y=Pct, fill=factor(DOM_BIO))) + geom_bar(stat="identity", colour="black", position="dodge") + scale_x_discrete(name="Deposit type") + scale_y_continuous(name="", limits=c(0,75)) + scale_fill_manual(values=c(col3,col4,col5,col6), guide=FALSE) + theme_bw() + theme(axis.text=element_text(size=8), axis.title=element_text(size=10), plot.margin = unit(c(1,0,0,0), "lines"))

p1<-ggplotGrob(plot1)
p2<-ggplotGrob(plot2)
p3<-ggplotGrob(plot3)
p4<-ggplotGrob(plot4)
p5<-ggplotGrob(plot5)
p6<-ggplotGrob(plot6)
p7<-ggplotGrob(plot7)
p8<-ggplotGrob(plot8)

maxWidth = grid::unit.pmax(p1$widths[2:5], p2$widths[2:5], p3$widths[2:5], p4$widths[2:5], p5$widths[2:5], p6$widths[2:5], p7$widths[2:5], p8$widths[2:5])

p1$widths[2:5] <- as.list(maxWidth)
p2$widths[2:5] <- as.list(maxWidth)
p3$widths[2:5] <- as.list(maxWidth)
p4$widths[2:5] <- as.list(maxWidth)
p5$widths[2:5] <- as.list(maxWidth)
p6$widths[2:5] <- as.list(maxWidth)
p7$widths[2:5] <- as.list(maxWidth)
p8$widths[2:5] <- as.list(maxWidth)

tiff("Categoricalpredictors.tif", width=5, height=5, units='in', res=300)
grid.arrange(p1, p2, p3, p4, p5, p6, p7,p8, ncol=2)
dev.off()


#########################
# Continuous predictors #
#########################

#Run pairwise comparisons for continuous variables
a9<-aov(dd5_Proj~DOM_BIO, data=alldata)
a10<-aov(mmin_tenth~DOM_BIO, data=alldata)
a11<-aov(gsp_Proj~DOM_BIO, data=alldata)

tuk9<-HSD.test(a9, "DOM_BIO", group=TRUE) 
tuk10<-HSD.test(a10, "DOM_BIO", group=TRUE) 
tuk11<-HSD.test(a11, "DOM_BIO", group=TRUE) 

#Plot continous variables
plot9<-ggplot(alldata, aes(x=factor(DOM_BIO), y=dd5_Proj, fill=factor(DOM_BIO))) + geom_violin(scale = "width") + geom_boxplot(width=.1, outlier.size=0) + scale_x_discrete(name="") + scale_y_continuous(name="Growing degree days") + scale_fill_manual(values=c(col3,col4,col5,col6), guide=FALSE) + theme_bw() + theme(plot.margin = unit(c(1,1,0,1), "lines"), axis.title=element_text(size=9)) # + annotate("text", x=c(1.15, 2.15, 3.15, 4.15), y=max(dd5_Proj), label=as.character(tuk9$groups$groups), size=3)
plot10<-ggplot(alldata, aes(x=factor(DOM_BIO), y=mmin_tenth, fill=factor(DOM_BIO))) + geom_violin(scale = "width") + geom_boxplot(width=.1, outlier.size=0) + scale_x_discrete(name="") + scale_y_continuous(name=expression(Minimum~~temperature~~"("~degree~C~")")) + scale_fill_manual(values=c(col3,col4,col5,col6), guide=FALSE) + theme_bw() + theme(plot.margin = unit(c(0,1,0,1), "lines"), axis.title=element_text(size=9)) # + annotate("text", x=c(1.15, 2.15, 3.15, 4.15), y=max(mmin_tenth), label=as.character(tuk10$groups$groups), size=3)
plot11<-ggplot(alldata, aes(x=factor(DOM_BIO), y=gsp_Proj, fill=factor(DOM_BIO))) + geom_violin(scale = "width") + geom_boxplot(width=.1, outlier.size=0) + scale_x_discrete(name="Bioclimatic domain")  + scale_y_continuous(name="Precipitation (mm)") + scale_fill_manual(values=c(col3,col4,col5,col6), guide=FALSE) + theme_bw() + theme(plot.margin = unit(c(0,1,0,1), "lines"), axis.title=element_text(size=9)) # + annotate("text", x=c(1.15, 2.15, 3.15, 4.15), y=max(gsp_Proj), label=as.character(tuk11$groups$groups), size=3)

tiff("Continuouspredictors.tif", width=5, height=5, units='in', res=300)
  grid.arrange(plot9, plot10, plot11, ncol=1) 
dev.off()



