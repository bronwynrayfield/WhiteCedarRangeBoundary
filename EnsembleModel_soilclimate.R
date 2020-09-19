##############################################################
# Bronwyn Rayfield, ApexRMS                                  #
# R version 3.5.1, biomod2 version 3.3-7                     #
#                                                            #
# Script to build an ensemble white cedar distribution model # 
# based on site-level habitat data (soil & topographic) and  #
# climate data.                                              #
#                                                            #
# Scripts from Rayfield et al. (in press)                    #
# Influence of habitat availability and fire disturbance on  #
# the northern range boundary of eastern white cedar (Thuja  #
# occidentalis L.). Journal of Biogeography.                 #
#                                                            #
# Table 2, 3, 4                                              #
#                                                            #
# Data available from Dryad:                                 # 
# https://doi.org/10.5061/dryad.sbcc2fr4j                    #
##############################################################

library(biomod2)

################
# 1. Prep data #
################
# Read in data, order categorical variables and normalize continuous variables
dataName <- "2020RayfieldForestSampleSiteData.csv"
alldata <- data.frame(read.table(file=dataName, header=TRUE))

# Re-order the factor levels
alldata$CL_OM=as.factor(alldata$CL_OM)
alldata$CL_DRAI=as.factor(alldata$CL_DRAI)
alldata$CL_DRAI<-factor(alldata$CL_DRAI, levels=c("Poor","Moderate","Good"))
alldata$HUMUS2=as.factor(alldata$HUMUS2)
alldata$HUMUS2<-factor(alldata$HUMUS2, levels=c("Org", "MD", "MR"))
alldata$EXPO2=as.factor(alldata$EXPO2)
alldata$EXPO2<-factor(alldata$EXPO2, levels=c("West", "North", "South", "East", "Total"))
alldata$SIT_PENTE=as.factor(alldata$SIT_PENTE)
alldata$SIT_PENTE<-factor(alldata$SIT_PENTE, levels=c("Bottom", "Middle", "Top"))
alldata$CL_PENT=as.factor(alldata$CL_PENT)
alldata$CL_dep_sur=as.factor(alldata$CL_dep_sur)
alldata$CL_dep_sur<-factor(alldata$CL_dep_sur, levels=c("Organic", "Fine", "Medium", "Coarse", "Rock"))
alldata$CL_textB=as.factor(alldata$CL_textB)
alldata$CL_textB<-factor(alldata$CL_textB, levels=c("Organic", "Fine", "Medium", "Coarse"))
alldata$TOUT_THO=as.numeric(alldata$TOUT_THO)

# Divide data into south (3 & 4) and north (5 & 6) bioclimatic domains
sol<-alldata[which(alldata$DOM_BIO==3 | alldata$DOM_BIO==4),c("TOUT_THO", "Lat", "Long", "CL_OM", "CL_DRAI", "HUMUS2", "EXPO2", "SIT_PENTE", "CL_PENT", "CL_textB", "CL_dep_sur", "dd5_Proj", "mmin_tenth", "gsp_Proj")]
sol_nord<-alldata[which(alldata$DOM_BIO==5 | alldata$DOM_BIO==6),c("TOUT_THO", "Lat", "Long", "CL_OM", "CL_DRAI", "HUMUS2", "EXPO2", "SIT_PENTE", "CL_PENT", "CL_textB", "CL_dep_sur", "dd5_Proj", "mmin_tenth", "gsp_Proj")]

# Normalize continuous variables in the south
meansol_dd5_Proj<-mean(sol$dd5_Proj)
meansol_mmin_tenth<-mean(sol$mmin_tenth)
meansol_gsp_Proj<-mean(sol$gsp_Proj)
sdsol_dd5_Proj<-sd(sol$dd5_Proj)
sdsol_mmin_tenth<-sd(sol$mmin_tenth)
sdsol_gsp_Proj<-sd(sol$gsp_Proj)
sol[,c("dd5_Proj","mmin_tenth","gsp_Proj")]<-scale(sol[,c("dd5_Proj","mmin_tenth","gsp_Proj")])
# Normalize continuous variables in the north based on the transformation done in the south
sol_nord[,"dd5_Proj"]<-(sol_nord[,"dd5_Proj"]-meansol_dd5_Proj)/sdsol_dd5_Proj
sol_nord[,"mmin_tenth"]<-(sol_nord[,"mmin_tenth"]-meansol_mmin_tenth)/sdsol_mmin_tenth
sol_nord[,"gsp_Proj"]<-(sol_nord[,"gsp_Proj"]-meansol_gsp_Proj)/sdsol_gsp_Proj


################################################
# 2. Prepare data in south ( 3 & 4) for BIOMOD #
################################################
# Withhold 30% of the data in the south for validation
# Identify rows to use for evaluating and rows for training
sol_rowsforevaluating<-sample(nrow(sol), round(nrow(sol)*.3))
sol_rowsfortraining<-c(1:nrow(sol))[-sol_rowsforevaluating]

# Presence/absences data
myResp <- as.numeric(sol[sol_rowsfortraining,"TOUT_THO"])
# XY coordinates
myRespXY <- sol[sol_rowsfortraining,c("Lat","Long")]
# Environmental variables
myExpl_soilclim <- sol[sol_rowsfortraining,c("CL_OM", "CL_DRAI", "HUMUS2", "EXPO2", "SIT_PENTE", "CL_PENT", "CL_textB", "CL_dep_sur", "dd5_Proj", "mmin_tenth", "gsp_Proj")]

evalResp <- as.numeric(sol[sol_rowsforevaluating,"TOUT_THO"])
evalRespXY <- sol[sol_rowsforevaluating,c("Lat","Long")]
evalExpl_soilclim <- sol[sol_rowsforevaluating,c("CL_OM", "CL_DRAI", "HUMUS2", "EXPO2", "SIT_PENTE", "CL_PENT", "CL_textB", "CL_dep_sur", "dd5_Proj", "mmin_tenth", "gsp_Proj")]

myBiomodData_soilclim <- BIOMOD_FormatingData(resp.var = myResp,
                                              expl.var = myExpl_soilclim,
                                              resp.xy = myRespXY,
                                              resp.name = "cedar.soilclimate",
                                              eval.resp.var = evalResp,
                                              eval.expl.var = evalExpl_soilclim,
                                              eval.resp.xy = evalRespXY)

myBiomodData_soilclim


########################
# 3. Run BIOMOD models #
########################

# Define Models Options
myBiomodOption <- BIOMOD_ModelingOptions(GAM = list(k=3))

# Compute the models
myBiomodModelOut_soilclim <- BIOMOD_Modeling(
  myBiomodData_soilclim,
  models = c('GLM','GBM','GAM','ANN','MARS','RF'), #,'MAXENT.Phillips','MAXENT.Tsuruoka'
  models.options = myBiomodOption,
  NbRunEval=100,
  DataSplit=70,
  Prevalence=NULL,
  VarImport=20,
  models.eval.meth = c('KAPPA','TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE)

save.image("EnsembleModel_soilclimate.RData")

myBiomodModelOut_soilclim

# Models evaluations
myBiomodModelEval_soilclim <- get_evaluations(myBiomodModelOut_soilclim)
# Print the dimnames of this object
dimnames(myBiomodModelEval_soilclim)

# Print evaluation scores of all models
myBiomodModelEval_soilclim[,"Testing.data",,,]
# Print ROC scores of all models
myBiomodModelEval_soilclim["ROC","Testing.data",,,]
# Summary of evaluation metrics
rowMeans(myBiomodModelEval_soilclim[,"Testing.data",,,])
apply(myBiomodModelEval_soilclim["TSS","Testing.data",,,],1, mean)
apply(myBiomodModelEval_soilclim["TSS","Testing.data",,,],1, sd)
#apply(myBiomodModelEval_soilclim["KAPPA","Testing.data",,,],1, mean)
apply(myBiomodModelEval_soilclim["ROC","Testing.data",,,],1, mean)
apply(myBiomodModelEval_soilclim["ROC","Testing.data",,,],1, sd)

apply(myBiomodModelEval_soilclim["TSS","Evaluating.data",,,],1, mean)
apply(myBiomodModelEval_soilclim["TSS","Evaluating.data",,,],1, sd)
#apply(myBiomodModelEval_soilclim["KAPPA","Evaluating.data",,,],1, mean)
apply(myBiomodModelEval_soilclim["ROC","Evaluating.data",,,],1, mean)
apply(myBiomodModelEval_soilclim["ROC","Evaluating.data",,,],1, sd)

# Variable importance
myBiomodVarImportance_soilclim<-get_variables_importance(myBiomodModelOut_soilclim)
apply(myBiomodVarImportance_soilclim[,,,],c(1,2), mean)
apply(myBiomodVarImportance_soilclim[,,,],c(1,2), sd)

#########################
# 4. Ensemble modelling #
#########################

myBiomodEM_soilclim <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut_soilclim,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('ROC'),
#  eval.metric.quality.threshold = c(0.6),
  models.eval.meth=c('TSS','ROC'),
  prob.mean = T,
  prob.cv = F,
  prob.ci = F,
  prob.median = F,
  committee.averaging = F,
  prob.mean.weight = F,
  prob.mean.weight.decay = 'proportional' )

myBiomodEM_soilclim


# Get evaluation scores
myBiomodEMEval_soilclim<-get_evaluations(myBiomodEM_soilclim)
myBiomodEMEval_soilclim

#####################################################
# 5. Forecast white cedar presence absence in north #
#####################################################
# Project individual models in the south on training data
myBiomodProj_south_train <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_soilclim,
  new.env = myExpl_soilclim,
  xy.new.env = myRespXY,
  proj.name = '',
  selected.models = 'all',
  binary.meth = 'ROC',
  compress = 'xz',
  clamping.mask = F)

myBiomodProj_south_train

myBiomodEF_south_train <- BIOMOD_EnsembleForecasting(EM.output= myBiomodEM_soilclim, projection.output=myBiomodProj_south_train)

myEMPred_south_train<-cbind(pred=get_predictions(myBiomodEF_south_train)/1000, obs=myResp)

myEMPred_south_train<-cbind(myEMPred_south_train,mse=abs(myEMPred_south_train[,1] - myEMPred_south_train[,2]))
mean(myEMPred_south_train[myEMPred_south_train[,2]==1,3])
mean(myEMPred_south_train[myEMPred_south_train[,2]==0,3])
mean(myEMPred_south_train[,3])

x11(); hist(myEMPred_south_train[myEMPred_south_train[,2]==1,1])
x11(); hist(myEMPred_south_train[myEMPred_south_train[,2]==0,1])

Cutoff<-myBiomodEMEval_soilclim[[1]][2,3]

myEMPred_south_train<-cbind(myEMPred_south_train,pred_cutoff=0)
myEMPred_south_train[,4][which(myEMPred_south_train[,1]>=Cutoff)]<-1

table(myEMPred_south_train[,4],myEMPred_south_train[,2])

myEMPred_south_train <- cbind(myEMPred_south_train, myBiomodEF_south_train@xy.coord)
myEMPred_south_train <- data.frame(myEMPred_south_train,obspred=0)
myEMPred_south_train[which(myEMPred_south_train$pred==0 & myEMPred_south_train$obs==1), "obspred"] <- 1
myEMPred_south_train[which(myEMPred_south_train$pred==1 & myEMPred_south_train$obs==0), "obspred"] <- 2
myEMPred_south_train[which(myEMPred_south_train$pred==1 & myEMPred_south_train$obs==1), "obspred"] <- 3

write.csv(myEMPred_south_train,"EnsembleModelPredictions_South_train_soilclim.csv", row.names=FALSE)


# Project ensemble model in the south on evaluation data
myBiomodProj_south_eval <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_soilclim,
  new.env = evalExpl_soilclim,
  xy.new.env = evalRespXY,
  proj.name = '',
  selected.models = 'all',
  binary.meth = 'ROC',
  compress = 'xz',
  clamping.mask = F)

myBiomodProj_south_eval

myBiomodEF_south_eval <- BIOMOD_EnsembleForecasting(EM.output= myBiomodEM_soilclim,
                                                    projection.output=myBiomodProj_south_eval)

myEMPred_south_eval<-cbind(pred=get_predictions(myBiomodEF_south_eval)/1000, obs=evalResp)
myEMPred_south_eval<-cbind(myEMPred_south_eval,mse=abs(myEMPred_south_eval[,1] - myEMPred_south_eval[,2]))
mean(myEMPred_south_eval[myEMPred_south_eval[,2]==1,3])
mean(myEMPred_south_eval[myEMPred_south_eval[,2]==0,3])
mean(myEMPred_south_eval[,3])

x11(); hist(myEMPred_south_eval[myEMPred_south_eval[,2]==1,1])
x11(); hist(myEMPred_south_eval[myEMPred_south_eval[,2]==0,1])

Cutoff<-myBiomodEMEval_soilclim[[1]][2,3]

myEMPred_south_eval<-cbind(myEMPred_south_eval,pred_cutoff=0)
myEMPred_south_eval[,4][which(myEMPred_south_eval[,1]>=Cutoff)]<-1

table_S_eval<-table(myEMPred_south_eval[,4],myEMPred_south_eval[,2])
