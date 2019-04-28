# Visualizing Relationship of Coal Quality to Surrogates

# Source Doc ----

setwd("G:/My Drive/R_Scripts/")
# setwd("C:/Users/jtuttle/Documents/FOPAM_CoalClassification")
source("PacifiCorp_Functions.R")


# *** ----




# Read In Chrono Data ----

rawData <- read_csv('G:/My Drive/KodyPowell_ResearchWork/PacifiCorp_Project/PlantData/Chronologic_All_5min.csv',col_types=cols(.default="d",CalendarDate="c"))

# rawData <- read_csv('Chronologic_All_5min.csv',col_types=cols(.default="d",CalendarDate="c"))

cleanedUp <- clean5minData(data=rawData)

coalDataInds <- !is.na(cleanedUp$AshSoft_Temp)
lagcoalDataInds <- !is.na(cleanedUp$AshSoft_Temp_Lagged)

stght_lag <- coalDataInds + lagcoalDataInds

keepInds <- stght_lag >= 1

wCoalData <- cleanedUp[keepInds,]

wCoalData <- wCoalData %>% filter(!is.na(Mill1_Amps))

# *** ----



# Create "Load-to-Coal" Surrogate ----


totCoalFlow <- rowSums(wCoalData[,466:470])

load2Coal <- wCoalData$MW_Gross/totCoalFlow

wCoalData <- cbind(wCoalData,load2Coal)
wCoalData <- cbind(wCoalData,totCoalFlow)

# ggplot() +
#   geom_line(aes(x=load2Coal,y=wCoalData$BTU_lb_AsRec_Lagged)) +
#   labs(x="Surrogate",y="Actual BTU")
# 
# ggplot() +
#   geom_line(aes(x=load2Coal,y=wCoalData$Ash_Per_Lagged)) +
#   labs(x="Surrogate",y="Actual Ash Percent")


# *** ----




# Create "drying air" surrogate ----

totHotAir   <- rowSums(wCoalData[,611:615])

# moisture_surr <- totCoalFlow/totHotAir

moisture_surr <- totHotAir/wCoalData$load2Coal

wCoalData <- cbind(wCoalData,moisture_surr)

# ggplot() +
#   geom_line(aes(x=moisture_surr,y=wCoalData$Moisture_AsRec_Lagged,color="Moisture_AsRec")) +
#   geom_line(aes(x=moisture_surr,y=wCoalData$H2O_Lagged,color="H2O")) +
#   labs(x="Surrogate",y="Actual Moisture")


# *** ----




# Create "hardness" surrogate ----

hardness1 <- wCoalData$Mill1_Amps/wCoalData$Mill_Flow_1_KLBH
hardness2 <- wCoalData$Mill2_Amps/wCoalData$Mill_Flow_2_KLBH
hardness3 <- wCoalData$Mill3_Amps/wCoalData$Mill_Flow_3_KLBH
hardness4 <- wCoalData$Mill4_Amps/wCoalData$Mill_Flow_4_KLBH
hardness5 <- wCoalData$Mill5_Amps/wCoalData$Mill_Flow_5_KLBH

hardness_mat <- cbind(hardness1,hardness2,hardness3,hardness4,hardness5)

hardness_surr <- vector(mode="numeric",length=nrow(wCoalData))

for(i in 1:length(hardness_surr)){
  avg_sum <- 0
  avg_cnt <- 0
  
  for(j in 466:470){
    if(wCoalData[i,j] > 20){
      avg_sum <- avg_sum + hardness_mat[i,j-465]
      avg_cnt <- avg_cnt + 1
    }
  }
  
  if(avg_cnt > 0){
    hardness_surr[i] <- avg_sum/avg_cnt  
  }else{
    hardness_surr[i] = 0
  }
  
  
  
}

wCoalData <- cbind(wCoalData,hardness_surr)

# ggplot() +
#   geom_line(aes(x=moisture_surr,y=wCoalData$Moisture_AsRec_Lagged,color="Moisture_AsRec")) +
#   geom_line(aes(x=moisture_surr,y=wCoalData$H2O_Lagged,color="H2O")) +
#   labs(x="Surrogate",y="Actual Moisture")

# *** ----




# Create "Sulfur" surrogate ----

sulfur_surr <- wCoalData$`1minLag_InSO2_PPM`/wCoalData$totCoalFlow


wCoalData <- cbind(wCoalData,sulfur_surr)

# plot_sulf <- wCoalData$`1minLag_InSO2_PPM` > 1 & wCoalData$`1minLag_InSO2_PPM` < 400
#
# ggplot() +
#   geom_line(aes(x=wCoalData$MW_Gross[plot_sulf],y=wCoalData$`1minLag_InSO2_PPM`[plot_sulf],color="InSO2")) +
#   geom_line(aes(x=wCoalData$MW_Gross,y=sulfur_surr,color="Surrogate")) +
#   labs(x="MW",y="Value")
# 
#
# 
# 
# 
# 
# # *** ----
# 
# 
# 
#
# # Density of each Coal Analyzer Value ----
# 
# ashDens <- ggplot() +
#   geom_density(aes(x=wCoalData$Ash_Per)) +
#   coord_cartesian(xlim=c(1,24)) +
#   scale_x_continuous(breaks=seq(from=1,to=24,by=2)) +
#   labs(title="Ash %")
# 
# BTUdens <- ggplot() +
#   geom_density(aes(x=wCoalData$BTU_lb_AsRec)) +
#   coord_cartesian(xlim=c(9000,14000)) +
#   scale_x_continuous(breaks=seq(from=9000,to=14000,by=500)) +
#   labs(title="BTU Content")
# 
# MoistDens <- ggplot() +
#   geom_density(aes(x=wCoalData$Moisture_AsRec)) +
#   coord_cartesian(xlim=c(5,12)) +
#   scale_x_continuous(breaks=seq(from=5,to=12,by=1)) +
#   labs(title="Moisture")
# 
# SulfDens <- ggplot() +
#   geom_density(aes(x=wCoalData$S_AsRec)) +
#   coord_cartesian(xlim=c(0,1)) +
#   scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) +
#   labs(title="Sulfur Content")
# 
# # ggplot() +
# #   geom_density(aes(x=wCoalData$`1minLag_InSO2_PPM`)) +
# #   # coord_cartesian(xlim=c(0,1)) +
# #   # scale_x_continuous(breaks=seq(from=0,to=1,by=0.1)) +
# #   labs(title="Sulfur Content")
# 
# gridExtra::grid.arrange(ashDens,BTUdens,MoistDens,SulfDens,ncol=1)
# 
# 
# 
# 
# 
# # *** ----
# 
# 
# 
# 
# # Weekly Plots of Coal Quality and Surrogates ----
# 
# starttime <- wCoalData$CalendarDate[1]
# endtime   <- wCoalData$CalendarDate[nrow(wCoalData)]
# 
# # Identify each week divding date
# curDate <- as.Date(starttime,format="%m/%d/%Y %H:%M")
# FinalDate <- as.Date(endtime,format="%m/%d/%Y %H:%M")
# 
# weekBrks <- vector(mode="character")
# weekBrks[1] <- format(curDate,format="%b/%d/%Y:%H:%M:%S")
# ind <- 2
# 
# curDate = curDate + 7
# 
# while(curDate < FinalDate){
#   
#   
#   weekBrks[ind] = format(curDate,format="%b/%d/%Y:%H:%M:%S")
#   
#   curDate = curDate + 7
#   ind = ind + 1
#   
# }
# 
# weekBrks[ind] = format(FinalDate,format="%b/%d/%Y:%H:%M:%S")
# 
# list_weekData <- vector(mode="list",length=length(weekBrks)-1)
# 
# for(i in 1:length(list_weekData)){
#   list_weekData[[i]] = timeRangeFind(weekBrks[i],weekBrks[i+1],wCoalData)
# }
# 
# # *** ----
# 
# 
# 
# 
# # Manually plot each week ----
# 
# i = 15
# 
# currWeek <- list_weekData[[i]]
#   
# ggplot() +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$MW_Gross,color="MW_Gross")) +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$totCoalFlow,color="Coal Flow")) +
#   geom_line(aes(x=currWeek$Date_Time,y=(currWeek$BTU_lb_AsRec_Lagged/100)+100,color="BTU Lagged")) +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$Ash_Per_Lagged*10,color="Ash Lagged")) +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$Moisture_AsRec_Lagged*10,color="Moisture Lagged")) +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$S_AsRec_Lagged*100,color="Sulfur Lagged")) +
#   geom_line(aes(x=currWeek$Date_Time,y=(currWeek$load2Coal*100)+70,color="BTU Surrogate")) +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$moisture_surr+50,color="Moisture Surrogate")) +
#   labs(x="Time",y="Values")
#   
# 
# 
# ggplot() +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$BTU_lb_AsRec_Lagged,color="BTU Lagged")) +
#   labs(x="Time",y="Values")
# 
# ggplot() +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$Ash_Per_Lagged,color="Ash Lagged")) +
#   labs(x="Time",y="Values")
# 
# ggplot() +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$Moisture_AsRec_Lagged,color="Moisture Lagged")) +
#   labs(x="Time",y="Values")
# 
# ggplot() +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$S_AsRec_Lagged,color="Sulfur Lagged")) +
#   labs(x="Time",y="Values")
# 
# ggplot() +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$`1minLag_InSO2_PPM`,color="Inlet SO2")) +
#   labs(x="Time",y="Values")
# 
# 
# *** ----




# Labelling Point System ----

# As explained in my notebook on April 8th (altered April 9), a point system is going to be used to label my data to start. I'll start here, get some results for the FOPAM abstract, then dial this in afterward.

# The labels will be good = 0, ok = 1, and poor = 2.

laggedData <- wCoalData %>% filter(!is.na(AshSoft_Temp_Lagged))

coalLabel <- vector(mode="numeric",length=nrow(laggedData))

for(i in 1:length(coalLabel)){
  
  curr_points = 0
  
  # Ash_Lagged Points
  if(laggedData$Ash_Per_Lagged[i] < 10.5){ }      # Good Qual
  else if(laggedData$Ash_Per_Lagged[i] < 15){   # Ok Qual
    curr_points = curr_points + 3
  }else{                                       # Bad Qual
    curr_points = curr_points + 6
  }
  
  # BTU_Lagged Points
  if(laggedData$BTU_lb_AsRec_Lagged[i] > 11750){ }      # Good Qual
  else if(laggedData$BTU_lb_AsRec_Lagged[i] > 11100){   # Ok Qual
    curr_points = curr_points + 1
  }else{                                               # Bad Qual
    curr_points = curr_points + 3
  }
  
  # Moisture_Lagged Points
  if(laggedData$Moisture_AsRec_Lagged[i] < 8.25){ }      # Good Qual
  else if(laggedData$Moisture_AsRec_Lagged[i] < 9.75){   # Ok Qual
    curr_points = curr_points + 2
  }else{                                                # Bad Qual
    curr_points = curr_points + 5
  }
  
  # Moisture_Lagged Points
  if(laggedData$S_AsRec_Lagged[i] < 0.325){ }      # Good Qual
  else if(laggedData$S_AsRec_Lagged[i] < 0.6){     # Ok Qual
    curr_points = curr_points + 1
  }else{                                          # Bad Qual
    curr_points = curr_points + 2
  }
  
  
  # Assign Label
  
  if(curr_points < 2.5){         # Good Label = 0
    coalLabel[i] = 0
  }else if(curr_points < 7.5){   # Ok Label = 1
    coalLabel[i] = 1
  }else(                         # Bad Label = 2
    coalLabel[i] = 2
  )

}

laggedData <- cbind(laggedData,coalLabel)


bad_cnt  <- sum(laggedData$coalLabel==2)

ok_cnt   <- sum(laggedData$coalLabel==1)

good_cnt <- sum(laggedData$coalLabel==0)


# ggplot() +
#   geom_line(aes(x=laggedData$Date_Time,y=laggedData$coalLabel))



# *** ----




# # Weekly Plots of Coal Quality and Labels ----
# 
# starttime <- laggedData$CalendarDate[1]
# endtime   <- laggedData$CalendarDate[nrow(laggedData)]
# 
# # Identify each week divding date
# curDate <- as.Date(starttime,format="%m/%d/%Y %H:%M")
# FinalDate <- as.Date(endtime,format="%m/%d/%Y %H:%M")
# 
# weekBrks <- vector(mode="character")
# weekBrks[1] <- format(curDate,format="%b/%d/%Y:%H:%M:%S")
# ind <- 2
# 
# curDate = curDate + 7
# 
# while(curDate < FinalDate){
# 
# 
#   weekBrks[ind] = format(curDate,format="%b/%d/%Y:%H:%M:%S")
# 
#   curDate = curDate + 7
#   ind = ind + 1
# 
# }
# 
# weekBrks[ind] = format(FinalDate,format="%b/%d/%Y:%H:%M:%S")
# 
# list_weekData <- vector(mode="list",length=length(weekBrks)-1)
# 
# for(i in 1:length(list_weekData)){
#   list_weekData[[i]] = timeRangeFind(weekBrks[i],weekBrks[i+1],laggedData)
# }
# 
# # *** ----
# 
# 
# 
# 
# # Manually plot each week ----
# 
# i = 8
# 
# currWeek <- list_weekData[[i]]
# 
# ggplot() +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$MW_Gross,color="MW_Gross")) +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$totCoalFlow,color="Coal Flow")) +
#   geom_line(aes(x=currWeek$Date_Time,y=(currWeek$BTU_lb_AsRec_Lagged/100)+100,color="BTU Lagged")) +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$Ash_Per_Lagged*10,color="Ash Lagged")) +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$Moisture_AsRec_Lagged*10,color="Moisture Lagged")) +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$S_AsRec_Lagged*100,color="Sulfur Lagged")) +
#   geom_line(aes(x=currWeek$Date_Time,y=currWeek$coalLabel*100,color="Coal Label")) +
#   # geom_line(aes(x=currWeek$Date_Time,y=(currWeek$load2Coal*100)+70,color="BTU Surrogate")) +
#   # geom_line(aes(x=currWeek$Date_Time,y=currWeek$moisture_surr+50,color="Moisture Surrogate")) +
#   labs(x="Time",y="Values")
# 
# 
# # 
# # ggplot() +
# #   geom_line(aes(x=currWeek$Date_Time,y=currWeek$BTU_lb_AsRec_Lagged,color="BTU Lagged")) +
# #   labs(x="Time",y="Values")
# # 
# # ggplot() +
# #   geom_line(aes(x=currWeek$Date_Time,y=currWeek$Ash_Per_Lagged,color="Ash Lagged")) +
# #   labs(x="Time",y="Values")
# # 
# # ggplot() +
# #   geom_line(aes(x=currWeek$Date_Time,y=currWeek$Moisture_AsRec_Lagged,color="Moisture Lagged")) +
# #   labs(x="Time",y="Values")
# # 
# # ggplot() +
# #   geom_line(aes(x=currWeek$Date_Time,y=currWeek$S_AsRec_Lagged,color="Sulfur Lagged")) +
# #   labs(x="Time",y="Values")
# # 
# # ggplot() +
# #   geom_line(aes(x=currWeek$Date_Time,y=currWeek$`1minLag_InSO2_PPM`,color="Inlet SO2")) +
# #   labs(x="Time",y="Values")
# 
# 
# 
# 
# 
# 
# 

# *** ----




# # Create 2 binary labeled sets ----
# 
# set.seed(7)
# 
# masterCoalSet <- na.omit(select(laggedData,load2Coal,moisture_surr,hardness_surr,sulfur_surr,coalLabel))
# 
# TrainSet_inds <- sample(1:nrow(masterCoalSet),floor(0.65*nrow(masterCoalSet)),replace=F)
# 
# coal_SVM_trnSet <- masterCoalSet[TrainSet_inds,]
# coal_SVM_tstSet <- masterCoalSet[-TrainSet_inds,]
# 
# 
# # Separate training set into poor/others and ok/good
# # Reset labels to {-1,1} for each set (-1 = bad, 1 = other) and (-1 = ok, 1 = good)
# 
# bad_oths_coal_SVM_trnSet <- coal_SVM_trnSet
# 
# badLabels <- bad_oths_coal_SVM_trnSet$coalLabel == 2
# bad_oths_coal_SVM_trnSet$coalLabel[badLabels] <- -1
# bad_oths_coal_SVM_trnSet$coalLabel[!badLabels] <- 1
# 
# 
# ok_good_coal_SVM_trnSet <- coal_SVM_trnSet %>% filter(coalLabel < 1.5)
# 
# okLabels <- ok_good_coal_SVM_trnSet$coalLabel == 1
# ok_good_coal_SVM_trnSet$coalLabel[okLabels] <- -1
# ok_good_coal_SVM_trnSet$coalLabel[!okLabels] <- 1
# 
# 
# 
# 
# # *** ----
# 
# 
# 
# 
# # Using package "e1071" to build binary SVMs ----
# 
# coalLabelInd <- ncol(masterCoalSet)
# 
# 
# # LINEAR SVM
# 
# svm_bad_linear <- svm(coalLabel ~ ., data = bad_oths_coal_SVM_trnSet,type="C-classification",kernel="linear")
# svm_ok_linear  <- svm(coalLabel ~ ., data = ok_good_coal_SVM_trnSet,type="C-classification",kernel="linear")
# 
# 
# svm_bad_linear_trnPreds <- predict(svm_bad_linear,bad_oths_coal_SVM_trnSet[,-coalLabelInd])
# svm_ok_linear_trnPreds  <- predict(svm_ok_linear,ok_good_coal_SVM_trnSet[,-coalLabelInd])
# 
# 
# acc_table_bad_linear    <- table(pred=svm_bad_linear_trnPreds,true=bad_oths_coal_SVM_trnSet[,coalLabelInd])
# acc_table_ok_linear     <- table(pred=svm_ok_linear_trnPreds,true=ok_good_coal_SVM_trnSet[,coalLabelInd])
# 
# 
# err_svm_bad_linear      <- sum(acc_table_bad_linear[1,2],acc_table_bad_linear[2,1])/sum(acc_table_bad_linear)
# err_svm_ok_linear       <- sum(acc_table_ok_linear[1,2],acc_table_ok_linear[2,1])/sum(acc_table_ok_linear)
# 
# 
# 
# # tuned_bad_linear        <- tune(svm,coalLabel~.,data=bad_oths_coal_SVM_trnSet,ranges=list(cost=c(0.01,0.1,1,5,10,50,100)))
# 
# tune_test <- tune(svm,coalLabel~.,data=bad_oths_coal_SVM_trnSet,ranges=list(kernel=c("linear","polynomial"),cost=c(1,5)),tunecontrol = tune.control(sampling="cross",cross=2))
# 
# 
# # summary(tune_test)
# # 
# #
# # Parameter tuning of 'svm':
# #   
# #   - sampling method: 2-fold cross validation 
# # 
# # - best parameters:
# #   kernel cost
# # polynomial    5
# # 
# # - best performance: 0.5990767 
# # 
# # - Detailed performance results:
# #   kernel cost     error  dispersion
# # 1     linear    1 0.6001212 0.004440948
# # 2 polynomial    1 0.6001543 0.004604509
# # 3     linear    5 0.6001316 0.004427734
# # 4 polynomial    5 0.5990767 0.004042746



# *** ----




# Create train and test datasets ----

set.seed(1234)

masterCoalSet <- na.omit(select(laggedData,load2Coal,moisture_surr,hardness_surr,sulfur_surr,coalLabel))

coalLabelInd <- ncol(masterCoalSet)

TrainSet_inds <- sample(1:nrow(masterCoalSet),floor(0.55*nrow(masterCoalSet)),replace=F)

coal_SVM_trnSet <- masterCoalSet[TrainSet_inds,]
coal_SVM_tstSet <- masterCoalSet[-TrainSet_inds,]

coal_multiSVM_trnSet <- coal_SVM_trnSet
coal_multiSVM_tstSet <- coal_SVM_tstSet


coal_multiSVM_trnSet[,coalLabelInd] <- as.factor(coal_multiSVM_trnSet[,coalLabelInd])
avgCol1 <- mean(coal_multiSVM_trnSet[,1])
stdCol1 <-   sd(coal_multiSVM_trnSet[,1])
avgCol2 <- mean(coal_multiSVM_trnSet[,2])
stdCol2 <-   sd(coal_multiSVM_trnSet[,2])
avgCol3 <- mean(coal_multiSVM_trnSet[,3])
stdCol3 <-   sd(coal_multiSVM_trnSet[,3])
avgCol4 <- mean(coal_multiSVM_trnSet[,4])
stdCol4 <-   sd(coal_multiSVM_trnSet[,4])

coal_multiSVM_trnSet[,1] <- (coal_multiSVM_trnSet[,1] - avgCol1)/stdCol1
coal_multiSVM_trnSet[,2] <- (coal_multiSVM_trnSet[,2] - avgCol2)/stdCol2
coal_multiSVM_trnSet[,3] <- (coal_multiSVM_trnSet[,3] - avgCol3)/stdCol3
coal_multiSVM_trnSet[,4] <- (coal_multiSVM_trnSet[,4] - avgCol4)/stdCol4


coal_multiSVM_tstSet[,coalLabelInd] <- as.factor(coal_multiSVM_tstSet[,coalLabelInd])

coal_multiSVM_tstSet[,1] <- (coal_multiSVM_tstSet[,1] - avgCol1)/stdCol1
coal_multiSVM_tstSet[,2] <- (coal_multiSVM_tstSet[,2] - avgCol2)/stdCol2
coal_multiSVM_tstSet[,3] <- (coal_multiSVM_tstSet[,3] - avgCol3)/stdCol3
coal_multiSVM_tstSet[,4] <- (coal_multiSVM_tstSet[,4] - avgCol4)/stdCol4

# sum(coal_multiSVM_trnSet$coalLabel==0)
# sum(coal_multiSVM_trnSet$coalLabel==1)
# sum(coal_multiSVM_trnSet$coalLabel==2)
# 
# sum(coal_multiSVM_tstSet$coalLabel==0)
# sum(coal_multiSVM_tstSet$coalLabel==1)
# sum(coal_multiSVM_tstSet$coalLabel==2)




# *** ----




# Using package "e1071" to build multi-classification SVMs ----


# # TEST MULTI-CLASS & TUNE FUNCTION
# testmulti_SVM <- masterCoalSet[1:500,]
# 
# testmulti_SVM[,5] <- as.factor(testmulti_SVM[,5])
# 
# multiSVM <- svm(coalLabel ~ ., data = testmulti_SVM,type="C-classification",kernel="linear")
# 
# multisVM_pred <- predict(multiSVM,testmulti_SVM[,-coalLabelInd])
# 
# 
# multiSVM_tbl <- table(pred=multisVM_pred,true=testmulti_SVM[,coalLabelInd])
# 
# err_multi <- (sum(multiSVM_tbl) - sum(diag(multiSVM_tbl)))/sum(multiSVM_tbl)
# 
# tune_test_multi <- tune(svm,coalLabel~.,data=testmulti_SVM,ranges=list(kernel=c("linear","polynomial","radial","sigmoid"),cost=c(0.01,1,100),gamma=c(1,3)),tolerance=0.005,tunecontrol = tune.control(sampling="cross",cross=2))
# 
# summary(tune_test_multi)
# 
# 
# 
# # FULL ON K-FOLD CROSS VALIDATION
# 
# strtTime <- Sys.time()
# 
# SVM_Coal_CrossVal <- tune(svm,coalLabel~.,data=coal_multiSVM_trnSet,ranges=list(kernel=c("linear","polynomial","radial","sigmoid"),cost=c(0.001,0.01,1,5,10,50,100,500,1000),gamma=c(0.001,0.01,0.05,1,5,10,50,100,500),coef0=c(0,1,5,8,10),degree=c(2,3,4,5,6,7,8)),tolerance=0.01,tunecontrol = tune.control(sampling="cross",cross=5))
# 
# summary(SVM_Coal_CrossVal)
# 
# saveRDS(SVM_Coal_CrossVal,"CoalClass_Cross_Val_Results.rds")
# 
# endtime <- Sys.time()
# 
# runDuration <- endtime - strtTime
# 
# print(runDuration)
# 
# 
# 
# 
### ----




# Parallel SVM Cross-Validation ----

strtTime <- Sys.time()


kernel_rng <- c("linear","polynomial","radial","sigmoid")

cost_rng   <- c(0.001,0.01,1,5,10,50,100,500,1000)

gamma_rng  <- c(0.001,0.01,0.05,1,5,10,50,100,500)

coef0_rng  <- c(0,1,5,8,10)

degree_rng <- c(2,3,4,5,6,7,8)

allFoldInds <- sample(1:nrow(coal_multiSVM_trnSet),nrow(coal_multiSVM_trnSet),replace=F)
nums        <- floor(nrow(coal_multiSVM_trnSet)/5)

foldInds    <- list("V1"=allFoldInds[1:nums],"V2"=allFoldInds[nums+1:nums*2],"V3"=allFoldInds[2*nums+1:nums*3],"V4"=allFoldInds[3*nums+1:nums*4],"V5"=allFoldInds[4*nums+1:nums*5])


poly_kFold_mat <- expand.grid(cost=cost_rng,gamma=gamma_rng,coef0=coef0_rng,degree=degree_rng)
rad_kFold_mat  <- expand.grid(cost=cost_rng,gamma=gamma_rng)
sig_kFold_mat  <- expand.grid(cost=cost_rng,gamma=gamma_rng,coef0=coef0_rng)



cl <- makeCluster(detectCores())
registerDoParallel(cl)




#  LINEAR KERNEL



count = 0

linval_Result <- foreach(i = 1:length(cost_rng), .combine = rbind) %do% {
  
  thisStrt <- Sys.time()
  
  cost   = cost_rng[i]
  
  out_curr <- foreach(j = 1:length(foldInds),.combine=rbind,.inorder=F) %dopar% {
    
    library(e1071)
    trn_cur = coal_multiSVM_trnSet[-foldInds[[j]],]
    tst_cur = coal_multiSVM_trnSet[foldInds[[j]],]
    
    svm_cur = svm(coalLabel ~ ., data=trn_cur, type="C-classification",kernel="linear",cost=cost,tolerance=0.01)
    cur_pred <- predict(svm_cur,tst_cur[,-coalLabelInd])
    
    cur_tbl <- table(pred=cur_pred,true=tst_cur[,coalLabelInd])
    
    cur_err <- (sum(cur_tbl) - sum(diag(cur_tbl)))/sum(cur_tbl)
    
  }
  
  
  data.frame("kernel"="linear","cost"=cost,"error"=mean(out_curr))
  
  count = count + 1
  
  paste("Linear ...",count,"/",length(cost_rng)," ... ",Sys.time() - thisStrt)
  
  
}

saveRDS(linval_Result,"CoalClass_Cross_linearVal_Results.rds")




# POLYNOMIAL KERNEL




count = 0

polyval_Result <- foreach(i = 1:nrow(poly_kFold_mat), .combine = rbind) %do% {
  
  thisStrt <- Sys.time()
  
  cost   = poly_kFold_mat$cost[i]
  gamma  = poly_kFold_mat$gamma[i] 
  coef0  = poly_kFold_mat$coef0[i] 
  degree = poly_kFold_mat$degree[i] 
  
  out_curr <- foreach(j = 1:length(foldInds),.combine=rbind,.inorder=F) %dopar% {
    
    library(e1071)
    trn_cur = coal_multiSVM_trnSet[-foldInds[[j]],]
    tst_cur = coal_multiSVM_trnSet[foldInds[[j]],]
    
    svm_cur = svm(coalLabel ~ ., data=trn_cur, type="C-classification",kernel="polynomial",cost=cost,gamma=gamma,coef0=coef0,degree=degree,tolerance=0.01)
    cur_pred <- predict(svm_cur,tst_cur[,-coalLabelInd])
    
    cur_tbl <- table(pred=cur_pred,true=tst_cur[,coalLabelInd])
    
    cur_err <- (sum(cur_tbl) - sum(diag(cur_tbl)))/sum(cur_tbl)
    
  }
  
  
  data.frame("kernel"="polynomial","cost"=cost,"gamma"=gamma,"coef0"=coef0,"degree"=degree,"error"=mean(out_curr))
  
  count = count + 1
  
  paste("Polynomial ...",count,"/",nrow(poly_kFold_mat)," ... ",Sys.time() - thisStrt)
  
  
}

saveRDS(polyval_Result,"CoalClass_Cross_polyVal_Results.rds")




# RADIAL BIAS KERNEL




count = 0

radval_Result <- foreach(i = 1:nrow(rad_kFold_mat), .combine = rbind) %do% {
  
  thisStrt <- Sys.time()
  
  cost   = rad_kFold_mat$cost[i]
  gamma  = rad_kFold_mat$gamma[i] 
  
  out_curr <- foreach(j = 1:length(foldInds),.combine=rbind,.inorder=F) %dopar% {
    
    library(e1071)
    trn_cur = coal_multiSVM_trnSet[-foldInds[[j]],]
    tst_cur = coal_multiSVM_trnSet[foldInds[[j]],]
    
    svm_cur = svm(coalLabel ~ ., data=trn_cur, type="C-classification",kernel="radial",cost=cost,gamma=gamma,tolerance=0.01)
    cur_pred <- predict(svm_cur,tst_cur[,-coalLabelInd])
    
    cur_tbl <- table(pred=cur_pred,true=tst_cur[,coalLabelInd])
    
    cur_err <- (sum(cur_tbl) - sum(diag(cur_tbl)))/sum(cur_tbl)
    
  }
  
  
  data.frame("kernel"="radial","cost"=cost,"gamma"=gamma,"error"=mean(out_curr))
  
  count = count + 1
  
  paste("Radial Bias ...",count,"/",nrow(rad_kFold_mat)," ... ",Sys.time() - thisStrt)
  
}

saveRDS(radval_Result,"CoalClass_Cross_radVal_Results.rds")




# SIGMOID KERNEL




count = 0

sigval_Result <- foreach(i = 1:nrow(sig_kFold_mat), .combine = rbind) %do% {
  
  thisStrt <- Sys.time()
  
  cost   = sig_kFold_mat$cost[i]
  gamma  = sig_kFold_mat$gamma[i] 
  coef0  = sig_kFold_mat$coef0[i] 
  
  out_curr <- foreach(j = 1:length(foldInds),.combine=rbind,.inorder=F) %dopar% {
    
    library(e1071)
    trn_cur = coal_multiSVM_trnSet[-foldInds[[j]],]
    tst_cur = coal_multiSVM_trnSet[foldInds[[j]],]
    
    svm_cur = svm(coalLabel ~ ., data=trn_cur, type="C-classification",kernel="sigmoid",cost=cost,gamma=gamma,coef0=coef0,tolerance=0.01)
    cur_pred <- predict(svm_cur,tst_cur[,-coalLabelInd])
    
    cur_tbl <- table(pred=cur_pred,true=tst_cur[,coalLabelInd])
    
    cur_err <- (sum(cur_tbl) - sum(diag(cur_tbl)))/sum(cur_tbl)
    
  }
  
  
  data.frame("kernel"="sigmoid","cost"=cost,"gamma"=gamma,"coef0"=coef0,"error"=mean(out_curr))
  
  count = count + 1
  
  paste("Sigmoid ...",count,"/",nrow(sig_kFold_mat)," ... ",Sys.time() - thisStrt)
  
}

saveRDS(sigval_Result,"CoalClass_Cross_sigVal_Results.rds")





stopCluster(cl)

endtime <- Sys.time()

runDuration <- endtime - strtTime


