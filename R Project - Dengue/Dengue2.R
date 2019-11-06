#### Load libraries & Setup ####
library(pacman)

p_load(doParallel, here, readr, rstudioapi,        #parallel computing, relative path
       caret, C50, caretEnsemble, mboost, mlr, Metrics, randomForest, party, MASS, gbm,
       #ParamHelpers, hydroGOF #Classification and Regression
       cluster, corrplot, fpc, e1071, recipes, Hmisc, #Clustering: Corr visulization, Clustering & SVM, 
       ggplot2, ggpubr, RColorBrewer, lattice, dygraphs, #Visualization
       ade4, inum, reshape,  #Cleaning, Preprocessing
       #FactoMineR, factoextra, #PCA, MCA
       plyr, dplyr, tidyr, tidyverse, textclean, 
       #arules, arulesViz,    # ASsociation Rules Mining: analyzing and visualize transactional data
       #markdown, shiny, tinytex, rmdformats, knitr #html docu, dashboard, Latex for PDF docu
       RMySQL, lubridate, # Time Series: MySQL, functions for DateTime
       BBmisc, #asbio, # BBmisc: normalization
       imputeTS, padr, #interpolate missing values (Time Series)
       forecast, tseries, prophet #aTSA #Time Series
)


## Enable parallel computing
cl <- makePSOCKcluster(4)
registerDoParallel(cl)


## Disable scientific notation
options(scipen = 999)


## File directory
current_path = getActiveDocumentContext()$path
setwd(dirname(current_path))
getwd()
setwd("..")
setwd("dataset")


## Load data
dengue_all_raw <- read.csv(file = "../dataset/dengue_all_raw.csv")
colnames(dengue_all_raw)[1] <- "city"

# Split per city
dengue_sj <- dengue_all_raw[which(dengue_all_raw$city == "sj"), ]
dengue_iq <- dengue_all_raw[which(dengue_all_raw$city == "iq"), ]

# Set timezone
dengue_sj$week_start_date <- as.POSIXct(dengue_sj$week_start_date, format = "%d/%m/%Y")
#dengue_sj$week_start_date <- as.POSIXct(dengue_sj$week_start_date, format = "%d/%m/%Y", tz = "America/Havana") --> one NA, Why?
tz(dengue_sj$week_start_date)

dengue_iq$week_start_date <- as.POSIXct(dengue_iq$week_start_date, format = "%d/%m/%Y", tz = "America/Lima")
tz(dengue_iq$week_start_date)

## Check data
# NAs7
sum(is.na(dengue_sj))
which(is.na(dengue_sj))

sum(is.na(dengue_iq))

#### Benchmark - Prophet ####
dengue_sj_p <- dengue_sj[, c(4, 25)]
colnames(dengue_sj_p) <- c("ds", "y")

dengue_iq_p <- dengue_iq[, c(4, 25)]
colnames(dengue_iq_p) <- c("ds", "y")

# Model
Proph1_sj <- prophet(dengue_sj_p, weekly.seasonality = TRUE, daily.seasonality = FALSE)
Proph1_iq <- prophet(dengue_iq_p, weekly.seasonality = TRUE, daily.seasonality = FALSE)

# Forecast
future1_sj <- make_future_dataframe(Proph1_sj, periods = 260, freq = 'week')
forecast1_sj <- predict(Proph1_sj, future1_sj)

future1_iq <- make_future_dataframe(Proph1_iq, periods = 156, freq = 'week')
forecast1_iq <- predict(Proph1_iq, future1_iq)

# Visualize Forecast
dyplot.prophet(Proph1_sj, forecast1_sj)
dyplot.prophet(Proph1_iq, forecast1_iq)

# Evaluation
proph1_sj_CV <- cross_validation(Proph1_sj, 
                                 initial = 675, 
                                 horizon = 260, units = "weeks") 
proph1_sj_Acc <- performance_metrics(proph1_sj_CV)
plot_cross_validation_metric(proph1_sj_CV, metric = "mae")

proph1_iq_CV <- cross_validation(Proph1_iq, 
                                 initial = 364, 
                                 horizon = 156, units = "weeks") 
proph1_iq_Acc <- performance_metrics(proph1_iq_CV)
plot_cross_validation_metric(proph1_iq_CV, metric = "mae")

# Save to CSV
city <- c(rep("sj", nrow(forecast1_sj)), rep("iq", nrow(forecast1_iq)))
forecast1_sj$ds <- forecast1_sj$ds + 8201 # due to tz, shift of 2 hours back, leading to wrong week

proph1_fc_compl <- cbind(rbind(forecast1_sj, forecast1_iq), city)
proph1_fc_compl$year <- year(proph1_fc_compl$ds)
proph1_fc_compl$weekofyear <- week(proph1_fc_compl$ds)
proph1_fc_compl$total_cases <- proph1_fc_compl$yhat
proph1_fc_compl <- proph1_fc_compl[c(937:1196, 1717:1872), c(20:23)]

write.csv(x = proph1_fc_compl, "proph1_fc_compl.csv")


#### Preprocessing ####
## Iteration 1 (rmvd: ndvi_ne, ndvi_nw; station_avg_temp_c(calc from min max) and min max flld w/ forecast("=reanalysis_min_air_temp_k,"=reanalysis_max_air_temp_k)
dengue_all_iter1 <- read.csv(file = "../dataset/dengue_all_raw_NA_flld.csv")
sum(is.na(dengue_all_iter1)) #263 NAs

# Pre-select columns and check NAs
dengue_sj_iter1 <- dengue_all_iter1[which(dengue_all_iter1$city == "sj"), c(1:4,7:26)]
sum(is.na(dengue_sj_iter1)) #102
dengue_sj_iter1$sumNAs <- rowSums(is.na(dengue_sj_iter1))
sum(dengue_sj_iter1$sumNAs)
dengue_sj_iter1 <- dengue_sj_iter1[which(dengue_sj_iter1$sumNAs == 0), c(1:24)] #936 - 927 NA rows rmvd (9)

dengue_iq_iter1 <- dengue_all_iter1[which(dengue_all_iter1$city == "iq"), c(1:4,7:26)]
dengue_iq_iter1$sumNAs <- rowSums(is.na(dengue_iq_iter1))
sum(dengue_iq_iter1$sumNAs) #117
dengue_iq_iter1 <- dengue_iq_iter1[which(dengue_iq_iter1$sumNA == 0), c(1:24)] #520 - 473 NA rows rmvd (47)

sum(is.na(dengue_sj_iter1)) #0
sum(is.na(dengue_iq_iter1)) #0

# Define df for next step
dengue_sj <- dengue_sj_iter1
dengue_iq <- dengue_iq_iter1

## Check correlation
# Create correlation-matrix
dengue_sj_corr <- rcorr(as.matrix(dengue_sj[,c(5:18, 21:24)]))           
dengue_iq_corr <- rcorr(as.matrix(dengue_iq[,c(5:18, 21:24)]))  

# Visualization of correlation matrix             
corrplot(dengue_sj_corr$r, method = "number",                                       
         addrect = 6, rect.col = "black",                                       
         type = "upper",                                                       
         tl.col="black", tl.srt=90, tl.cex = 0.65, tl.pos = "t",               
         addgrid.col = "lightgray")

corrplot(dengue_iq_corr$r, method = "number",                                       
         addrect = 6, rect.col = "black",                                          
         type = "upper",
         tl.col="black", tl.srt=90, tl.cex = 0.65, tl.pos = "t",
         addgrid.col = "lightgray")

## Calculate and check lagged correlation (cross-correlation)
i <- NULL
lagCorr_sj <- NULL
lagCorr_iq <- NULL
p <- ccf(dengue_sj[4+1], dengue_sj$total_cases, lag.max = 10)
q <- ccf(dengue_iq[4+1], dengue_iq$total_cases, lag.max = 10)
lagCorr_sj <- data.frame(p$lag)
lagCorr_iq <- data.frame(q$lag)
colnames(lagCorr_sj)[1] <- "lag"
colnames(lagCorr_iq)[1] <- "lag"

for (i in c(1:ncol(dengue_sj[,c(5:24)]))){
  p <- ccf(dengue_sj[4+i], dengue_sj$total_cases, lag.max = 10, plot = FALSE)
  q <- ccf(dengue_iq[4+i], dengue_iq$total_cases, lag.max = 10, plot = FALSE)
  lagCorr_sj[1+i] <- p$acf
  lagCorr_iq[1+i] <- q$acf
  colnames(lagCorr_sj)[1+i] <- colnames(dengue_sj)[4+i]
  colnames(lagCorr_iq)[1+i] <- colnames(dengue_iq)[4+i]
  print(plot(lagCorr_sj[,c(1,i+1)]))
  print(plot(lagCorr_iq[,c(1,i+1)]))
}

# Identify and select lag of highest cross-correlation
j <- NULL
max <- NULL
lag <- NULL
varLag_sj <- NULL
for (j in c(2:length(lagCorr_sj))){
  max <- c(max, max(lagCorr_sj[,j]))
  lag <- c(lag, lagCorr_sj[which(lagCorr_sj[,j] == max[j-1]), 1])
}
varLag_sj <- data.frame(cbind(lag, max, colnames(lagCorr_sj)[c(2:21)]))
varLag_sj$lag <- as.character(varLag_sj$lag)
varLag_sj$lag <- as.numeric(varLag_sj$lag)
class(varLag_sj$lag)

jj <- NULL
max1 <- NULL
lag1 <- NULL
varLag_iq <- NULL
for (jj in c(2:length(lagCorr_iq))){
  max1 <- c(max1, max(lagCorr_iq[,jj]))
  lag1 <- c(lag1, lagCorr_iq[which(lagCorr_iq[,jj] == max1[jj-1]), 1])
}
varLag_iq <- data.frame(cbind(lag1, max1, colnames(lagCorr_iq)[c(2:21)]))
colnames(varLag_iq) <- colnames(varLag_sj)
varLag_iq$lag <- as.character(varLag_iq$lag)
varLag_iq$lag <- as.numeric(varLag_iq$lag)

#### Feature Engineer - Create lagged features ####
# Define names of new features
lag_names <- colnames(lagCorr_sj)[c(2:21)]
lag_names <- c("precipitation_amt_mm_lag", "reanalysis_air_temp_k_lag", "reanalysis_avg_temp_k_lag", "reanalysis_dew_point_temp_k_lag", 
               "reanalysis_max_air_temp_k_lag", "reanalysis_min_air_temp_k_lag", "reanalysis_precip_amt_kg_per_m2_lag", "reanalysis_relative_humidity_percent_lag", 
               "reanalysis_sat_precip_amt_mm_lag", "reanalysis_specific_humidity_g_per_kg_lag", "reanalysis_tdtr_k_lag", "station_avg_temp_c_lag", 
               "avg_calc_lag", "station_diur_temp_rng_c_lag", "reanalysis_max_air_temp_k.1_lag", "reanalysis_min_air_temp_k.1_lag", 
               "station_max_temp_c_lag", "station_min_temp_c_lag", "station_precip_mm_lag", "total_cases_lag")

class(varLag_sj$lag)
class(varLag_iq$lag)

# Create and store lagged features in df
for (k in c(1:ncol(dengue_sj[,c(5:24)]))){
  variabel <- dengue_sj[,4+k]
  if (varLag_sj$lag[k] < 0){
  nas <- c(rep(NA, abs(as.numeric(varLag_sj$lag[k]))))
  rest <- variabel[c(1:(length(dengue_sj[,4+k])-length(nas)))]
  df.lagged <- data.frame(c(nas, rest))
  } else {
  nas <- c(rep(NA, abs(as.numeric(varLag_sj$lag[k]))))
  rest <- variabel[c((length(nas)+1):(length(dengue_sj[,4+k])))]
  df.lagged <- data.frame(c(rest, nas))
  }
  colnames(df.lagged)[1] <- lag_names[k]
  dengue_sj <- cbind(dengue_sj, df.lagged)
}

for (k in c(1:ncol(dengue_iq[,c(5:24)]))){
  variabel <- dengue_iq[,4+k]
  if (varLag_iq$lag[k] < 0){
    nas <- c(rep(NA, abs(as.numeric(varLag_iq$lag[k]))))
    rest <- variabel[c(1:(length(dengue_iq[,4+k])-length(nas)))]
    df.lagged <- data.frame(c(nas, rest))
  } else {
    nas <- c(rep(NA, abs(as.numeric(varLag_iq$lag[k]))))
    rest <- variabel[c((length(nas)+1):(length(dengue_iq[,4+k])))]
    df.lagged <- data.frame(c(rest, nas))
  }
  colnames(df.lagged)[1] <- lag_names[k]
  dengue_iq <- cbind(dengue_iq, df.lagged)
}

# Calculate corr matrix for lagged data
lag_corr_sj <- dengue_sj[c(11:917), c(25:38, 41:44)]
sum(is.na(lag_corr_sj))

lag_corr_iq <- dengue_iq[c(11:463), c(25:38, 41:44)]
sum(is.na(lag_corr_iq))

dengue_sj_lag_corr <- rcorr(as.matrix(lag_corr_sj))
dengue_iq_lag_corr <- rcorr(as.matrix(lag_corr_iq))  

# Visualization of correlation matrix 
corrplot(dengue_sj_lag_corr$r, method = "number",                                       
         addrect = 6, rect.col = "black",                                 
         type = "upper",                                                  
         tl.col="black", tl.srt=90, tl.cex = 0.65, tl.pos = "t",          
         addgrid.col = "lightgray")

corrplot(dengue_iq_lag_corr$r, method = "number",        
         addrect = 6, rect.col = "black",                
         type = "upper",                                 
         tl.col="black", tl.srt=90, tl.cex = 0.65, tl.pos = "t",      
         addgrid.col = "lightgray")


#### Normalization and Data-splitting
## Normalization
dengue_sj_n <- dengue_sj[c(11:917), c(3, 5:44)]
dengue_iq_n <- dengue_iq[c(11:453), c(3, 5:44)]
str(dengue_sj_n)
for (i in c(1:40)) {
  dengue_sj_n[, i] <- normalize(dengue_sj_n[, i], method = "standardize", range = c(0, 1), margin = 1L)
}
for (i in c(1:40)) {
  dengue_iq_n[, i] <- normalize(dengue_iq_n[, i], method = "standardize", range = c(0, 1), margin = 1L)
}

sj <- dengue_sj_n[, c(1:21, 41)]
sj_lag <- dengue_sj_n[, c(1, 22:41)]

iq <- dengue_iq_n[, c(1:21, 41)]
iq_lag <- dengue_iq_n[, c(1, 22:41)]


## Data splitting
set.seed(123)
inTrain_sj <- createDataPartition(y = sj_lag$total_cases_lag,      
                                  p = .80,                         
                                  list = FALSE)                    

inTrain_iq <- createDataPartition(y = iq_lag$total_cases_lag,      
                                  p = .80,                       
                                  list = FALSE)      

trainS_sj <- sj_lag[inTrain_sj,]
testS_sj  <- sj_lag[-inTrain_sj,]

trainS_iq <- iq_lag[inTrain_iq,]
testS_iq  <- iq_lag[-inTrain_iq,]

## Create dfs w/ independent(x) and independent(y) variable
x_sj <- trainS_sj[, c(1, 3, 5, 7, 11, 13, 18:19)]
y_sj <- trainS_sj[, 21]
y_sj <- as.numeric(y_sj)
class(y_sj)

x_sj_test <- testS_sj[, c(1, 3, 5, 7, 11, 13, 18:19)]
y_sj_test <- testS_sj[, 21]


x_iq <- trainS_iq[, c(1, 3, 5, 7, 11, 13, 18:19)]
y_iq <- trainS_iq[, 21]
y_iq <- as.numeric(y_iq)
class(y_iq)

x_iq_test <- testS_iq[, c(1, 3, 5, 7, 11, 13, 18:19)]
y_iq_test <- testS_iq[, 21]


#### Train, evaluate, forecast ####
## Train model - Selection
ctrl <- trainControl(method = "repeatedcv",        
                     repeats = 3,
                     number = 10,                  
                     allowParallel = TRUE,
                     savePredictions = "final"
)

models_list_sj <- caretEnsemble::caretList(x = x_sj,
                            y = y_sj,
                            trControl = ctrl,
                            methodList = c("rf", "glm", "gbm", "glmboost",  
                                           "treebag", 
                                           "parRF",                                 #pack: e1071
                                           "extraTrees",                            #pack: extraTrees 
                                           "svmLinear", "svmRadial", "svmRadial",
                                           "rotationForestCp",                      #pack: rotationForest 
                                           "rpart",
                                           "C5.0",                                  #pack: C50 
                                           "LogitBoost",                            #pack: caTools 
                                           "adaboost",                              #pack: fastAdaboost
                                           "randomGLM",                             #pack: randomGLM 
                                           "xgbTree"                                #pack: xgboost
                                           #"elm",                                  #pack:elmNN
                                           #"monmlp",
                                           #"dnn",                                  #deepnet
                                           #"nnet",
                                           #"pcaNNet"                               #nnet 
                                           ),     
                            preProc = c("center", "scale"),    # Center and scale predictors for training
                            metric = "MAE",
                            tuneList = NULL,
                            continue_on_fail = TRUE,
                            tuneLength = 3                     # Automatic Tuning Grid
)

models_list_iq <- caretEnsemble::caretList(x = x_iq,
                                           y = y_iq,
                                           trControl = ctrl,
                                           methodList = c("rf", 
                                                          #"glm", 
                                                          "gbm", #"glmboost",  
                                                          "treebag", 
                                                          #"parRF",                                 #pack: e1071
                                                          "extraTrees",                            #pack: extraTrees 
                                                          "svmLinear", "svmRadial", "svmRadial",
                                                          "rotationForestCp",                      #pack: rotationForest 
                                                          #"rpart",
                                                          "C5.0",                                  #pack: C50 
                                                          "LogitBoost",                            #pack: caTools 
                                                          "adaboost",                              #pack: fastAdaboost
                                                          "randomGLM",                             #pack: randomGLM 
                                                          "xgbTree"                                #pack: xgboost
                                                          #"elm",                                  #pack:elmNN
                                                          #"monmlp",
                                                          #"dnn",                                  #deepnet
                                                          #"nnet",
                                                          #"pcaNNet"                               #nnet 
                                           ),     
                                           preProc = c("center", "scale"),    # Center and scale predictors for training
                                           metric = "MAE",
                                           tuneList = NULL,
                                           continue_on_fail = TRUE,
                                           tuneLength = 5                     # Automatic Tuning Grid
)


## Train Ensemble
# Define trainControl and train selected models
stack_models_ctrl <- trainControl(method = "repeatedcv",        
                     repeats = 3,
                     number = 10,                  
                     allowParallel = TRUE,
                     savePredictions = "final"
)

stack_models_list_sj <- caretEnsemble::caretList(x = x_sj,
                                           y = y_sj,
                                           trControl = stack_models_ctrl,
                                           methodList = c("rf",
                                                          "svmLinear",
                                                          "svmRadial",
                                                          "svmPoly",
                                                          "rpart"),     
                                           preProc = c("center", "scale"),
                                           metric = "MAE",
                                           tuneList = NULL,
                                           continue_on_fail = TRUE
                                           #tuneLength = 2          
                                           
)

# Train caret Stack - Ensemble
stackControl <- trainControl(method = "repeatedcv", 
                             number = 10, repeats = 3, 
                             savePredictions = TRUE, 
                             classProbs = TRUE, 
                             verboseIter = TRUE,
                             tuneL)

stack <- caretStack(stack_models_list_sj, method = "rf", metric = "MAE", trControl = stackControl)


## Iter w/ stat._avg only
x_sj_2 <- x_sj[, c(6, 8)]
x_sj_test_2 <- x_sj_test[, c(6, 8)]

models_list_sj2 <- caretEnsemble::caretList(x = x_sj_2,
                                            y = y_sj,
                                            trControl = ctrl,
                                            methodList = c("rf", "glm", "gbm", "glmboost", "nnet", 
                                                           "treebag", 
                                                           "parRF",     
                                                           "extraTrees",
                                                           "svmLinear", "svmRadial", "svmRadial",
                                                           "rotationForestCp", 
                                                           "rpart",
                                                           "C5.0",          
                                                           "LogitBoost",     
                                                           "adaboost",      
                                                           "randomGLM",     
                                                           "xgbTree", 
                                                           "elm",
                                                           "monmlp",
                                                           "dnn",
                                                           "pcaNNet"), 
                                            preProc = c("center", "scale"),
                                            metric = "MAE",
                                            tuneList = NULL,
                                            continue_on_fail = TRUE
                                            #tuneLength = 2
                                           
)
# Results: No improvement, stick to the first selection of features


## Predict
l_pred_mlist_sj <- lapply(models_list_sj, function(x) predict(x, newdata = x_sj_test))
l_pred_mlist_iq <- lapply(models_list_iq, function(x) predict(x, newdata = x_iq_test))

l_pred_stack_mlist_sj <- lapply(stack_models_list_sj, function(x) predict(x, newdata = x_sj_test))

l_pred_mlist_sj2 <- lapply(models_list_sj2, function(x) predict(x, newdata = x_sj_test_2))

stack_pred <- predict(stack, newdata = x_sj_test)


## Evaluate
l_Acc_sj <- lapply(l_pred_mlist_sj, function(x) postResample(x, obs = y_sj_test))
#svmRadial
#RMSE   Rsquared        MAE 
#50.5499016  0.2531654 21.2976794 
l_Acc_iq <- lapply(l_pred_mlist_iq, function(x) postResample(x, obs = y_iq_test))
#svmRadial
#RMSE    Rsquared         MAE 
#10.77216434  0.04130063  5.85101611
#RMSE    Rsquared         MAE 
#10.57428669  0.09552265  5.64616292 

l_Acc_sj_mstack <- lapply(l_pred_stack_mlist_sj, function(x) postResample(x, obs = y_sj_test))

l_Acc_sj2 <- lapply(l_pred_mlist_sj2, function(x) postResample(x, obs = y_sj_test))   

postResample(stack_pred, y_sj_test)

## Gather performance metrics in one df
l <- NULL
df_Acc_sj <- NULL
for (l in c(1:length(l_Acc_sj))){
  model <- l_Acc_sj[[l]]
  df_Acc_sj <- rbind(df_Acc_sj, model)
}
colnames(df_Acc_sj) <- c("RMSE", "Rsquared", "MAE")
row.names(df_Acc_sj) <- names(l_Acc_sj)

m <- NULL
df_Acc_iq <- NULL
for (m in c(1:length(l_Acc_iq))){
  model_iq <- l_Acc_iq[[m]]
  df_Acc_iq <- rbind(df_Acc_iq, model_iq)
}
colnames(df_Acc_iq) <- c("RMSE", "Rsquared", "MAE")
row.names(df_Acc_iq) <- names(l_Acc_iq)

#### Preparing testset ####
## Load and check NAs
dengue_test_raw <- read.csv(file = "../dataset/dengue_test_raw_NA_flld.csv")
sum(is.na(dengue_test_raw))


## Split for each city 
dengue_test <- dengue_test_raw[, -c(5:6)]
dengue_test_sj <- dengue_test[c(which(dengue_test$city == "sj")), ]
dengue_test_iq <- dengue_test[c(which(dengue_test$city == "iq")), ]

## Select variable lags and store seperately
varLag_test_sj <- varLag_sj[-c(13:16, 19, 20) ,]
varLag_test_iq <- varLag_iq[-c(13:16, 19, 20) ,]

lag_names_test <- c("precipitation_amt_mm_lag", "reanalysis_air_temp_k_lag", "reanalysis_avg_temp_k_lag", "reanalysis_dew_point_temp_k_lag", 
               "reanalysis_max_air_temp_k_lag", "reanalysis_min_air_temp_k_lag", "reanalysis_precip_amt_kg_per_m2_lag", "reanalysis_relative_humidity_percent_lag", 
               "reanalysis_sat_precip_amt_mm_lag", "reanalysis_specific_humidity_g_per_kg_lag", "reanalysis_tdtr_k_lag", "station_avg_temp_c_lag", 
               "station_max_temp_c_lag", "station_min_temp_c_lag")
length(lag_names_test)

## Match datasets - Using same features as for training
lag_den_sj <- dengue_sj[, c(which(colnames(dengue_sj) %in% colnames(dengue_test_raw)))]
lag_den_iq <- dengue_iq[, c(which(colnames(dengue_iq) %in% colnames(dengue_test_raw)))]

length(which(colnames(lag_den_sj) %in% colnames(dengue_test)))

## Lag features of testset
k <- NULL
nas <- NULL
rest <- NULL
for (k in c(1:ncol(dengue_test_sj[,c(5:18)]))){
  variabel <- dengue_test_sj[,4+k]
  variabel2 <- lag_den_sj[,4+k]
  if (varLag_test_sj$lag[k] < 0){
    nas <- c(rep(NA, abs(as.numeric(varLag_test_sj$lag[k]))))
    lag_train <- variabel2[c((length(lag_den_sj[,4+k])-length(nas)+1):length(lag_den_sj[,4+k]))]
    rest <- variabel[c(1:(length(dengue_test_sj[,4+k])-length(nas)))]
    df.lagged <- data.frame(c(lag_train, rest))
  } else {                                                           # no else as there is no positiv lag
  }
  colnames(df.lagged)[1] <- lag_names_test[k]
  dengue_test_sj <- cbind(dengue_test_sj, df.lagged)
}

k <- NULL
nas <- NULL
rest <- NULL
df.lagged <- NULL
for (k in c(1:ncol(dengue_test_iq[,c(5:18)]))){
  variabel <- dengue_test_iq[,4+k]
  variabel2 <- lag_den_iq[,4+k]
  if (varLag_test_iq$lag[k] < 0){
    nas <- c(rep(NA, abs(as.numeric(varLag_test_iq$lag[k]))))
    lag_train <- variabel2[c((length(lag_den_iq[,4+k])-length(nas)+1):length(lag_den_iq[,4+k]))]
    rest <- variabel[c(1:(length(dengue_test_iq[,4+k])-length(nas)))]
    df.lagged <- data.frame(c(lag_train, rest))
  } else {
  }
  colnames(df.lagged)[1] <- lag_names_test[k]
  dengue_test_iq <- cbind(dengue_test_iq, df.lagged)
}


## Normalization
x_sj_test_test <- dengue_test_sj[, c(which(colnames(dengue_test_sj) %in% colnames(x_sj_test)))][, c(1:8)]
for (i in c(1:8)) {
  x_sj_test_test[, i] <- normalize(x_sj_test_test[, i], method = "standardize", range = c(0, 1), margin = 1L)
}

x_iq_test_test <- dengue_test_iq[, c(which(colnames(dengue_test_iq) %in% colnames(x_iq_test)))]
for (i in c(1:8)) {
  x_iq_test_test[, i] <- normalize(x_iq_test_test[, i], method = "standardize", range = c(0, 1), margin = 1L)
}

#### Model selected ####
# Define trainControl and train selected model
ctrl_fin <- trainControl(method = "repeatedcv",        
                         repeats = 5,
                         number = 10,                          
                         allowParallel = TRUE,
                         savePredictions = "final"
)

svmRad_sj <- caret::train(x = x_sj, 
                          y = y_sj, 
                          method = "svmRadial",
                          #tuneLength = 5,                
                          trControl = ctrl_fin,
                          preProc = c("center", "scale"),
                          metric = "MAE",
                          tuneList = NULL,
                          continue_on_fail = TRUE, 
                          tuneGrid = data.frame(.C = c(2.0), .sigma = c(.2))        # Identified set of parameters
)

svmRad_iq <- caret::train(x = x_iq, 
                          y = y_iq, 
                          method = "svmRadial",
                          #tuneLength = 5,
                          trControl = ctrl_fin,
                          preProc = c("center", "scale"), 
                          metric = "MAE",
                          tuneList = NULL,
                          continue_on_fail = TRUE,
                          tuneGrid = data.frame(.C = c(.5), .sigma = c(0.1356796))  # Identified set of parameters
)

## Predict
fin_pred_sj <- predict(svmRad_sj, newdata = x_sj_test)
#C=2.00, Sigma = const.:  41.32787  0.3903652  20.85323
#sigma = 0.2738562 and C = 2.

fin_pred_iq <- predict(svmRad_iq, newdata = x_iq_test)
#C=0.50. sigma = const.:   8.802692  0.15467538   5.367500
#sigma = 0.1356796 and C = 0.5.


## Evaluate
fin_Acc_svmRad_sj <- postResample(fin_pred_sj, obs = y_sj_test)
#C = 2, sigma = 0.2: 50.3051413  0.2458799 21.5154618 

fin_Acc_svmRad_iq <- postResample(fin_pred_iq, obs = y_iq_test)
#C = .5, sigma = 0.1356796: 10.77661834  0.04061223  5.84964049

#### Predict Test Set ####
fin_pred_sj_test_test <- predict(svmRad_sj, newdata = x_sj_test_test)
fin_pred_iq_test_test <- predict(svmRad_iq, newdata = x_iq_test_test)

PREDICTION <- data.frame(cbind(dengue_test_raw[, c(1:3)], c(fin_pred_sj_test_test, fin_pred_iq_test_test)))
write.csv(PREDICTION, file = "final_prediction_svmRad_compl2.csv")

#### End ####  
  