# RipRAM Random Forest Draft
# January 29, 2021
# Heili Lowman

# The following script will walk through a random forest created to predict state-wide RipRAM scores, with datasets from Kevin O'Connor, SMC, and StreamCat databases. The dependent variable in this case will be the Riparian Rapid Assessment Method (RipRAM) index state-wide.

# Thresholds from Kevin O'Connor 4/7/2021 (40, 60, 80) as default values

# Step One - Load In ------------------------------------------------------
setwd("L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/")

# install.packages("remotes")                                        # One way to get USGS package
# remotes::install_github("USGS-R/dataRetrieval")
# 
# library(remotes)
# install_github("USGS-R/dataRetrieval",                             # Another way to get it
#                build_opts = c("--no-resave-data", "--no-manual"),
#                build_vignettes = TRUE)

# Load packages. Which one prevents 'st_sfc(st_point' from working???  The abled packages below play nice with Spatial Join Method # 1-Alpha
library(quantregForest)
library(caret)
library(tidyverse)
library(tidymodels)
library(skimr)
library(sf)
# library(ggspatial)
library(nhdplusTools)
# library(patchwork)
library(Metrics)
# library(gt)
# library(sp)
# library(maptools)
# library(rgdal)
# library(dataRetrieval)



# Load datasets.

# Need to bind RipRAM data to NHD CA dataset to get COMIDs.
# Load in cleaned dataset from Kevin O'Connor's spreadsheet.
ripram_df <- read_csv("RipRAM_clean_012621.csv")
# Load in NHD_Plus_CA dataset from Annie as well as watersheds from Jeff.
# Full state of California
#nhd_ca <- read_sf("/Users/heilil/Desktop/hw_datasets/NHD_Plus_CA/NHDPlus_V2_FLowline_CA.shp") %>%


####
# # Spatial Join Method # 1-Alpha
# Use USGS 'discover_nhdplus_id' (not 'findNLDI' as in Option #3 below)

library(tidyverse)
library(sf)
library(nhdplusTools)
ripram_df2  <- ripram_df %>%
  rename(longitude=D_long, latitude=D_lat) %>%
  filter(!is.na(longitude)) %>%
  filter(!is.na(latitude)) %>%
  #arrange(SiteTag) %>%
  mutate(OneNum = 1,
         SiteNum = cumsum(OneNum)) %>%  # Create a running incremental count of sites.  Used as unique identifier in loop below.
  select(-OneNum)


ripram_df2 <- ripram_df2[!duplicated(ripram_df2$SiteNum), ]
Alist <- ripram_df2[!duplicated(ripram_df2$SiteNum), ]
Alist <- Alist$SiteNum

FirstTime <- 1
for (i in Alist){
  CurrentSite <- ripram_df2[ripram_df2$SiteNum == i, ]
  test <- st_sfc(st_point(c(CurrentSite$longitude, CurrentSite$latitude)), crs = 4326)  # Note, one of the packages prevents this from running
  test_comid <- discover_nhdplus_id(test)
  test_comid <- ifelse(is_empty(test_comid), -9999, test_comid) # this solves error when integer(empty)
  Score       <- data.frame(i, test_comid)
  if (FirstTime == 1){
    Score1 <- Score
  } else {
    Score1 <- rbind(Score1, Score)
  }
  FirstTime <- FirstTime + 1
}

#library (reshape)
#Score1  <-rename(Score1, c(i="SiteNum", test_comid="COMID"))       #old=new, using rename from reshape
Score1  <-dplyr::rename(Score1, c(SiteNum = i, COMID = test_comid)) # new=old, using rename from dplyr (opened with tidyverse)
#Score1 <- Score1[order(Score1$SiteNum), ]
#Bad.COMID <- Score1[Score1$COMID == -9999, ] # See how many bad sites (might be marine, with no actual stream)

ripram_comid <- ripram_df2 %>%
  inner_join(Score1, by=c("SiteNum"="SiteNum")) %>%
  mutate(COMID = as.numeric(COMID),  # bind said column to a new dataframe that will be used below
         RIPRAM = Idx_Score) %>% # and create newly named score column
  select(-SiteNum)

# Watershed characteristics' data available from StreamCat.
ca <- read_csv("streamcat_params.csv")
skim(ca) # uses library(skimr)
str(ca) # Checking to be sure COMID is numeric in both datasets.
# Perennial stream assessment data available from SCCWRP server.
ps6 <- read_csv("ps6_params.csv")
# In the ps6_rf_data script, if there are multiple Length_Fin measures for a given COMID, I have chosen the maximum of them and the associated PSA6 designation with that maximum.

# Bind the datasets together.
mydf <- ripram_comid %>%
  select(SiteTag, COMID, RIPRAM) %>% 
  inner_join(ca) %>% # Join with StreamCat watershed characteristics.
  inner_join(ps6) %>% # Join with PSA region dataset.
  select(-c(PctOpCat, PctOpWs, PctOpCatRp100, PctOpWsRp100, NPDESDensCat, 
            NPDESDensWs, TRIDensCat, TRIDensWs, SuperfundDensCat, SuperfundDensWs)) # Remove "open" land use and discharge site columns.
skim(mydf) # Examing completeness of this joined dataset.
length(unique(mydf$COMID)) # Checking for duplicates. 219 unique COMIDs. Not many...

# Pull out only one instance of each COMID.
set.seed(1) # Every time I run the code below, it's based on the same random pull of data.
mydf2 <- mydf %>% 
  group_by(COMID) %>%
  sample_n(size = 1) %>% 
  ungroup()

skim(mydf2) # Checking to make sure the dataset is complete.
# Important to have complete datasets for training data. For testing data, it's less critical.
#write.csv(ripram_comid, "L:/RipRAM_ES/Data/Working/RipRAM with COMID.csv")
#### END Spatial Join Method # 1-Alpha
library(patchwork)


# # # Spatial Join Method #1:
# nhd_ca <- read_sf("L:/RipRAM_ES/Data/Working/MapStuff/NHDPlus_NAD83.shp") %>%
#   mutate(COMID = as.numeric(COMID))
# 
# ripram_sf <- st_as_sf(ripram_df, # create sf compatible dataframe
#     coords = c("D_long", "D_lat"), # identify lon & lat
#     remove = F, # do not remove lat/lon columns
#     crs = 4269) # use NAD83 projection
# 
# nhd_z <- st_zm(nhd_ca)
# #
# EPSG <- make_EPSG() # create data frame of available EPSG codes.   Requires rgdal package
# EPSG[grepl("NAD83$", EPSG$note), ] # search for NAD 83 code
# # Units are in meters
# #
# # Join NHD dataset to RipRAM dataset.
# # Using this method, provided by O. Liu, since it's closest to the ArcGIS workflow.
# stream_samples_join <- st_join(nhd_ca, ripram_sf,                                          # error
#   join = st_is_within_distance, # the predicate to st_join
#   dist = 40) # joins samples within 40 m distance
# stream_remaining <- stream_samples_join %>%
#   filter(SiteTag != "NA")
# 
# # Spatial Join Method #2:
# # Another option is to make polygon buffer from lines.
# # This yielded a lot of duplicate matches (150+) so I chose to go with the above workflow.
# nhd_z <- st_zm(nhd_ca)
# streams_buff <- st_buffer(nhd_z, dist = 0.001) # width/thickness = 0.001 degrees
# ripram_sf <- st_as_sf(ripram_df, # create sf compatible dataframe
#     coords = c("D_long", "D_lat"), # identify lon & lat
#     remove = F, # do not remove lat/lon columns
#     crs = 4269) # use NAD83 projection
# streams_samples_join <- st_join(streams_buff, ripram_sf) # Joins points to polygons made in the line above.
# streams_remaining <- streams_samples_join %>%
#  filter(SiteTag != "NA")
# #
# # Exporting dataset for Annie to double-check in ArcGIS. It matched pretty well but went with method #3 to avoid issues in small catchments that may arise from snapping to flowlines.
# as_tibble(stream_remaining) %>% # remove geometry
#   select(c(COMID, D_lat, D_long)) %>% # select necessary columns
#   write_csv("ripram_sites.csv") # export as .csv
# 
# # Spatial Join Method #3:
# # I'll be using the USGS' dataRetrieval package function findNLDI(), which uses the NHD catchment boundaries instead of flowlines.
# test_comid <- findNLDI(location = c(-89.362239, 43.090266)) # testing to be sure this generates COMID 13293750
# 
# test <- findNLDI(location = c(ripram_df$D_long[1], ripram_df$D_lat[1])) # this works
# 
# # After a lot of testing, it seems the findNLDI function can only take one input at a time, so I need to figure out a way to map all of the rows into it one by one.
# 
# long <- ripram_df$D_long # vector of longitudes
# lat <- ripram_df$D_lat # vector of latitudes
# new_df <- data.frame(matrix(ncol = 360, nrow = 1)) # empty vector to receive the data
# 
# # for loop
# for(k in 1:360) { # iterate from 1 to 360
#   new_df[k] <- as.data.frame.list(unlist(findNLDI(location = c(long[k], lat[k]), no_sf = TRUE))) %>%
#     pivot_wider(names_from = sourceName, values_from = comid)
#   # save input of findNLDI to output vector new_df for every iteration of k
#   # needs to be unlisted because of findNLDI output structure
# }
# 
# # I couldn't quite figure out how to get the dataset to format properly, so I'm creating a WIDE dataframe in the for loop above. I QAQC'ed a few by hand, and they are indeed printing out the COMIDs.
# 
# ripram_comids_beta <- new_df %>% # take the wide dataset
#   pivot_longer(cols = starts_with("X"), 
#     names_to = "index",
#     values_to = "comid") # pivot this into a column
# 
# ripram_comid <- ripram_df %>%
#   mutate(COMID = as.numeric(ripram_comids_beta$comid),  # bind said column to a new dataframe that will be used below
#         RIPRAM = Idx_Score) # and create newly named score column
# # Watershed characteristics' data available from StreamCat.
# ca <- read_csv("streamcat_params.csv")
# skim(ca)
# str(ca) # Checking to be sure COMID is numeric in both datasets.
# # Perennial stream assessment data available from SCCWRP server.
# ps6 <- read_csv("ps6_params.csv")
# # In the ps6_rf_data script, if there are multiple Length_Fin measures for a given COMID, I have chosen the maximum of them and the associated PSA6 designation with that maximum.
# 
# # Bind the datasets together.
# mydf <- ripram_comid %>%
#   select(SiteTag, COMID, RIPRAM) %>% 
#   inner_join(ca) %>% # Join with StreamCat watershed characteristics.
#   inner_join(ps6) %>% # Join with PSA region dataset.
#   select(-c(PctOpCat, PctOpWs, PctOpCatRp100, PctOpWsRp100, NPDESDensCat, 
#     NPDESDensWs, TRIDensCat, TRIDensWs, SuperfundDensCat, SuperfundDensWs)) # Remove "open" land use and discharge site columns.
# skim(mydf) # Examing completeness of this joined dataset.
# length(unique(mydf$COMID)) # Checking for duplicates. 219 unique COMIDs. Not many...
# 
# # Pull out only one instance of each COMID.
# set.seed(1) # Every time I run the code below, it's based on the same random pull of data.
# mydf2 <- mydf %>% 
#   group_by(COMID) %>%
#   sample_n(size = 1) %>% 
#   ungroup()
# 
# skim(mydf2) # Checking to make sure the dataset is complete.
# # Important to have complete datasets for training data. For testing data, it's less critical.

# Step Two - Training Data ------------------------------------------------

# Create calibration and validation splits with tidymodels initial_split() function.

set.seed(4)
mydf2_split <- mydf2 %>%
  initial_split(prop = 0.75, strata = PSA6) # splits data into training and testing set.
# default is 3/4ths split (but 75% training, 25% testing).
# Stratification (strata) = grouping training/testing sets by region, state, etc.
# Using the "strata" call ensures the number of data points in the training data is equivalent to the proportions in the original data set. (Strata below 10% of the total are pooled together.)

# Create a training data set with the training() function
# Pulls from training and testing sets created by initial_split()
mydf2_train <- training(mydf2_split)
mydf2_test <- testing(mydf2_split)
# Examine the environment to be sure # of observations looks like the 75/25 split. 165:54.

# Create a separate dataset of available COMIDS that were not used in the training dataset.
nottrain <- ca %>% # all COMIDS from StreamCat data, sampled or not
  filter(!COMID %in% mydf2_train$COMID) # Removing sites used to train the model. n = 140,545

# Step Three - Kitchen Sink model -----------------------------------------

# Create finalized training dataset and include all possible variables. 
rf_dat <- mydf2_train %>%
  select(-SiteTag, -COMID, -PSA6, -Length_Fin)

# Random forest -- 
# a decision tree model, using predictors to answer dichotomous questions to create nested splits.
# no pruning happens - rather, multiple trees are built (the forest) and then you are looking for consensus across trees
# training data goes down the tree and ends up in a terminal node.
# if testing data goes down the same route, then this upholds our conclusions. Or, if it goes awry, this allows us to look for patterns in how it goes awry.

set.seed(2) # assures the data pulled is random, but sets it for the run below (makes outcome stable)
myrf <- randomForest(y = rf_dat$RIPRAM, # dependent variable
  x = rf_dat %>%
    select(-RIPRAM), # selecting all predictor variables
  importance = T, # how useful is a predictor in predicting values (nothing causal)
  proximity = T, 
  ntrees = 500) # 500 trees. 

myrf # examine the results.
# 37.69% variance explained.

summary(myrf)
# mtry allows you to parameterize the number of splits

plot(myrf)
# model performance appears to improve most at ~300 trees

varImpPlot(myrf)
# displays which variables are most important
# helps to winnow down list of predictors
# recommended to weigh left pane more
# right pane also shows how evenly things split based on the list of predictors
# values close to 0 can be dropped, but don't have to be
# road, mine, and dam density appear to have the greatest impact

importance <- myrf$importance
View(importance)
# displays the data plotted in the plot above

# predict()
# returns out of bag predictions for training data
# in the bag: every time a tree is built, it uses ~80% of the original 75% we set aside from the original dataset used to create a tree to assure random data selection
# out of bag: looking at the remaining 20% of the training data to predict, when you want to know what your model does at the training location sites

# Predict RIPRAM scores state-wide for all COMIDs.
nottrain_prediction <- nottrain %>% # taking all available COMIDS, that haven't been used in training
  na.omit() %>% # remove NAs
  mutate(ripram_predicted = predict(myrf, newdata = nottrain %>% na.omit())) # using developed model (myrf), inputting predictor variables (nottrain - which contains COMIDs and associated StreamCat data) to predict output/dependent variable (ripram_predicted a.k.a. RIPRAM).

# rePredict RIPRAM scores for training data.
mydf2_train$ripram_predicted <- predict(myrf) # Add column of predicted RIPRAM values to training dataset.

# Creates new dataset of bound rows for both ...
ca_predictions <- bind_rows(nottrain_prediction %>%
    mutate(Set = "Non-training"), # statewide COMIDs (not used for training)
  mydf2_train %>%
    mutate(Set = "Training")) # COMIDS from training dataset
# This creates the dataset that will be plotted to create a state-wide plot of predicted CRAM scores.

# Plot the data.
rf_plot1 <- ggplot(ca_predictions, aes(x = PctImp2011CatRp100, y = ripram_predicted)) +
  geom_point(alpha = 0.1) +
  labs(x = "Mean % imperviousness within catchment and within a 100-m buffer of NHD stream lines",
    y = "Predicted RIPRAM Score") +
  theme_classic() +
  facet_wrap(.~Set)

rf_plot1

# Step Four - Predictor Selection -----------------------------------------

# Using caret to select the best predictors
# What are the parameters you want to use to run recursive feature elimination (rfe)?
my_ctrl <- rfeControl(functions = rfFuncs,
  method = "cv",
  verbose = FALSE,
  returnResamp = "all")

# rfe = recursive feature elimination
# THIS STEP TAKES FOR-EV-ER!!!
set.seed(22)
my_rfe <- rfe(y = rf_dat$RIPRAM, # set dependent variable
  x = rf_dat %>% select(-RIPRAM), # set predictor variables
  size = c(3:10, 15, 20, 25, 30), # sets how many variables are in the overall model
  # I have 34 total possible variables, so I've chosen increments of 5 to look at.
  rfeControl = my_ctrl) # pull in control from above

# can you make your model even simpler?
# the following will pick a model with the smallest number of predictor variables based on the tolerance ("tol") that you specify (how much less than the best are you willing to tolerate?)
my_size <- pickSizeTolerance(my_rfe$results, metric = "RMSE", tol = 1, maximize = F)
# higher tol (~10) gives you less variables
# lower tol (~1) gives you more variables - "I'd like the simplest model within 1% of the best model."
pickVars(my_rfe$variables, size = my_size)

# pickVars (20):
# "PctImp2011WsRp100"  "PctImp2011CatRp100" "PctAgCat"           "FertWs"             "CanalDensWs"        "PctAgWs"           
# "PctImp2011Ws"       "PctImp2011Cat"      "FertCat"            "PctUrbCatRp100"     "PctAgCatRp100"      "PctAgWsRp100"      
# "CBNFCat"            "RdCrsWs"            "ManureWs"           "AgKffactCat"        "AgKffactWs"         "PctUrbWsRp100"     
# "CBNFWs"             "PctUrbCat"

# Proceed with a regular RF that yields mean weighted values and fit those into the following classification scheme:

#Likely condition approach: Compare mean to three RIPRAM thresholds (??, ??, ??) based on condition classes.
# Very likely altered: mean < ??
# Likely altered: mean < ??
# Possibly altered: mean < ??
# Likely unaltered: mean >= ??

# Predict scores using the above 20 variables:

# Create re-finalized training dataset and include all possible variables. 
rf_dat2 <- mydf2_train %>%
  select(RIPRAM, PctImp2011WsRp100, PctImp2011CatRp100, PctAgCat, FertWs, CanalDensWs, PctAgWs, PctImp2011Ws, PctImp2011Cat, FertCat, PctUrbCatRp100, PctAgCatRp100, PctAgWsRp100, CBNFCat, RdCrsWs, ManureWs, AgKffactCat, AgKffactWs, PctUrbWsRp100, CBNFWs, PctUrbCat)

set.seed(4) # assures the data pulled is random, but sets it for the run below (makes outcome stable)
myrf2 <- randomForest(y = rf_dat2$RIPRAM, # dependent variable
  x = rf_dat2 %>%
    select(-RIPRAM),
  importance = T, 
  proximity = T, 
  ntrees = 500)  

myrf2 # examine the results. 
# 38.97% variance explained.
summary(myrf2)
plot(myrf2) # need min of 200 trees.
varImpPlot(myrf2)

importance2 <- as.data.frame(as.table(myrf2$importance))
View(importance2) # displays the data plotted in the plot above

# Nicer ggplot variable importance plot.
vip_plot_a <- importance2 %>%
  filter(Var2 == "%IncMSE") %>%
  mutate(Var1 = factor(Var1)) %>%
  mutate(Var1_f = fct_reorder(Var1, Freq)) %>%
  ggplot(aes(x = Freq, y = Var1_f)) +
  geom_point(size = 3, alpha = 0.75) +
  labs(x = "% Importance (MSE)",
    y = "Variables") +
  theme_bw()

vip_plot_b <- importance2 %>%
  filter(Var2 == "IncNodePurity") %>%
  mutate(Var1 = factor(Var1)) %>%
  mutate(Var1_f = fct_reorder(Var1, Freq)) %>%
  ggplot(aes(x = Freq, y = Var1_f)) +
  geom_point(size = 3, alpha = 0.75) +
  labs(x = "Node Purity",
    y = "Variables") +
  theme_bw()

vip_plot <- vip_plot_a + vip_plot_b

vip_plot

# png(file="ripram_vip_plot.png", units="in", width=8, height=5, res=300)
# vip_plot
# dev.off()

# ggsave("ripram_vip_plot.png",
#      path = "/Users/heilil/Desktop/R_figures",
#      width = 25,
#      height = 10,
#      units = "cm"
#    )

# predict(myrf2) # returns out of bag predictions for training data

# Predict RIPRAM scores state-wide.
nottrain_prediction2 <- nottrain %>% # taking all COMIDS that haven't been used in training
  na.omit() %>% # remove NAs
  mutate(ripram_predicted = predict(myrf2, newdata = nottrain %>% na.omit())) # using developed model (myrf2), inputting predictor variables (nottrain - COMIDs and associated StreamCat data) to predict output/dependent variable (ripram_predicted a.k.a. RIPRAM).

# rePredict RIPRAM scores for training and testing data (to be used in validation below).
mydf2_train2 <- mydf2_train
mydf2_train2$ripram_predicted <- predict(myrf2) # Add column of predicted RIPRAM scores to training dataset.

mydf2_test2 <- mydf2_test %>%
  mutate(ripram_predicted = predict(myrf2, newdata = mydf2_test %>% select(-c(SiteTag, RIPRAM, PSA6, Length_Fin)))) # Adds column of predicted RIPRAM values to testing dataset.

# Creates new dataset of bound rows for both ...
ca_predictions2 <- bind_rows(nottrain_prediction2 %>%
    mutate(Set = "Non-training"), # statewide COMIDs (not used for training data)
  mydf2_train2 %>%
    mutate(Set = "Training")) # COMIDS from training dataset (used for training the model).
# This creates the dataset that will be plotted.

# Create table of number of sites that fall into each category.
# NEEDS TO BE REDONE WITH ACTUAL CLASSIFICATIONS!
# Add classification column.
# ca_predictions2 <- ca_predictions2 %>%
#   mutate(classification = case_when(ripram_predicted < 60 ~"Very Likely Altered",
#     ripram_predicted >= 60 & ripram_predicted < 75 ~"Likely Altered",
#     ripram_predicted >= 75 & ripram_predicted < 90 ~"Possibly Altered",
#     ripram_predicted >= 90 ~"Likely Unaltered")) %>%
#   mutate(class_f = factor(classification, levels = c("Very Likely Altered", "Likely Altered", "Possibly Altered", "Likely Unaltered"))) # relevel classifications
ca_predictions2 <- ca_predictions2 %>%
  mutate(classification = case_when(round(ripram_predicted, digits = 0) < 40 ~"Very Likely Altered",
                                    round(ripram_predicted, digits = 0) < 60 ~"Likely Altered",
                                    round(ripram_predicted, digits = 0) < 80 ~"Possibly Altered",
                                    round(ripram_predicted, digits = 0) >= 80 ~"Likely Unaltered")) %>%
  mutate(class_f = factor(classification, levels = c("Very Likely Altered", "Likely Altered", "Possibly Altered", "Likely Unaltered"))) # relevel classifications

#### Results .csv ####
# Export results.
#write_csv(ca_predictions2, "ripram_rf_results.csv")
#write_csv(ca_predictions2, "L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/ripram_rf_results.csv")

# Summary table by site #.
ca_summary <- ca_predictions2 %>%
  count(class_f) # count sites statewide by classification
# The numbering is greatly skewed to the "possibly altered" classification, so perhaps other thresholds are necessary.

# Summary table by stream length (m)
ca_summary_length <- ca_predictions2 %>%
  group_by(class_f) %>% # group by classification
  summarize(length = sum(Length_Fin, na.rm=TRUE)) # sum stream lengths

# Join and export.
ca_sum <- full_join(ca_summary, ca_summary_length)
#write_csv(ca_sum, "ripram_rf_results_summary.csv")
#write_csv(ca_sum, "L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/ripram_rf_results_summary.csv")

# Step Five - Quantile Regression model -----------------------------------

# Note - for the Healthy Watersheds Project, I did not pursue this structure, but I've kept some example code below in case future iterations call for it.

# Quantile random forest regression mode, instead of looking at the mode of trees, can compare to 10th, 50th, 90th percentiles etc.

# Need to make a new dataset taking the above results of pickVars into account.
# Create finalized training dataset and include all possible variables. 
# qrf_dat <- mydf2_train %>%
#   select(asci, RdCrsWs, PctAgWs, PctUrbWsRp100, PctOpWsRp100, PctOpWs, DamDensWs, RdDensWs, NABD_DensWs, PctUrbWs, PctUrbCatRp100, RdDensWsRp100, PctOpCat, PctUrbCat, RdDensCat, CBNFWs, PctOpCatRp100, PctAgWsRp100, TRIDensWs, AgKffactWs, FertWs) 

# set.seed(20)
# myqrf <- quantregForest(y = qrf_dat$asci, # dependent variable
#               x = qrf_dat %>%
#                   select(-asci),
#               importance = T, 
#               proximity = T,
#               keep.inbag=T,
#               ntrees = 500) 

#predict(myqrf) # automatically presents 10th %tile, median, and 90th %tile
#predict(myqrf, what=c(0.2, 0.3, 0.999)) # to print specific quantiles

#plot(myqrf) # plots the results.
# Again appears to improve after ~100 trees.

# Step Six - Model validation ---------------------------------------------

# Compare predicted vs. actual results, including by PSA region.
# Adding lines of slope=1 and linear models to each plot.
val1 <- ggplot(mydf2_train2, aes(x = ripram_predicted, y = RIPRAM)) +
  geom_point(color = "#2A3927", alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#2A3927") +
  labs(x = "RIPRAM predicted",
    y = "RIPRAM measured",
    title = "Training Data\nn=165") +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()

val1

lm1 <- lm(RIPRAM~ripram_predicted, data = mydf2_train2)
summary(lm1)

val2 <- ggplot(mydf2_test2, aes(x = ripram_predicted, y = RIPRAM)) +
  geom_point(color = "#3793EC", alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#3793EC") +
  scale_x_continuous(breaks = c(0.5, 0.7, 0.9)) +
  labs(x = "RIPRAM predicted",
    y = "RIPRAM measured",
    title = "Testing Data\nn=54") +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()

val2

lm2 <- lm(RIPRAM~ripram_predicted, data = mydf2_test2)
summary(lm2)

# Create the full testing + training dataset to plot together.
mydf2_test2$set <- "Testing"
mydf2_train2$set <- "Training"
full_train_test <- bind_rows(mydf2_test2, mydf2_train2) %>%
  mutate(set_f = factor(set, levels = c("Training", "Testing")))

val3 <- ggplot(full_train_test, aes(x = ripram_predicted, y = RIPRAM, color = set_f)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(name = "Set", values = c("#2A3927", "#3793EC"), drop = FALSE) +
  labs(x = "RIPRAM predicted",
    y = "RIPRAM measured",
    title = "All Data\nn=219") +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  facet_wrap(~PSA6) +
  theme_bw()

val3

val_fig <- (val1 + val2) /
  (val3)

val_fig + plot_annotation(
  title = 'RIPRAM Random Forest Results',
  subtitle = 'All modeling performed using StreamCAT datasets.',
  caption = 'Linear models are colored according to dataset. Lines of slope = 1 are denoted in black.'
)

# Save figure.
# ggsave("ripram_rfmodel_validation.png",
#      path = "/Users/heilil/Desktop",
#      width = 35,
#      height = 25,
#      units = "cm"
#    )

# Chose not to compute confusion matrix / accuracy score since this is more applicable to categorical ouputs from random forest models -
# Instead, calculated Root Mean Squared Error (RMSE) of both training and test datasets.
# If test RMSE values are much greater than training, then possible the model has been over fit.

predtest <- predict(myrf2, mydf2_test2)
rmse(mydf2_test2$RIPRAM, predtest)
# 17.63

predtrain <- predict(myrf2, mydf2_train2)
rmse(mydf2_train2$RIPRAM, predtrain)
# 8.07

# Double checking using the original random forest dataset (rf_dat) with all 35 possible variables included to see where the error in number of predictors starts to increase dramatically (to help double check our decision to include 25 parameters).
dc <- rfcv(rf_dat %>%
    select(-RIPRAM), 
  rf_dat$RIPRAM,
  step = 0.7, # default is 0.5
  scale="log")

dc$error.cv
#     34       24       17       12        8        6        4        3        2        1 
# 370.0889 366.7509 386.6936 409.5337 441.2518 460.4203 513.3549 500.2896 536.4003 666.5010

# Appears between 24 and 17 variables, there is a significant increase in error. This model is roughly the same size as the CSCI (20) model.

# Step Seven - Map results state-wide -------------------------------------

# Using ca_predictions2 dataset generated above. But need to first associate lat/lon with each COMID.

# Load in NHD_Plus_CA dataset from Annie as well as watersheds from Jeff.
# Full state of California
#nhd_ca <- read_sf("/Users/heilil/Desktop/hw_datasets/NHD_Plus_CA/NHDPlus_V2_FLowline_CA.shp") %>%
nhd_ca <- read_sf("L:/RipRAM_ES/Data/Working/MapStuff/NHDPlus_NAD83.shp") %>%
  mutate(COMID = as.numeric(COMID))

# South Coast watersheds - Ventura River, San Juan Creek, San Diego River
#nhd_vr <- read_sf("/Users/heilil/Desktop/hw_datasets/NHD_Watersheds/VenturaRiver_NHD_Clip.shp") %>%
nhd_vr <- read_sf("L:/RipRAM_ES/Data/Working/MapStuff/VenturaRiver_NHD_Clip.shp") %>%
  mutate(COMID = as.numeric(COMID))

#nhd_sjc <- read_sf("/Users/heilil/Desktop/hw_datasets/NHD_Watersheds/SanJuanCreek_NHD_Clip.shp") %>%
nhd_sjc <- read_sf("L:/RipRAM_ES/Data/Working/MapStuff/SanJuanCreek_NHD_Clip.shp") %>%
  mutate(COMID = as.numeric(COMID))

#nhd_sdr <- read_sf("/Users/heilil/Desktop/hw_datasets/NHD_Watersheds/SanDiegoRiver_NHD_Clip.shp") %>%
nhd_sdr <- read_sf("L:/RipRAM_ES/Data/Working/MapStuff/SanDiegoRiver_NHD_Clip.shp") %>%
  mutate(COMID = as.numeric(COMID))

# Assign modeled COMIDs to mcomid.
mcomid <- ca_predictions2$COMID

# Filter by and plot only modeled stream reaches.

modeled_ripram_map <- nhd_ca %>%
  filter(COMID %in% mcomid) %>%
  inner_join(ca_predictions2) %>%
  ggplot() +
  geom_sf(aes(color = class_f)) +
  scale_color_manual(name = "Condition", values = c("red2", "lightpink", "lightskyblue2", "steelblue"), drop = FALSE) +
  theme_bw()

modeled_ripram_map
# Note, sometimes this takes forever to render in the "plot" pane.
# Best to just save to your machine (below) and then take a look.

# png(file="ripram_modeled_CA.png", units="in", width=8, height=5, res=300)
# png(file="L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/ripram_modeled_CA.png", units="in", width=8, height=5, res=300)
# modeled_ripram_map
# dev.off()

# ggsave("ripram_modeled_CA.png",
#        width = 35,
#        height = 35,
#        units = "cm"
# )


# Additional Notes - Healthy Watersheds project ---------------------------

# Condition approach favored by Anna per meeting on 11/3/2020.

# Works Cited:

# Hill, Ryan A., Marc H. Weber, Scott G. Leibowitz, Anthony R. Olsen, and Darren J. Thornbrugh, 2016. The Stream-Catchment (StreamCat) Dataset: A Database of Watershed Metrics for the Conterminous United States. Journal of the American Water Resources Association (JAWRA) 52:120-128. DOI: 10.1111/1752-1688.12372.

# End of R script.
