# CSCI Random Forest Draft
# October 21, 2020
# Heili Lowman

# Note '_V2' of this script uses published thresholds from Mazor et al 2016.

# The following script will walk through a random forest created to predict state-wide CSCI scores, with datasets from SMC and StreamCat databases. The dependent variable in this case will be the California Stream Condition Index (CSCI) state-wide.

# Step One - Load In ------------------------------------------------------
setwd("L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/")             # saved locally

# Load packages.
library(quantregForest)
library(caret)
library(tidyverse)
library(tidymodels)
library(skimr)
library(sf)
library(ggspatial)
library(nhdplusTools)
library(patchwork)
library(Metrics)
library(gt)

# Load datasets.
# CSCI data available from SCCWRP database.
csci_df <- read_csv("csci_rf_data1.csv") %>% # Loads in dataset pulled 12/3/20
  rename(COMID = comid) %>%
  drop_na(csci) # Drop the NA values.
skim(csci_df) # Examine the dataset.
str(csci_df)

# Watershed characteristics' data available from StreamCat.
ca <- read_csv("streamcat_params.csv")
skim(ca)
str(ca) # Checking to be sure COMID is numeric in both datasets.
# Perennial stream assessment data available from SCCWRP server.
ps6 <- read_csv("ps6_params.csv")
# In the ps6_rf_data script, if there are multiple Length_Fin measures for a given COMID, I have chosen the maximum of them and the associated PSA6 designation with that maximum.

# Bind the datasets together.
# From here on out, naming conventions will be identical to the asci script (for ease of duplication). But be aware of this when running both scripts.
mydf <- csci_df %>%
  select(stationcode, COMID, csci) %>%
  inner_join(ca) %>%
  inner_join(ps6) %>%
  select(-c(PctOpCat, PctOpWs, PctOpCatRp100, PctOpWsRp100, NPDESDensCat, NPDESDensWs, TRIDensCat, TRIDensWs, SuperfundDensCat, SuperfundDensWs)) # Remove "open" land use and discharge site columns.
skim(mydf) # Examing completeness of this joined dataset.
length(unique(mydf$COMID)) # Checking for duplicates. 2020 unique COMIDs.

# Pull out only one instance of each COMID.
set.seed(1) # Ensures each run of code below is based on same random pull of data.
mydf2 <- mydf %>% 
  filter(stationcode!="109PS0162") %>% #There's only one site missing RdDensCatRp100. Better to drop the site than to drop the metric
  group_by(COMID) %>%
  sample_n(size = 1) %>% 
  ungroup()

# Important to have complete datasets for training data. For testing data, it's less critical.

# Step Two - Training Data ------------------------------------------------

# Create calibration and validation splits with tidymodels initial_split() function.

set.seed(4)
mydf2_split <- mydf2 %>%
  initial_split(prop = 0.75, strata = PSA6) # splits data into training and testing set.
# default is 3/4ths split (but 75% training, 25% testing).
# Stratification (strata) = grouping training/testing sets by region, state, etc.

# Create a training data set with the training() function
# Pulls from training and testing sets created by initial_split()
mydf2_train <- training(mydf2_split)
mydf2_test <- testing(mydf2_split)
# Examine the environment to be sure # of observations looks like the 75/25 split. 1516:503.

# Create a separate dataset of available COMIDS that were not used in the training dataset.
nottrain <- ca %>% # all COMIDS from StreamCat data, sampled or not
  filter(!COMID %in% mydf2_train$COMID) # Removing sites used to train the model. n = 139194

# Step Three - Kitchen Sink model -----------------------------------------
# Checked to be sure all variables used by Marcus Beck to create the scape tool are included in this study:
# "CanalDensCat","CanalDensWs", "PctImp2011Cat", "PctImp2011Ws", "PctImp2011CatRp100", "PctImp2011WsRp100", "PctUrbCat","PctUrbWs","PctAgCat","PctAgWs", "PctUrbCatRp100","PctUrbWsRp100","PctAgCatRp100","PctAgWsRp100", "RdDensCat","RdDensWs","RdDensCatRp100","RdDensWsRp100", "RdCrsCat","RdCrsWs"

# Create finalized training dataset and include all possible variables. 
rf_dat <- mydf2_train %>%
  select(-stationcode, -COMID, -PSA6, -Length_Fin)

# Random forest -- 
# a decision tree model, using predictors to answer dichotomous questions to create nested splits.
# no pruning happens - rather, multiple trees are built (the forest) and then you are looking for consensus across trees
# training data goes down the tree and ends up in a terminal node.
# if testing data goes down the same route, then this upholds our conclusions. Or, if it goes awry, this allows us to look for patterns in how it goes awry.

set.seed(2) # assures the data pulled is random, but sets it for the run below (makes outcome stable)
myrf <- randomForest(y = rf_dat$csci, # dependent variable
  x = rf_dat %>%
    select(-csci), # selecting all predictor variables
  importance = T, # how useful is a predictor in predicting values (nothing causal)
  proximity = T, 
  ntrees = 500) # 500 trees. 

myrf # examine the results.
# 54.37% variance explained.

summary(myrf)
# mtry allows you to parameterize the number of splits

plot(myrf)
# model performance appears to improve most at ~75 trees

varImpPlot(myrf)
# displays which variables are most important
# helps to winnow down list of predictors
# recommended to weigh left pane more
# right pane shows how evenly things split based on the list of predictors
# values close to 0 can be dropped, but don't have to be

# In both panes - urban/ag/impervious land cover and stream-road crossings appear important.

importance <- myrf$importance
View(importance)
# displays the data plotted in the plot above

predict(myrf)
# returns out of bag predictions for training data
# in the bag: every time a tree is built, it uses ~80% of the original 75% we set aside from the original dataset used to create a tree to assure random data selection
# out of bag: looking at the remaining 20% of the training data to predict, when you want to know what your model does at the training location sites

# Predict CSCI scores state-wide.
nottrain_prediction <- nottrain %>% # taking all available COMIDS, that haven't been used to train the model
  na.omit() %>% # remove NAs
  mutate(csci_predicted = predict(myrf, newdata = nottrain %>% na.omit())) # using developed model (myrf), inputting predictor variables (nottrain - which contains COMIDs and associated StreamCat data) to predict output/dependent variable (csci_predicted a.k.a. CSCI).

# rePredict CSCI scores for training data.
mydf2_train$csci_predicted <- predict(myrf) # Adds column of predicted CSCI values to training dataset.

# Creates new dataset of bound rows for both ...
ca_predictions <- bind_rows(nottrain_prediction %>%
                            mutate(Set = "Non-training"), # statewide COMIDs (not used for training data)
                            mydf2_train %>%
                            mutate(Set = "Training")) # COMIDS from training dataset (that were used for training the random forest model).
# This creates the dataset that will be plotted (i.e. you're trying to create a state-wide plot of predicted CSCI scores).

# Plot the data.
rf_plot1 <- ggplot(ca_predictions, aes(x = RdCrsWs, y = csci_predicted)) +
  geom_point(alpha = 0.1) +
  labs(x = "Watershed Density of Road-Stream Intersections",
    y = "Predicted CSCI Score") +
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
my_rfe <- rfe(y = rf_dat$csci, # set dependent variable
              x = rf_dat %>% select(-csci), # set predictor variables
              size = c(3:10, 15, 20, 25, 30), # sets how many variables are in the overall model - I have 34 total possible variables.
              rfeControl = my_ctrl) # pull in control from above

# can you make your model even simpler?
# the following will pick a model with the smallest number of predictor variables based on the tolerance ("tol") that you specify (how much less than the best are you willing to tolerate?)
my_size <- pickSizeTolerance(my_rfe$results, metric = "RMSE", tol = 1, maximize = F)
# higher tol (~10) gives you less variables
# lower tol (~1) gives you more variables - "I'm taking the simplest model that's within 1% of the best model."
pickVars(my_rfe$variables, size = my_size)

# pickVars (20): "PctImp2011Cat" "PctAgWs" "PctImp2011Ws" "PctImp2011CatRp100"
# "PctUrbWsRp100" "PctImp2011WsRp100" "RdDensWs" "PctAgWsRp100"      
# "PctUrbCatRp100" "PctUrbWs" "AgKffactWs" "RdDensCatRp100"    
# "FertWs" "CanalDensWs" "RdDensCat" "PctUrbCat"         
# "RdDensWsRp100" "RdCrsWs" "CBNFWs" "ManureWs" 

# Proceed with a regular RF that yields mean weighted values and fit those into the following classification scheme:

# Mazor et al. 2016 1st percentile = 0.63, not 0.61 !!!!

#Likely condition approach: Compare mean to three CSCI thresholds (0.61, 0.79, 0.92) @ reference sites (1st, 10th, 30th percentiles)
# Very likely altered: mean < 0.61
# Likely altered: mean < 0.79
# Possibly altered: mean < 0.92
# Likely unaltered: mean >= 0.92

# Predict scores using the above 20 variables:

# Create re-finalized training dataset and include all possible variables. 
rf_dat2 <- mydf2_train %>%
  select(csci, PctImp2011Cat, PctAgWs, PctImp2011Ws, PctImp2011CatRp100, PctUrbWsRp100, PctImp2011WsRp100, RdDensWs, PctAgWsRp100, PctUrbCatRp100, PctUrbWs, AgKffactWs, RdDensCatRp100, FertWs, CanalDensWs, RdDensCat, PctUrbCat, RdDensWsRp100, RdCrsWs, CBNFWs, ManureWs)

set.seed(4) # assures the data pulled is random, but sets it for the run below (makes outcome stable)
myrf2 <- randomForest(y = rf_dat2$csci, # dependent variable
  x = rf_dat2 %>%
    select(-csci),
  importance = T, 
  proximity = T, 
  ntrees = 500)  

myrf2 # examine the results. 
# 53.93% variance explained.
summary(myrf2)
plot(myrf2) # need min of 100 trees.
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

# ggsave("csci_vip_plot.png",
#      width = 25,
#      height = 10,
#      units = "cm"
#    )

# png(file="csci_vip_plot.png", units="in", width=8, height=5, res=300)
# png(file="L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_vip_plot_Mazor.png", units="in", width=8, height=5, res=300)
# vip_plot
# dev.off()

predict(myrf2) # returns out of bag predictions for training data

# Predict CSCI scores state-wide.
nottrain_prediction2 <- nottrain %>% # taking all COMIDS that haven't been used to train the model
  na.omit() %>% # remove NAs
  mutate(csci_predicted = predict(myrf2, newdata = nottrain %>% na.omit())) # using developed model (myrf2), inputting predictor variables (nottrain - which contains COMIDs and associated StreamCat data) to predict output/dependent variable (csci_predicted a.k.a. CSCI).

# rePredict CSCI scores for training and testing data (to be used for validation).
mydf2_train2 <- mydf2_train
mydf2_train2$csci_predicted <- predict(myrf2) # Adds column of predicted CSCI values to training dataset.

mydf2_test2 <- mydf2_test %>%
  mutate(csci_predicted = predict(myrf2, newdata = mydf2_test %>% select(-c(stationcode, csci, PSA6, Length_Fin)))) # Adds column of predicted CSCI values to testing dataset.

# Creates new dataset of bound rows for both ...
ca_predictions2 <- bind_rows(nottrain_prediction2 %>%
    mutate(Set = "Non-training"), # statewide COMIDs (not used for training data)
  mydf2_train2 %>%
    mutate(Set = "Training")) # COMIDS from training dataset (used for training the model).
# This creates the dataset that will be plotted.

# Create table of number of sites that fall into each category.

# # Add classification column.  3 thresholds & 4 categories
# ca_predictions2 <- ca_predictions2 %>%
#   #mutate(classification = case_when(csci_predicted < 0.61~"Very Likely Altered",   # !!! Mazor et al. 2016 1st percentile = 0.63 !!!
#   mutate(classification = case_when(round(csci_predicted, digits = 2) < 0.63~"Very Likely Altered",
#                                     round(csci_predicted, digits = 2) < 0.79~"Likely Altered",
#                                     round(csci_predicted, digits = 2) < 0.92~"Possibly Altered",
#                                     round(csci_predicted, digits = 2) >= 0.92~"Likely Unaltered")) %>%
#   mutate(class_f = factor(classification, levels = c("Very Likely Altered", "Likely Altered", "Possibly Altered", "Likely Unaltered"))) # relevel classifications

# Add classification column.  1 thresholds & 2 categories
ca_predictions2 <- ca_predictions2 %>%
  #mutate(classification = case_when(csci_predicted < 0.61~"Very Likely Altered",   # !!! Mazor et al. 2016 1st percentile = 0.63 !!!
  mutate(classification = case_when(round(csci_predicted, digits = 2) < 0.79~"Degraded",
                                    round(csci_predicted, digits = 2) >= 0.79~"Intact")) %>%
  mutate(class_f = factor(classification, levels = c("Degraded", "Intact"))) # relevel classifications

#### Results .csv ####
# Export results.
# write_csv(ca_predictions2, "csci_rf_results.csv")
# write_csv(ca_predictions2, "L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_rf_results_Mazor.csv")
# write_csv(ca_predictions2, "L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_rf_results_Mazor_IntactDegraded.csv")

# Summary table by site #.
ca_summary <- ca_predictions2 %>%
  count(class_f) # count sites statewide by classification

# Summary table by stream length (m)
ca_summary_length <- ca_predictions2 %>%
  group_by(class_f) %>% # group by classification
  summarize(length = sum(Length_Fin, na.rm=TRUE)) # sum stream lengths

# Join and export.
ca_sum <- full_join(ca_summary, ca_summary_length)
# write_csv(ca_sum, "csci_rf_results_summary.csv")
# write_csv(ca_sum, "L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_rf_results_summary_Mazor.csv")
# write_csv(ca_sum, "L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_rf_results_summary_Mazor_IntactDegraded.csv")

# Step Five - Quantile Regression model -----------------------------------
# Note - for the Healthy Watersheds Project, I did not pursue this structure.
# For additional information on how to pursue this method, see "asci_rf.R" script in this project.

# Step Six - Model validation ---------------------------------------------

# Compare predicted vs. actual results, including by PSA region.
# Adding lines of slope=1 and linear models to each plot.
val1 <- ggplot(mydf2_train2, aes(x = csci_predicted, y = csci)) +
  geom_point(color = "#2A3927", alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#2A3927") +
  labs(x = "CSCI predicted",
    y = "CSCI measured",
    title = "Training Data\nn=1,516") +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()

val1

lm1 <- lm(csci~csci_predicted, data = mydf2_train2)
summary(lm1)

val2 <- ggplot(mydf2_test2, aes(x = csci_predicted, y = csci)) +
  geom_point(color = "#3793EC", alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#3793EC") +
  scale_x_continuous(breaks = c(0.5, 0.7, 0.9)) +
  labs(x = "CSCI predicted",
    y = "CSCI measured",
    title = "Testing Data\nn=503") +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()

val2

lm2 <- lm(csci~csci_predicted, data = mydf2_test2)
summary(lm2)

# Create the full testing + training dataset to plot together.
mydf2_test2$set <- "Testing"
mydf2_train2$set <- "Training"
full_train_test <- bind_rows(mydf2_test2, mydf2_train2) %>%
  mutate(set_f = factor(set, levels = c("Training", "Testing")))

val3 <- ggplot(full_train_test, aes(x = csci_predicted, y = csci, color = set_f)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(name = "Set", values = c("#2A3927", "#3793EC"), drop = FALSE) +
  labs(x = "CSCI predicted",
    y = "CSCI measured",
    title = "All Data\nn=2,019") +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  facet_wrap(~PSA6) +
  theme_bw()

val3

val_fig <- (val1 + val2) /
  (val3)

val_fig + plot_annotation(
  title = 'CSCI Random Forest Results',
  subtitle = 'All modeling performed using StreamCAT datasets.',
  caption = 'Linear models are colored according to dataset. Lines of slope = 1 are denoted in black.'
)

# Save figure.
# ggsave("csci_rfmodel_validation.png",
#      path = "/Users/heilil/Desktop/R_figures",
#      width = 35,
#      height = 25,
#      units = "cm"
#    )

# png(file="CSCI Random Forest Results.png", units="in", width=8, height=5, res=300)
# png(file="L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/CSCI Random Forest Results_Mazor.png", units="in", width=8, height=5, res=300)
# val_fig + plot_annotation(
#   title = 'CSCI Random Forest Results',
#   subtitle = 'All modeling performed using StreamCAT datasets.',
#   caption = 'Linear models are colored according to dataset. Lines of slope = 1 are denoted in black.'
#   )
# dev.off()


lm3 <- lm(csci~csci_predicted, 
  data = full_train_test %>%
    filter(PSA6 == "Central_Valley") %>%
    filter(set_f == "Training"))

lm4 <- lm(csci~csci_predicted, 
  data = full_train_test %>%
    filter(PSA6 == "Central_Valley") %>%
    filter(set_f == "Testing"))

lm5 <- lm(csci~csci_predicted, 
  data = full_train_test %>%
    filter(PSA6 == "Chaparral") %>%
    filter(set_f == "Training"))

lm6 <- lm(csci~csci_predicted, 
  data = full_train_test %>%
    filter(PSA6 == "Chaparral") %>%
    filter(set_f == "Testing"))

lm7 <- lm(csci~csci_predicted, 
  data = full_train_test %>%
    filter(PSA6 == "Deserts_Modoc") %>%
    filter(set_f == "Training"))

lm8 <- lm(csci~csci_predicted, 
  data = full_train_test %>%
    filter(PSA6 == "Deserts_Modoc") %>%
    filter(set_f == "Testing"))

lm9 <- lm(csci~csci_predicted, 
  data = full_train_test %>%
    filter(PSA6 == "North_Coast") %>%
    filter(set_f == "Training"))

lm10 <- lm(csci~csci_predicted, 
  data = full_train_test %>%
    filter(PSA6 == "North_Coast") %>%
    filter(set_f == "Testing"))

lm11 <- lm(csci~csci_predicted, 
  data = full_train_test %>%
    filter(PSA6 == "Sierra") %>%
    filter(set_f == "Training"))

lm12 <- lm(csci~csci_predicted, 
  data = full_train_test %>%
    filter(PSA6 == "Sierra") %>%
    filter(set_f == "Testing"))

lm13 <- lm(csci~csci_predicted, 
  data = full_train_test %>%
    filter(PSA6 == "South_Coast") %>%
    filter(set_f == "Training"))

lm14 <- lm(csci~csci_predicted, 
  data = full_train_test %>%
    filter(PSA6 == "South_Coast") %>%
    filter(set_f == "Testing"))

Rsq1      <- summary(lm1)$r.squared # get r-squared value
Slp1      <- lm1$coefficients[2] # get the slope
Pval1     <- summary(lm1)$coefficients[2,4] # get p-value also anova(lm1)$'Pr(>F)'[1]
Int1      <- lm1$coefficients[1] # get the y-intercept
PInt1     <- summary(lm1)$coefficients[1,4] # get the Intercept p-value

Rsq2      <- summary(lm2)$r.squared # get r-squared value
Slp2      <- lm2$coefficients[2] # get the slope
Pval2     <- summary(lm2)$coefficients[2,4] # get p-value also anova(lm2)$'Pr(>F)'[1]
Int2      <- lm2$coefficients[1] # get the y-intercept
PInt2     <- summary(lm2)$coefficients[1,4] # get the Intercept p-value

Rsq3      <- summary(lm3)$r.squared # get r-squared value
Slp3      <- lm3$coefficients[2] # get the slope
Pval3     <- summary(lm3)$coefficients[2,4] # get p-value also anova(lm3)$'Pr(>F)'[1]
Int3      <- lm3$coefficients[1] # get the y-intercept
PInt3     <- summary(lm3)$coefficients[1,4] # get the Intercept p-value

Rsq4      <- summary(lm4)$r.squared # get r-squared value
Slp4      <- lm4$coefficients[2] # get the slope
Pval4     <- summary(lm4)$coefficients[2,4] # get p-value also anova(lm4)$'Pr(>F)'[1]
Int4      <- lm4$coefficients[1] # get the y-intercept
PInt4     <- summary(lm4)$coefficients[1,4] # get the Intercept p-value

Rsq5      <- summary(lm5)$r.squared # get r-squared value
Slp5      <- lm5$coefficients[2] # get the slope
Pval5     <- summary(lm5)$coefficients[2,4] # get p-value also anova(lm5)$'Pr(>F)'[1]
Int5      <- lm5$coefficients[1] # get the y-intercept
PInt5     <- summary(lm5)$coefficients[1,4] # get the Intercept p-value

Rsq6      <- summary(lm6)$r.squared # get r-squared value
Slp6      <- lm6$coefficients[2] # get the slope
Pval6     <- summary(lm6)$coefficients[2,4] # get p-value also anova(lm6)$'Pr(>F)'[1]
Int6      <- lm6$coefficients[1] # get the y-intercept
PInt6     <- summary(lm6)$coefficients[1,4] # get the Intercept p-value

Rsq7      <- summary(lm7)$r.squared # get r-squared value
Slp7      <- lm7$coefficients[2] # get the slope
Pval7     <- summary(lm7)$coefficients[2,4] # get p-value also anova(lm7)$'Pr(>F)'[1]
Int7      <- lm7$coefficients[1] # get the y-intercept
PInt7     <- summary(lm7)$coefficients[1,4] # get the Intercept p-value

Rsq8      <- summary(lm8)$r.squared # get r-squared value
Slp8      <- lm8$coefficients[2] # get the slope
Pval8     <- summary(lm8)$coefficients[2,4] # get p-value also anova(lm8)$'Pr(>F)'[1]
Int8      <- lm8$coefficients[1] # get the y-intercept
PInt8     <- summary(lm8)$coefficients[1,4] # get the Intercept p-value

Rsq9      <- summary(lm9)$r.squared # get r-squared value
Slp9      <- lm9$coefficients[2] # get the slope
Pval9     <- summary(lm9)$coefficients[2,4] # get p-value also anova(lm9)$'Pr(>F)'[1]
Int9      <- lm9$coefficients[1] # get the y-intercept
PInt9     <- summary(lm9)$coefficients[1,4] # get the Intercept p-value

Rsq10      <- summary(lm10)$r.squared # get r-squared value
Slp10      <- lm10$coefficients[2] # get the slope
Pval10     <- summary(lm10)$coefficients[2,4] # get p-value also anova(lm10)$'Pr(>F)'[1]
Int10      <- lm10$coefficients[1] # get the y-intercept
PInt10     <- summary(lm10)$coefficients[1,4] # get the Intercept p-value

Rsq11      <- summary(lm11)$r.squared # get r-squared value
Slp11      <- lm11$coefficients[2] # get the slope
Pval11     <- summary(lm11)$coefficients[2,4] # get p-value also anova(lm11)$'Pr(>F)'[1]
Int11      <- lm11$coefficients[1] # get the y-intercept
PInt11     <- summary(lm11)$coefficients[1,4] # get the Intercept p-value

Rsq12      <- summary(lm12)$r.squared # get r-squared value
Slp12      <- lm12$coefficients[2] # get the slope
Pval12     <- summary(lm12)$coefficients[2,4] # get p-value also anova(lm12)$'Pr(>F)'[1]
Int12      <- lm12$coefficients[1] # get the y-intercept
PInt12     <- summary(lm12)$coefficients[1,4] # get the Intercept p-value

Rsq13      <- summary(lm13)$r.squared # get r-squared value
Slp13      <- lm13$coefficients[2] # get the slope
Pval13     <- summary(lm13)$coefficients[2,4] # get p-value also anova(lm13)$'Pr(>F)'[1]
Int13      <- lm13$coefficients[1] # get the y-intercept
PInt13     <- summary(lm13)$coefficients[1,4] # get the Intercept p-value

Rsq14      <- summary(lm14)$r.squared # get r-squared value
Slp14      <- lm14$coefficients[2] # get the slope
Pval14     <- summary(lm14)$coefficients[2,4] # get p-value also anova(lm14)$'Pr(>F)'[1]
Int14      <- lm14$coefficients[1] # get the y-intercept
PInt14     <- summary(lm14)$coefficients[1,4] # get the Intercept p-value


csci_lms <- data.frame("Region" = c("Statewide", "Statewide", "Central_Valley", "Central_Valley", "Chaparral", "Chaparral", "Deserts_Modoc", "Deserts_Modoc", "North_Coast", "North_Coast", "Sierra", "Sierra", "South_Coast", "South_Coast"),
                        "Dataset" = c("Training", "Testing", "Training", "Testing", "Training", "Testing", "Training", "Testing", "Training", "Testing", "Training", "Testing", "Training", "Testing"),
                        "R2" = c(Rsq1, Rsq2, Rsq3, Rsq4, Rsq5, Rsq6, Rsq7, Rsq8, Rsq9, Rsq10, Rsq11, Rsq12, Rsq13, Rsq14),
                        "Slope" = c(Slp1, Slp2, Slp3, Slp4, Slp5, Slp6, Slp7, Slp8, Slp9, Slp10, Slp11, Slp12, Slp13, Slp14),
                        "Slope_p" = c(Pval1, Pval2, Pval3, Pval4, Pval5, Pval6, Pval7, Pval8, Pval9, Pval10, Pval11, Pval12, Pval13, Pval14),
                        "Intercept" = c(Int1, Int2, Int3, Int4, Int5, Int6, Int7, Int8, Int9, Int10, Int11, Int12, Int13, Int14),
                        "Intercept_p" = c(PInt1, PInt2, PInt3, PInt4, PInt5, PInt6, PInt7, PInt8, PInt9, PInt10, PInt11, PInt12, PInt13, PInt14))
csci_lms <- csci_lms %>%
  mutate(Slope_p = round(Slope_p, digits=6)) %>%
  mutate(Slope_p = ifelse(Slope_p < 0.0001, "<0.0001", Slope_p))
#write_csv(csci_lms,"L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci.lms.csv" )

# # Import the results of these linear models to generate summary table.
# 
# csci_lms <- read_csv(L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_lms.csv")
#csci_lms <- read_csv("L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/csci_lms.csv")

# Run the code and save table to a png file: 
summary_table <- csci_lms %>%
  gt(groupname_col = "Region", rowname_col = "Dataset") %>%
  fmt_number(columns = vars(R2, Slope, Intercept, Intercept_p), decimals = 4) %>%
  tab_header(title = "CSCI Results Validation",
    subtitle = "All modeling performed using StreamCAT datasets.") %>%
  cols_label(R2 = html("R<sup>2</sup"),
    Slope_p = html("<i>p</i>"),
    Intercept_p = html("<i>p</i>")) %>%
  cols_align(
    align = "left",
    columns = vars(R2, Slope, Slope_p, Intercept, Intercept_p))

summary_table

# Save table.
# gtsave(summary_table,
#   "csci_rfmodel_lms.png",
#   path = "/Users/heilil/Desktop/R_figures")


#webshot::install_phantomjs()
summary_table <- csci_lms %>%
  gt(groupname_col = "Region", rowname_col = "Dataset") %>%
  fmt_number(columns = vars(R2, Slope, Intercept, Intercept_p), decimals = 4) %>%
  tab_header(title = "CSCI Results Validation",
             subtitle = "All modeling performed using StreamCAT datasets.") %>%
  cols_label(R2 = html("R<sup>2</sup"),
             Slope_p = html("<i>p</i>"),
             Intercept_p = html("<i>p</i>")) %>%
  cols_align(
    align = "left",
    columns = vars(R2, Slope, Slope_p, Intercept, Intercept_p)) %>%
  gtsave(
    "csci_rfmodel_lms.png",
    path = "L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds"
    )
# 

# gtsave(file="L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_rfmodel_lms.png")

# png("L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_rfmodel_lms.png")
# summary_table
# dev.off()

process_summary <- data.frame("Dataframe" = c("ca", "ca_predictions", "ca_predictions2", "full_train_test", "mydf", "mydf2",
                                              "mydf2_test", "mydf2_test2", "mydf2_train", "mydf2_train2", "nottrain",
                                              "nottrain_prediction", "nottrain_prediction2", "ps6", "rf_dat", "rf_dat2"),
                              "Count" = c(nrow(ca), nrow(ca_predictions), nrow(ca_predictions2), nrow(full_train_test), nrow(mydf),
                                          nrow(mydf2), nrow(mydf2_test), nrow(mydf2_test2), nrow(mydf2_train), nrow(mydf2_train2),
                                          nrow(nottrain), nrow(nottrain_prediction), nrow(nottrain_prediction2), nrow(ps6),
                                          nrow(rf_dat), nrow(rf_dat2)))
#write.csv(process_summary, "L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_process_summary.csv")




# Chose not to compute confusion matrix / accuracy score since this is more applicable to categorical ouputs from random forest models -
# Instead, calculated Root Mean Squared Error (RMSE) of both training and test datasets.
# If test RMSE values are much greater than training, then possible the model has been over fit.

predtest <- predict(myrf2, mydf2_test2)
rmse(mydf2_test2$csci,predtest)
# 0.18

predtrain <- predict(myrf2, mydf2_train2)
rmse(mydf2_train2$csci,predtrain)
# 0.08

# Double checking using the original random forest dataset (rf_dat) with all 34 possible variables included to see where the error in number of predictors starts to increase dramatically (to help double check our decision to include only 20 parameters).
dc <- rfcv(rf_dat %>%
    select(-csci), 
  rf_dat$csci,
  step = 0.7, # default is 0.5
  scale="log")

dc$error.cv
#34         24         17         12          8          6          4          3          2          1 
#0.03113729 0.03154227 0.03202414 0.03273603 0.03362253 0.03452958 0.03714524 0.03682255 0.04077565 0.05935556

# Appears between 34 and 17 variables, there is an insignificant increase in error.

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

# #Map using 3 thresholds & 4 categories
# modeled_csci_map <- nhd_ca %>%
#   filter(COMID %in% mcomid) %>%
#   inner_join(ca_predictions2) %>%
#   ggplot() +
#   geom_sf(aes(color = class_f)) +
#   scale_color_manual(name = "Condition", values = c("red2", "lightpink", "lightskyblue2", "steelblue"), drop = FALSE) +
#   theme_bw()


# ggsave("csci_modeled_CA.png",
#      path = "/Users/heilil/Desktop/R_figures",
#      width = 35,
#      height = 35,
#      units = "cm"
#    )

# png(file="csci_modeled_CA.png", units="in", width=8, height=5, res=300)
# png(file="L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_modeled_CA_Mazor.png", units="in", width=8, height=5, res=300)
# modeled_csci_map
# dev.off()

#Map using 1 threshold & 2 categories
modeled_csci_map <- nhd_ca %>%
  filter(COMID %in% mcomid) %>%
  inner_join(ca_predictions2) %>%
  ggplot() +
  geom_sf(aes(color = class_f)) +
  scale_color_manual(name = "Condition", values = c("red2", "steelblue"), drop = FALSE) +
  theme_bw()
# png(file="csci_modeled_CA.png", units="in", width=8, height=5, res=300)
# png(file="L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_modeled_CA_Mazor_2.png", units="in", width=8, height=5, res=300)
# modeled_csci_map
# dev.off()


# Ventura River inset

ventura_csci_map <- nhd_vr %>%
  filter(COMID %in% mcomid) %>%
  inner_join(ca_predictions2) %>%
  ggplot() +
  geom_sf(aes(color = class_f)) +
  scale_color_manual(name = "Condition", values = c("red2", "lightpink", "lightskyblue2", "steelblue"), drop = FALSE) +
  labs(title = "Ventura River") +
  theme_bw() +
  theme(legend.position = "none")

ventura_csci_map

# png(file="csci_modeled_ventura.png", units="in", width=8, height=5, res=300)
# png(file="L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_modeled_ventura_Mazor.png", units="in", width=8, height=5, res=300)
# ventura_csci_map
# dev.off()

# ggsave("csci_modeled_Ventura.png",
#      path = "/Users/heilil/Desktop/R_figures",
#      width = 15,
#      height = 15,
#      units = "cm"
#    )

# San Juan Creek inset

sanjuan_csci_map <- nhd_sjc %>%
  filter(COMID %in% mcomid) %>%
  inner_join(ca_predictions2) %>%
  ggplot() +
  geom_sf(aes(color = class_f)) +
  scale_color_manual(name = "Condition", values = c("red2", "lightpink", "lightskyblue2", "steelblue"), drop = FALSE) +
  labs(title = "San Juan Creek") +
  theme_bw() +
  theme(legend.position = "none")

sanjuan_csci_map

# png(file="csci_modeled_sanjuan.png", units="in", width=8, height=5, res=300)
# png(file="L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_modeled_sanjuan_Mazor.png", units="in", width=8, height=5, res=300)
# sanjuan_csci_map
# dev.off()


# San Diego River inset

sandiego_csci_map <- nhd_sdr %>%
  filter(COMID %in% mcomid) %>%
  inner_join(ca_predictions2) %>%
  ggplot() +
  geom_sf(aes(color = class_f)) +
  scale_color_manual(name = "Condition", values = c("red2", "lightpink", "lightskyblue2", "steelblue"), drop = FALSE) +
  labs(title = "San Diego River") +
  theme_bw() +
  theme(legend.position = "none")

sandiego_csci_map

# png(file="csci_modeled_sandiego.png", units="in", width=8, height=5, res=300)
# png(file="L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_modeled_sandiego_Mazor.png", units="in", width=8, height=5, res=300)
# sandiego_csci_map
# dev.off()


# South coast sites inset figures

scoast <- (ventura_csci_map) /
  (sanjuan_csci_map) /
  (sandiego_csci_map)

scoast

# png(file="csci_modeled_scoast.png", units="in", width=8, height=5, res=300)
# png(file="L:/RipRAM_ES/Data/Working/healthy_watershed_random_forest/Results using published thresholds/csci_modeled_scoast_Mazor.png", units="in", width=8, height=5, res=300)
# scoast
# dev.off()

# ggsave("csci_modeled_SCoast.png",
#      path = "/Users/heilil/Desktop/R_figures",
#      width = 20,
#      height = 40,
#      units = "cm"
#    )

# Additional Notes - Healthy Watersheds project ---------------------------

#Classification options:
#"Constrained" approach, following Beck et al. 2019: Compare q10, q50, and q90 to ASCI 10th percentile threshold (i.e., 0.82)
#Likely constrained: q90 < 0.82
#Possibly constrained: q50 < 0.82
#Possibly unconstrained: q50 >= 0.82 and q10 < 0.82
#Likely unconstrained: q10 > 0.82

#"Likely condition approach: Compare q50 to three ASCI thresholds (0.67, 0.82, 0.93) @ reference sites (1st, 10th, 30th percentiles)
# Very likely altered: q50 < 0.67
# Likely altered: q50 < 0.82
# Possibly altered: q50 < 0.93
# Likely unaltered: q50 >= 0.93

# Condition approach favored by Anna per meeting on 11/3/2020.

# Works Cited:

# Hill, Ryan A., Marc H. Weber, Scott G. Leibowitz, Anthony R. Olsen, and Darren J. Thornbrugh, 2016. The Stream-Catchment (StreamCat) Dataset: A Database of Watershed Metrics for the Conterminous United States. Journal of the American Water Resources Association (JAWRA) 52:120-128. DOI: 10.1111/1752-1688.12372.
# Mazor, RD, AC Rehn, PR Ode, M Engeln, KC Schiff, ED Stein, DJ Gillett, DB Herbst, CP Hawkins. 2016. Bioassessment in complex environments: designing an index for consistent meaning in different settings

# End of R script.