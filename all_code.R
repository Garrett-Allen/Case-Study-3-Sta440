# Overall RMSE across visual fields for nonspatial model: 7.65.
# Overall RMSE across visual fields for spatial model: 14.53.
library(patchwork)
library(tidymodels)
library(spBFA)
library(yardstick)
library(ggplot2)
library(dplyr)
library(rstan)
library(loo)
#generate-data
source(file = "case3_convert_to_right_eye.R")

dat_car <- dat

#adding variable indicating visit number
dat <- dat %>% 
  group_by(patient, location, eye) %>% 
  mutate(visit = 1:n()) %>% 
  ungroup()

blind_spot <- c(26, 35) # define blind spot
W <- HFAII_Queen[-blind_spot, -blind_spot] # HFA-II visual field adjacency matrix
library(igraph)
adj.graph <- graph.adjacency(W, mode = "undirected") 
plot(adj.graph)

#filtering for early-glaucoma patients
dat <- dat %>%
  group_by(patient) %>%
  mutate(min_age = min(years)) %>%
  filter(min_age <= 45)

#creating test and train for nonspatial pointwise linear regression
train_data <- dat %>% 
  group_by(patient) %>% 
  filter(visit <= 5)

test_data <- dat %>% 
  group_by(patient) %>% 
  filter(visit > 5)


#fitting a separate model for each location
models <- train_data %>% 
  group_by(location) %>% 
  do(model = lm(dls ~ (eye + age + sex + md + time), data = .))

#collecting predictions for first model 
location_x <- filter(test_data, location == 1) %>% 
  ungroup()

model <- models$model[[1]]

location_x <- location_x %>% 
  mutate(predict_dls = predict(model, newdata = location_x))

results_sep <- location_x

#collecting predictions for all models
for(x in 2:54){
  location_x <- filter(test_data, location == x) %>% 
    ungroup()
  model <- models$model[[x]]
  location_x <- location_x %>% 
    mutate(predict_dls = predict(model, newdata = location_x))
  
  results_sep <- add_row(results_sep, location_x)
}

patients <- c("4", "12", "30", "34", "43", "72", "76", "90", "125", "128", "130")

# all_pred is for spatial model predictions
all_pred <- data.frame(patient = integer(),
                       eye = character(),
                       location = integer(),
                       visit = integer(),
                       pred = numeric())
for(patient in patients) {
  for(which_eye in c("OD", "OS")) {
    patientInt <- as.integer(patient)
    dat_filter <- dat %>%
      filter(patient==patientInt & eye==which_eye & visit > 5 & location != 26 & location != 35) %>%
      select(location, visit, dls, eye)
    
    pred <- read.csv(paste0("patient_predictions/ppd_pred_",which_eye,"_", patient, ".csv"), sep="\t")
    
    pred <- pred %>%
      mutate(location=seq(1,52))
    
    pred_long <- pred %>%
      gather(key = "visit", value = "value", -location) %>%
      mutate(visit = gsub("V", "", visit),
             visit = as.numeric(visit) + 4)
    pred_long <- pred_long[, c("location", "visit", "value")]
    colnames(pred_long) <- c("location", "visit", "pred_spatial")
    pred_long <- pred_long %>%
      filter(visit > 5 & visit <= max(dat_filter$visit))
    pred_long <- pred_long %>% 
      mutate(
        location = location + ifelse(location >= 34, 2, ifelse(location >= 26, 1, 0))
      )
    
    all_pred <- rbind(all_pred, data.frame(patient=as.integer(patient), eye=which_eye, location=as.integer(pred_long$location), visit=as.integer(pred_long$visit), pred=pred_long$pred_spatial))
  }
}

dat_filter <- dat %>%
  filter(visit > 5 & location != 26 & location != 35)

spatial_rmse <- rmse_vec(all_pred$pred, dat_filter$dls) # Final RMSE = 14.52834
#calculating rmse for nonspatial approach across visual fields

print(paste0(paste0("RMSE across visual fields for spatial: ",spatial_rmse)))
rmse_sep_models <- results_sep %>% 
  group_by(location) %>% 
  yardstick::rmse(predict_dls,dls) %>% 
  pull(.estimate) 

#removing blindspots
rmse_sep_models <- rmse_sep_models[!(seq(1:54) %in% c(26,35))] 

rmse_nonspatial <- results_sep %>% 
  filter(patient != 26 & patient != 35) %>% 
  summarize(rmse_nonspatial = rmse_vec(predict_dls, dls))

#prediction + real_dls dataframe
pred_actual <- dat_filter %>%
  ungroup() %>%
  select(patient, eye, location, visit, dls) %>%
  left_join(all_pred, by = c("patient", "eye", "location", "visit"))

rmse_spatial <- pred_actual %>%
  group_by(location) %>%
  summarize(rmse=rmse_vec(dls, pred))

rmse_spatial <- rmse_spatial$rmse

print(paste0("RMSE across visual fields for nonspatial: ",rmse_nonspatial$rmse_nonspatial))

#heatmaps of RMSE for nonspatial and spatial approaches
PlotSensitivity(Y = rmse_sep_models,
                main = "RMSE for DLS using PWLR",
                legend.lab = "RMSE for DLS", legend.round = 1)

PlotSensitivity(Y = rmse_spatial,
                main = "RMSE for DLS using Spatial Model",
                legend.lab = "RMSE for DLS", legend.round = 1)

#rmse grouped by visit 
results_visit <- results_sep %>% 
  group_by(visit) %>%
  summarize(rmse_visit = rmse_vec(predict_dls, dls))

#rmse by visit visualization
spatial_results_visit <- pred_actual %>%
  group_by(visit) %>%
  summarize(rmse_visit=rmse_vec(dls, pred))

results_visit <- results_sep %>% 
  group_by(visit) %>%
  summarize(rmse_visit = rmse_vec(predict_dls, dls))

results_visit$source <- "non_spatial"
spatial_results_visit$source <- "spatial"
combined_visit <- rbind(results_visit, spatial_results_visit)
ggplot(data = combined_visit, aes(x = visit, y = rmse_visit, color = source)) +
  geom_line(size = 1.2) +
  geom_point(size = 3, shape = 21, fill = "#FFFFFF") +
  labs(
    title = "RMSE by Visit",
    x = "Visit",
    y = "RMSE"
  ) +
  scale_color_manual(values = c("non_spatial" = "#0072B2", "spatial" = "red")) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "#D3D3D3"),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

patients <- c("4", "12", "30", "34", "43", "72", "76", "90", "125", "128", "130")

#calculating spatial rmse for each patient
spatial_rmse_patient <- data.frame(patient = integer(),
                                   eye = character(),
                                   rmse_patient = numeric(),
                                   stringsAsFactors = FALSE)
for(patient in patients) {
  for(which_eye in c("OD", "OS")) {
    patientInt <- as.integer(patient)
    dat_filter <- dat %>%
      filter(patient==patientInt & eye==which_eye & visit > 5 & location != 26 & location != 35) %>%
      select(location, visit, dls, eye)
    
    pred <- read.csv(paste0("patient_predictions/ppd_pred_",which_eye,"_", patient, ".csv"), sep="\t")
    
    pred <- pred %>%
      mutate(location=seq(1,52))
    
    pred_long <- pred %>%
      gather(key = "visit", value = "value", -location) %>%
      mutate(visit = gsub("V", "", visit),
             visit = as.numeric(visit) + 4)
    pred_long <- pred_long[, c("location", "visit", "value")]
    colnames(pred_long) <- c("location", "visit", "pred_spatial")
    pred_long <- pred_long %>%
      filter(visit > 5 & visit <= max(dat_filter$visit))
    pred_long <- pred_long %>% 
      mutate(
        location = location + ifelse(location >= 34, 2, ifelse(location >= 26, 1, 0))
      )
    
    
    rmse_val <- rmse_vec(pred_long$pred_spatial, dat_filter$dls)
    
    spatial_rmse_patient <- rbind(spatial_rmse_patient, data.frame(patient=as.integer(patient), eye=which_eye, rmse_patient=rmse_val, stringsAsFactors = FALSE))
  }
}

#rmse summarized by patient/eye
results_patient <- results_sep %>%
  group_by(patient, eye) %>%
  summarize(rmse_patient=rmse_vec(predict_dls, dls))

# bar plot for rmse by patient/eye
results_patient$source <- "rmse_patient"
spatial_rmse_patient$source <- "spatial_rmse_patient"
combined_data <- rbind(results_patient, spatial_rmse_patient)
ggplot(combined_data, aes(x = interaction(patient, eye), y = rmse_patient, fill = source)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, color = "black") +
  labs(title = "RMSE by Patient-Eye Combination", x = "Patient-Eye", y = "RMSE") +
  scale_fill_manual(values = c("#0072B2", "red")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 18, face = "bold", margin = margin(b = 10)),
        axis.title = element_text(size = 14, margin = margin(t = 10)),
        axis.text = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom")

#see QMD for how we fit models/extracted results, otherwise .R file would take
#very long to run. 

early_onset_glaucoma <- unique(train_data$patient)

dir.create('fitted_models')
#setting options
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
blind_spot <- c(26, 35) # define blind spot
dat_car <- dat_car[!dat_car$location %in% blind_spot, ] # remove blind spot locations
dat_car$LocationNew <- as.numeric(as.factor(dat_car$location)) # reorder locations to be ordered integers
dat_car <- dat_car %>%
  group_by(patient) %>%
  mutate(visit = match(time, unique(time)))

#fitting OD eyes for all early onset patients
for(x in early_onset_glaucoma){
  data = dat_car %>%
    filter(patient == x , eye == "OD")
  
  time <- unique(data$time)
  T <- length(unique(data$visit))
  W <- HFAII_Queen[-blind_spot, -blind_spot] # visual field adjacency matrix
  T_train <- 3
  dataTrain <- data[data$visit <= T_train, ]
  n <- length(unique(dataTrain$location))
  stan_data <- list(
    y = dataTrain$dls,
    x = dataTrain$time,
    s = dataTrain$LocationNew,
    N = n * T_train,
    n = n,
    T = n,
    W = W
  )
  iter = 100
  n_iter = iter/100 
  model_compiled <- stan_model("spatial.stan")
  fit <- sampling(model_compiled, data = stan_data, chains = 1, iter = iter,
                  refresh = n_iter)
  
  saveRDS(fit, file = paste0("fitted_models/fit_",x,"_od.rds")) # an option for saving the stan model fit
}


#fitting all OS eyes
for(x in early_onset_glaucoma){
  data = dat_car %>%
    filter(patient == x , eye == "OS")
  
  time <- unique(data$time)
  T <- length(unique(data$visit))
  W <- HFAII_Queen[-blind_spot, -blind_spot] # visual field adjacency matrix
  T_train <- 3
  dataTrain <- data[data$visit <= T_train, ]
  n <- length(unique(dataTrain$location))
  stan_data <- list(
    y = dataTrain$dls,
    x = dataTrain$time,
    s = dataTrain$LocationNew,
    N = n * T_train,
    n = n,
    T = n,
    W = W
  )
  iter = 4000
  n_iter = iter/100 
  model_compiled <- stan_model("spatial.stan")
  fit <- sampling(model_compiled, data = stan_data, chains = 1, iter = iter,
                  refresh = n_iter)
  
  saveRDS(fit, file = paste0("fitted_models/fit_",x,"_os.rds")) # an option for saving the stan model fit
}





# sensitivity analysis

# OD: Beta Distribution from 1 to 5, saves file and runs it
data = dat_car %>%
  filter(patient == 4 , eye == "OD")

time <- unique(data$time)
T <- length(unique(data$visit))
W <- HFAII_Queen[-blind_spot, -blind_spot] # visual field adjacency matrix
T_train <- 3
dataTrain <- data[data$visit <= T_train, ]
n <- length(unique(dataTrain$location))
stan_data <- list(
  y = dataTrain$dls,
  x = dataTrain$time,
  s = dataTrain$LocationNew,
  N = n * T_train,
  n = n,
  T = n,
  W = W
)
iter = 4000
n_iter = iter/100 
model_compiled <- stan_model("15spatial.stan")
fit <- sampling(model_compiled, data = stan_data, chains = 1, iter = iter,
                refresh = n_iter)

saveRDS(fit, file = "fit_4_right_sensitivity.rds")
fit_4_right_sensitivity <- readRDS("~/case3-team06/fit_4_right_sensitivity.rds")
# OS: Beta Distribution from 1 to 5, saves file and runs it
data = dat_car %>%
  filter(patient == 4 , eye == "OS")

time <- unique(data$time)
T <- length(unique(data$visit))
W <- HFAII_Queen[-blind_spot, -blind_spot] # visual field adjacency matrix
T_train <- 3
dataTrain <- data[data$visit <= T_train, ]
n <- length(unique(dataTrain$location))
stan_data <- list(
  y = dataTrain$dls,
  x = dataTrain$time,
  s = dataTrain$LocationNew,
  N = n * T_train,
  n = n,
  T = n,
  W = W
)
iter = 4000
n_iter = iter/100 
model_compiled <- stan_model("15spatial.stan")
fit <- sampling(model_compiled, data = stan_data, chains = 1, iter = iter,
                refresh = n_iter)

saveRDS(fit, file = "fit_4_left_sensitivity.rds")
fit_4_left_sensitivity <- readRDS("~/case3-team06/fit_4_left_sensitivity.rds")

#OD: Beta Distribution from 5 to 1, saves file and runs it
data = dat_car %>%
  filter(patient == 4 , eye == "OD")

time <- unique(data$time)
T <- length(unique(data$visit))
W <- HFAII_Queen[-blind_spot, -blind_spot] # visual field adjacency matrix
T_train <- 3
dataTrain <- data[data$visit <= T_train, ]
n <- length(unique(dataTrain$location))
stan_data <- list(
  y = dataTrain$dls,
  x = dataTrain$time,
  s = dataTrain$LocationNew,
  N = n * T_train,
  n = n,
  T = n,
  W = W
)
iter = 4000
n_iter = iter/100 
model_compiled <- stan_model("51spatial.stan")
fit <- sampling(model_compiled, data = stan_data, chains = 1, iter = iter,
                refresh = n_iter)

saveRDS(fit, file = "fit_4_right_sensitivity51.rds")
fit_4_right_sensitivity51 <- readRDS("~/case3-team06/fit_4_right_sensitivity51.rds")

#OS: Beta Distribution from 5 to 1, saves file and runs it
data = dat_car %>%
  filter(patient == 4 , eye == "OS")

time <- unique(data$time)
T <- length(unique(data$visit))
W <- HFAII_Queen[-blind_spot, -blind_spot] # visual field adjacency matrix
T_train <- 3
dataTrain <- data[data$visit <= T_train, ]
n <- length(unique(dataTrain$location))
stan_data <- list(
  y = dataTrain$dls,
  x = dataTrain$time,
  s = dataTrain$LocationNew,
  N = n * T_train,
  n = n,
  T = n,
  W = W
)
iter = 4000
n_iter = iter/100 
model_compiled <- stan_model("51spatial.stan")
fit <- sampling(model_compiled, data = stan_data, chains = 1, iter = iter,
                refresh = n_iter)

saveRDS(fit, file = "fit_4_left_sensitivity51.rds")
fit_4_left_sensitivity51 <- readRDS("~/case3-team06/fit_4_left_sensitivity51.rds")


