####### The following script contains all the REVISED code used for the
## statistical analyses and data visualization for the submitted manuscript:


## Drivers of mosquito indoor free-flight and resting

#####################################################

### load packages
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(knitr)
library(visreg)
library(sjPlot)
library(effects)
library(lme4)
library(Matrix)
library(xtable)
library(texreg)
library(visreg)
library(bayesplot)
library(broom.mixed)
library(merTools)
library(DHARMa)
library(influence.ME)
library(glmmTMB)
library(MuMIn)
library(ggeffects)
library(RColorBrewer)
library(emmeans)
library(performance)
library(FSA)
library(purrr)
library(dunn.test)
library(ggpattern)
library(rlang)
library(broom)
library(viridis)
library(reshape2)
library(mgcv)
library(plotly)
library(ggsignif)

#####################################################
#####################################################

#### Statistical summaries of collected data 

## PFMD Flight Trajectories

# read in the PFMD trajectory CSV file
tracks <- read_csv("track_hobo_merged_CLEAN_8_15_2026.csv")

# set categorical variables as factors ahead of analysis
tracks$id <- as.factor(tracks$id)
tracks$strain <- as.factor(tracks$strain)
tracks$sex <- as.factor(tracks$sex)
tracks$color <- as.factor(tracks$color)
tracks$tod <- as.factor(tracks$tod)
tracks$unique.id <- as.factor(tracks$unique.id)

# recode FIELD for WT
tracks <- tracks %>% 
  mutate(strain = if_else(strain == "FIELD", "WT", strain))

# set strain as factor again
tracks$strain <- as.factor(tracks$strain)


# summarizing total number of trajectories
# stratifying by strain, sex, and color (diel replication collapses to n_sessions)
trajectory_summary <- tracks %>% 
  group_by(strain, sex, color) %>% 
  summarise(n_sessions = n_distinct(unique.id),
            n_tracks = n_distinct(id),
            .groups = "drop")


# Including diel replication as a grouping factor (SI table)
trajectory_diel_summary <- tracks %>% 
  group_by(strain, sex, color, tod) %>% 
  summarise(n_tracks = n_distinct(id),
            .groups = "drop")


# Obtain average number of tracks for each strain x sex x color combination
# Average across three replicates
avg_tracks_per_combi <- tracks %>% 
  group_by(strain, sex, color, unique.id) %>% 
  summarize(n_tracks = n_distinct(id), .groups = "drop") %>%
  group_by(strain, sex, color) %>% 
  summarize(mean_tracks = mean(n_tracks),
            sd_tracks = sd(n_tracks),
            n_sessions = n(),
            .groups = "drop")



## Sticky Trap Landing Counts
# first read in landing data
sticky <- read_csv("Processed_Data/sticky_hobo_model_ready_8_19_2024.csv")

# set appropriate variables as factor
sticky$group <- as.factor(sticky$group) # sex and physiological status
sticky$color <- as.factor(sticky$color) # color
sticky$replicate <- as.factor(sticky$replicate) # diel replicate / Time of Day
sticky$replicate.id <- as.factor(sticky$replicate.id) # unique replicate id

# recode MID for WT, and ROCK for Lab
sticky <- sticky %>% 
  mutate(strain = if_else(strain == "MID", "WT", strain))

sticky <- sticky %>% 
  mutate(strain = if_else(strain == "ROCK", "LAB", strain))

# set strain as factor again
sticky$strain <- as.factor(sticky$strain)

# recode group column to sex column
sticky <- sticky %>% 
  rename(sex = group) %>% 
  mutate(sex = recode(sex, "BF" = "Blood-Fed",
                      "UF" = "Unfed",
                      "M" = "Male"))

# set sex as factor
sticky$sex <- as.factor(sticky$sex)

# summarizing total number of landing counts
# stratifying by strain, sex, and color (diel replication collapses to n_sessions)\
landing_count_summary <- sticky %>% 
  group_by(strain, sex, color) %>% 
  summarize(n_sessions = n_distinct(replicate.id),
            n_landings = n(), 
            .groups = "drop")

# Including diel replication as a grouping factor (SI Table)
landing_count_diel_summary <- sticky %>% 
  group_by(strain, sex, color, replicate) %>% 
  summarize(n_landings = n(), 
            .groups = "drop")

# Obtain average number of tracks for each strain x sex x color combination
# Average across three replicates
avg_landings_per_combi <- sticky %>% 
  group_by(strain, sex, color, replicate.id) %>% 
  summarize(n_landings = n(), .groups = "drop") %>%
  group_by(strain, sex, color) %>% 
  summarize(mean_landings = mean(n_landings),
            sd_landings = sd(n_landings),
            n_sessions = n(),
            .groups = "drop")


## Combining landings and trajectories into one table
# merge these initial summaries into a singular table
main_summary_table <- full_join(trajectory_summary, landing_count_summary,
                                by = c("strain", "sex", "color"))

# rename columns
colnames(main_summary_table) <- c("Strain", "Sex", "Color", "PFMD_Sessions", 
                                  "N_Tracks", "Landing_Sessions", "N_Landings")

# print table
print(main_summary_table, n = Inf)


# Getting subtotals (by Strain, Sex, and Color) and grand total

# Strain-level subtotals across all replicates
strain_subtotals <- main_summary_table %>% 
  group_by(Strain) %>% 
  summarize(Sex = "Subtotal",
            Color = "",
            PFMD_Sessions = sum(PFMD_Sessions),
            N_Tracks = sum(N_Tracks),
            Landing_Sessions = sum(Landing_Sessions),
            N_Landings = sum(N_Landings),
            .groups = "drop")

# Sex-level subtotals across all replicates
sex_subtotals <- main_summary_table %>% 
  group_by(Strain, Sex) %>% 
  summarize(Color = "Subtotal",
            PFMD_Sessions = sum(PFMD_Sessions),
            N_Tracks = sum(N_Tracks),
            Landing_Sessions = sum(Landing_Sessions),
            N_Landings = sum(N_Landings),
            .groups = "drop")

# Color-level subtotals across all replicates (collapsing sex)
color_subtotals <- main_summary_table %>% 
  group_by(Strain, Color) %>% 
  summarize(PFMD_Sessions = sum(PFMD_Sessions),
            N_Tracks = sum(N_Tracks),
            Landing_Sessions = sum(Landing_Sessions),
            N_Landings = sum(N_Landings),
            .groups = "drop")

# Grand total across all replicates
grand_total <- main_summary_table %>% 
  summarize(Strain = "Total",
            Sex = "",
            Color = "",
            PFMD_Sessions = sum(PFMD_Sessions),
            N_Tracks = sum(N_Tracks),
            Landing_Sessions = sum(Landing_Sessions),
            N_Landings = sum(N_Landings))

# Interleaving all subtotals to grand total
# make sure columns are consistent all across
col_order <- c("Strain", "Sex", "Color", "PFMD_Sessions",
               "N_Tracks", "Landing_Sessions", "N_Landings")


# assemble summaries
totals_assembled <- main_summary_table %>% 
  group_split(Strain, Sex) %>% 
  map_dfr(function(chunk) {
    st <- chunk$Strain[1] ; sx <- chunk$Sex[1]
    sub <- sex_subtotals %>%  filter(Strain == st, Sex == sx)
    bind_rows(chunk, sub)
  })

# insert subtotals after ea Strain, and grand total at the end
final_totals_table <- totals_assembled %>% 
  group_split(Strain) %>% 
  map_dfr(function(block) {
    st <- block$Strain[1]
    bind_rows(block, strain_subtotals %>% filter(Strain == st))
  }) %>% 
  bind_rows(grand_total) %>% 
  dplyr::select(all_of(col_order))



## HOBO microclimate data (Temperature, RH, VPD)

# read in the data
hobo <- read_csv("Processed_Data/hobo_all_data.csv")

# recode each cobo's height with corresponding height value (cm)
hobo <- hobo %>% 
  mutate(HOBO = case_when(HOBO == "One" ~ 192,
         HOBO == "Four" ~ 128,
         HOBO == "Five" ~ 64,
         HOBO == "Three" ~ 0))

# set HOBO as a factor
hobo$HOBO <- as.factor(hobo$HOBO)

# calculate VPD from Temp and RH (Tetens Equation)
hobo <- hobo %>% 
  mutate(es = 0.6108 * exp(17.27 * Temp / (Temp + 237.3)),
         ea = es * (RH / 100),
         vpd = es - ea)

# Now we can summarize all variables (Tempp, RH, and VPD)
hobo_summary_height <- hobo %>% 
  group_by(HOBO) %>% 
  summarize(n = n(),
            mean_temp = mean(Temp), sd_temp = sd(Temp), 
            min_temp = min(Temp), max_temp = max(Temp),
            mean_rh = mean(RH), sd_rh = sd(RH), 
            min_rh = min(RH), max_rh = max(RH),
            mean_vpd = mean(vpd), sd_vpd = sd(vpd), 
            min_vpd = min(vpd), max_vpd = max(vpd),
            .groups = "drop")


#####################################################
#####################################################


### Correlated Random Walk Model Results and Impact of Strain on Mosquito Height Preferences

## Visual representation of 25 tracks picked at random from the entire PFMD dataset (Panel A, Figure 1)
# Quick per-track step lengths on full PFMD dataset
tracks_sl <- tracks %>%
  arrange(id, datetime) %>%
  group_by(id) %>%
  mutate(
    sl_quick = sqrt((x_cm - lag(x_cm))^2 + (y_cm - lag(y_cm))^2 + (z_cm - lag(z_cm))^2)
  ) %>%
  ungroup()

# random 10 clean tracks (still exclude jump/stitched tracks)
clean_ids <- tracks_sl %>%
  group_by(id) %>%
  summarise(has_jump = any(sl_quick > 3, na.rm = TRUE), .groups = "drop") %>%
  filter(!has_jump) %>%
  pull(id)

set.seed(123)
sample_ids <- sample(clean_ids, 10)

top10_obs <- tracks_sl %>% filter(id %in% sample_ids) %>%
  mutate(
    x_centered = x_cm - 151)

# extract start and end points of each track
starts <- top10_obs %>% group_by(id) %>% slice_head(n = 1) %>% ungroup()
ends   <- top10_obs %>% group_by(id) %>% slice_tail(n = 1) %>% ungroup()

# plotting Panel A
plot_ly() %>%
  # the track lines
  add_trace(data = top10_obs, x = ~x_centered, y = ~z_cm, z = ~y_cm,
    color = ~id, colors = "plasma", type = "scatter3d", mode = "lines",
    line = list(width = 7), showlegend = FALSE) %>%
  # start points — circles
  add_trace(data = starts, x = ~x_centered, y = ~z_cm, z = ~y_cm, 
            type = "scatter3d", mode = "markers", 
            marker = list(symbol = "circle",size = 6, color = "green"),
            name = "Start", showlegend = TRUE) %>%
  # end points — triangles (diamond used, see note)
  add_trace(data = ends, x = ~x_centered, y = ~z_cm, z = ~y_cm,
            type = "scatter3d", mode = "markers", 
            marker = list(symbol = "diamond", size = 6, color = "red"),
            name = "End", showlegend = TRUE) %>%
  layout(font = list(size = 45), legend = list(font = list(size = 45)),
         scene = list(aspectmode = "cube",
                      xaxis = list(title = "X (cm)", backgroundcolor = "lightgrey",
                                   gridcolor = "black", gridwidth = 3, 
                                   showbackground = TRUE, color = "black", 
                                   titlefont = list(size = 45), 
                                   tickfont = list(size = 15)),
                      yaxis = list(title = "Z (cm)", backgroundcolor = "lightgrey",
                                   gridcolor = "black", gridwidth = 3, 
                                   showbackground = TRUE, color = "black",
                                   titlefont = list(size = 45),
                                   tickfont = list(size = 15)),
                      zaxis = list(title = "Height (cm)", range = c(0, 204.5992),
                                   backgroundcolor = "lightgrey", gridcolor = "black",
                                   gridwidth = 3, showbackground = TRUE,
                                   color = "black", titlefont = list(size = 45),
                                   tickfont = list(size = 15)))) %>%
  config(toImageButtonOptions = list(format = "png", width = 2400, height = 2000))



## Obtaining distributions of Step Length and Turning Angles from LAB, UF, W subset
lab_uf_w_tracks <- subset(tracks, strain == "LAB" & sex == "Unfed" & color == "W")
lab_uf_w_tracks <- droplevels(lab_uf_w_tracks)

# arrange data by ID and datetime
lab_uf_w_tracks <- lab_uf_w_tracks %>% 
  arrange(id, datetime)

## Source functions for processing Correalted Random Walk Model workflow 
source("/Volumes/External/DATA/AIM_1/PFMD/Track_Projection_Code/SourceFile_3D_4_Troubleshooting.R")

### Now let's produce the empirical distributions of trajectory step lengths and turning angles
# Read in movement data 
data.lab <- lab_uf_w_tracks


# Initialize dx, dy, dz columns
data.lab <- as.data.frame(data.lab)
data.lab <- data.lab[,names(data.lab) %in% c("id", "datetime", "x_cm", "y_cm", "z_cm")]
head(data.lab)


# Step 1. Get dx, dy, dz and step-length columns using pyt_3d function (3D pythagorean)
##### GENERATE STEPLENGTHS AND DX, DY, DZ
data.lab2 <- pyt_3d(data.lab)
# Now you have dx, dy, dz columns and a 3D step-length column
head(data.lab2)


# Step 2. Calculate absolute and relative turning angles in both the XY and XZ planes
###### GET TURNING ANGLE DISTRIBUTIONS IN XY XZ
data.lab3 <- get_turning_angles(data.lab2)
head(data.lab3)


# Step 3. Some SLs and TAs are quite extreme and will be excluded
data.lab4 <- data.lab3 %>%
  arrange(id, datetime) %>%
  group_by(id) %>%
  mutate(
    jump = !is.na(sl) & sl > 3,
    sl = ifelse(jump, NA, sl),
    rel.angle.xy = ifelse(jump | lag(jump, default = FALSE), NA, rel.angle.xy),
    rel.angle.xz = ifelse(jump | lag(jump, default = FALSE), NA, rel.angle.xz)
  ) %>%
  ungroup()

# Step 4. Fit (population-level) empirical step-lengths to theoretical distribution 
# Need to enter your data set (your_df = ), your distribution (either dist = "gamma", or dist = "weibull", and whether or not you want to plot them (plot = T or plot = F))
# Output is a list: 1st slot is the fit itself, 2nd slot is the est mean, 3rd slot is the est SD
# A Warning message about NaNs being produced occurs here sometimes, unclear as to why yet

## The plot produced here is Panel B in Figure 1
save_sl_fit <- fit_steplength(your_df = data.lab4, dist = "gamma", plot = T)
save_sl_fit

# getting SL plot - Panel B in Figure 1
sl_df <- data.frame(sl = na.omit(data.lab4$sl))

shape_sl <- save_sl_fit[[2]] # gamma shape
rate_sl <- save_sl_fit[[3]] # gamma rate

xx_sl <- seq(0, max(sl_df$sl), length.out = 400)
gam_df <- data.frame(x = xx_sl,
                     d = dgamma(xx_sl, shape = shape_sl, rate = rate_sl))

ggplot() +
  geom_histogram(data = sl_df, aes(x = sl, y = after_stat(density)),
                 bins = 60, fill = "#66C2A5", color = "black", linewidth = 1) +
  geom_line(data = gam_df, aes(x = x, y = d), color = "red", linewidth = 1.2, alpha = 0.7) +
  annotate("text", x = 3, y = Inf,
           label = sprintf("Gamma\nShape (\u03B1) = %.2f\nRate (\u03B2) = %.2f", shape_sl, rate_sl),
           hjust = 1, vjust = 1.2, size = 10, color = "red") +
  scale_x_continuous(breaks = seq(0, ceiling(max(sl_df$sl)), by = 0.5)) +
  labs(title = "Step Lengths", x = "Step Length (cm)", y = "Density") +
  theme_bw(base_size = 45) +
  theme(plot.title = element_text(hjust = 0.5))

# Step 4. Fit (population-level) empirical relative step-lengths in both XY and XZ planes to theoretical distributions
# Need to enter your data set (your_df = )
# You will get a warning, about coercion to class circular, not an issue 
# Output is a list: 1st slot is the fit for XY, 2nd slot is the fit for XZ

save_ta_fit <- fit_turning_angles(your_df = data.lab4)
save_ta_fit

# getting TA plots - these are Panel C in Figure 1
library(circular)

# empirical angles from the data
ang_xy <- na.omit(data.lab4$rel.angle.xy)
ang_xz <- na.omit(data.lab4$rel.angle.xz)
# wrap to [-pi, pi] for display
ang_xy[ang_xy > pi] <- ang_xy[ang_xy > pi] - 2*pi
ang_xz[ang_xz > pi] <- ang_xz[ang_xz > pi] - 2*pi

# fitted parameters from save_ta_fit
mu_xy    <- as.numeric(save_ta_fit$xy_distribution$mu)
kappa_xy <- save_ta_fit$xy_distribution$kappa
mu_xz    <- as.numeric(save_ta_fit$xz_distribution$mu)
kappa_xz <- save_ta_fit$xz_distribution$kappa

par(mfrow = c(1, 2))

## --- XY ---
hist(ang_xy, prob = TRUE, breaks = 60, col = "lightblue", border = "white",
     main = "XY Turning Angles", xlab = "Relative Turning Angle (rad)", xlim = c(-pi, pi))
xx <- seq(-pi, pi, length.out = 400)
dd_xy <- exp(kappa_xy * cos(xx - mu_xy)) / (2 * pi * besselI(kappa_xy, 0))   # von Mises density
lines(xx, dd_xy, col = "red", lwd = 2.5)
legend("topright", sprintf("von Mises\nκ = %.2f", kappa_xy), col = "red", lwd = 2, bty = "n")

# ggplot equivalent - Figure 1 Panel C-a
# data frames for the histogram and the von Mises curve
ang_df_xy <- data.frame(angle = ang_xy)
xx <- seq(-pi, pi, length.out = 400)
vm_df_xy <- data.frame(
  x = xx,
  d = exp(kappa_xy * cos(xx - mu_xy)) / (2 * pi * besselI(kappa_xy, 0))
)

ggplot() +
  geom_histogram(data = ang_df_xy, aes(x = angle, y = after_stat(density)),
                 bins = 60, fill = "#9E7BB5", color = "black", linewidth = 1) +
  geom_line(data = vm_df_xy, aes(x = x, y = d), color = "red", linewidth = 1.2, alpha = 0.7) +
  annotate("text", x = pi, y = Inf,
           label = sprintf("von Mises\nκ = %.2f", kappa_xy),
           hjust = 1, vjust = 1.2, size = 12, color = "red") +
  scale_x_continuous(limits = c(-pi, pi)) +
  labs(title = "XY Turning Angles",
       x = "Relative Turning Angle (rad)", y = "Density") +
  theme_bw(base_size = 45) +
  theme(plot.title = element_text(hjust = 0.5))



# --- XZ ---
hist(ang_xz, prob = TRUE, breaks = 60, col = "lightgreen", border = "white",
     main = "XZ Turning Angles", xlab = "Relative Turning Angle (rad)", xlim = c(-pi, pi))
dd_xz <- exp(kappa_xz * cos(xx - mu_xz)) / (2 * pi * besselI(kappa_xz, 0))
lines(xx, dd_xz, col = "red", lwd = 2.5)
legend("topright", sprintf("von Mises\nκ = %.2f", kappa_xz), col = "red", lwd = 2, bty = "n")

par(mfrow = c(1, 1))

# ggplot equivalent - Figure 1 Panel C-b
# data frames for the histogram and the von Mises curve
ang_df_xz <- data.frame(angle = ang_xz)
xx <- seq(-pi, pi, length.out = 400)
vm_df_xz <- data.frame(
  x = xx,
  d = exp(kappa_xz * cos(xx - mu_xz)) / (2 * pi * besselI(kappa_xz, 0))
)

ggplot() +
  geom_histogram(data = ang_df_xz, aes(x = angle, y = after_stat(density)),
                 bins = 60, fill = "#C39BD3", color = "black", linewidth = 1) +
  geom_line(data = vm_df_xz, aes(x = x, y = d), color = "red", linewidth = 1.2, alpha = 0.7) +
  annotate("text", x = pi, y = Inf,
           label = sprintf("von Mises\nκ = %.2f", kappa_xz),
           hjust = 1, vjust = 1.2, size = 12, color = "red") +
  scale_x_continuous(limits = c(-pi, pi)) +
  labs(title = "XZ Turning Angles",
       x = "Relative Turning Angle (rad)", y = "Density") +
  theme_bw(base_size = 45) +
  theme(plot.title = element_text(hjust = 0.5))

# Step 5. Simulate trajectory from empirically fitted theoretical distributions
# define bounds: Height uses the physical tracking volume (204.5992 - max observed value across full PFMD dataset)
# Width includes entire 3m of hut, and Depth (z) only the final 4 to 5m of the hut (most complete active region subset)
my_bounds <- list(x = c(0, 302),
                  y = c(0, 204.5992),
                  z = c(400, 501))

# set number of trajectories to be simulated
num_trajectories <- 1000

# Initialize an empty list to store each simulated trajectory
simulated_trajectories <- list()

# set seed then run loop to construct the trajectories
set.seed(123)

for (i in 1:num_trajectories) {
  # Simulate trajectories
  sim_df <- simulate_3d(your_df = data.lab4,
                        num_points = 365,
                        sl_dist = save_sl_fit,
                        ta_dist = save_ta_fit,
                        distribution = "gamma",
                        bounds = my_bounds)
  
  # Remove rows where 'y' column is NA
  sim_df <- sim_df %>%
    filter(!is.na(y))
  
  # Add a unique trajectory ID
  sim_df$trajectory_id <- i
  
  # Store in the list
  simulated_trajectories[[i]] <- sim_df
}

# Combine all trajectories into one data frame
all_trajectories <- do.call(rbind, simulated_trajectories)

# set trajectory id as factor
all_trajectories$trajectory_id <- as.factor(all_trajectories$trajectory_id)

# now lets randomly pick 10 simulated tracks 
set.seed(123)

sim_ids <- sample(unique(all_trajectories$trajectory_id), 10)

top10_sim <- all_trajectories %>% 
  filter(trajectory_id %in% sim_ids) %>%
  mutate(x_centered = x - 151)

# extract start and end points of each track
starts_sim <- top10_sim %>% group_by(trajectory_id) %>% slice_head(n = 1) %>% ungroup()
ends_sims   <- top10_sim %>% group_by(trajectory_id) %>% slice_tail(n = 1) %>% ungroup()

# plotting - Panel D Figure 1
plot_ly() %>%
  # the track lines
  add_trace(
    data = top10_sim,
    x = ~x_centered, y = ~z, z = ~y,
    color = ~trajectory_id, colors = "plasma",
    type = "scatter3d", mode = "lines",
    line = list(width = 7),
    showlegend = FALSE
  ) %>%
  # start points — circles
  add_trace(
    data = starts_sim,
    x = ~x_centered, y = ~z, z = ~y,
    type = "scatter3d", mode = "markers",
    marker = list(symbol = "circle", size = 6, color = "green"),
    name = "Start", showlegend = TRUE
  ) %>%
  # end points — triangles (diamond used, see note)
  add_trace(
    data = ends_sims,
    x = ~x_centered, y = ~z, z = ~y,
    type = "scatter3d", mode = "markers",
    marker = list(symbol = "diamond", size = 6, color = "red"),
    name = "End", showlegend = TRUE
  ) %>%
  layout(
    font = list(size = 45),
    legend = list(font = list(size = 45)),
    scene = list(
      aspectmode = "cube",
      xaxis = list(title = "X (cm)", range = c(-150, 150), backgroundcolor = "lightgrey",
                   gridcolor = "black", gridwidth = 3, showbackground = TRUE, 
                   color = "black", titlefont = list(size = 45),
                   tickfont = list(size = 15)),
      yaxis = list(title = "Z (cm)", backgroundcolor = "lightgrey",
                   gridcolor = "black", gridwidth = 3, showbackground = TRUE, 
                   color = "black", titlefont = list(size = 45),
                   tickfont = list(size = 15)),
      zaxis = list(title = "Height (cm)", range = c(0, 204.5992),
                   backgroundcolor = "lightgrey", gridcolor = "black",
                   gridwidth = 3, showbackground = TRUE, color = "black",
                   titlefont = list(size = 45),
                   tickfont = list(size = 15))
    )
  )


## Panel E Figure 1 - y coordinate density curves
# This panel will show the density curves of LAB vs WT UF mosquitoes under W
# color, compared visually against the the density curve of bootsrapped SIM data


# Bootstrapping SIM height density for a 95% CI band
sim_by_trajectory <- split(all_trajectories$y, all_trajectories$trajectory_id)
grid <- seq(0, 204.5992, length.out = 200)
n_boot <- 500

sim_trajectories_boot <- replicate(n_boot, {
  samp <- sample(names(sim_by_trajectory), replace = TRUE)
  vals <- unlist(sim_by_trajectory[samp])
  density(vals, from = 0, to = 204.5992, n = 200)$y
})

sim_trajectory_CI <- data.frame(height = grid, 
                                median = apply(sim_trajectories_boot, 1, median),
                                lower = apply(sim_trajectories_boot, 1, quantile, 0.025),
                                upper = apply(sim_trajectories_boot, 1, quantile, 0.975))


# Obtain observed flight heights for LAB and WT
# LAB
lab_flight <- lab_uf_w_tracks

# WT 
wt_flight <- subset(tracks, strain == "WT" & sex == "Unfed" & color == "W")
wt_flight <- droplevels(wt_flight)

# bind LAB and WT
flight_densitities <- rbind(lab_flight, wt_flight)

  
# With flight done, now its landing density curves from sticky data
landing_densities <- sticky %>%
  filter(color == "W" & sex == "Unfed")

landing_densities <- droplevels(landing_densities)

# un-normalize landing height values
landing_densities <- landing_densities %>%
  mutate(resting.height = resting.height * 192)

# combine flight and landing into one df with a behavior column
flight_densitities <- flight_densitities %>%
  transmute(height = y_cm, strain = strain, behavior = "FLIGHT")

landing_densities <- landing_densities %>%
  transmute(height = resting.height, strain = strain, behavior = "LANDING")

behavior_density <- rbind(flight_densitities, landing_densities)

# SIM only applies for flight
sim_trajectory_CI$behavior <- "FLIGHT"

# now plot
ggplot() +
  geom_ribbon(data = sim_trajectory_CI,
              aes(x = height, ymin = lower, ymax = upper),
              fill = "grey70", alpha = 0.4) +
  geom_line(data = sim_trajectory_CI,
            aes(x = height, y = median, color = "SIM"), linewidth = 3) +
  geom_density(data = behavior_density,
               aes(x = height, color = strain), linewidth = 3) +
  geom_vline(xintercept = 100, color = "red", linewidth = 3) +
  annotate("text", x = 100, y = 0.014, label = "Hut Midpoint", angle = -360, vjust = -1, color = "red") +
  scale_color_manual(values = c("LAB" = "skyblue", "WT" = "darkorange", "SIM" = "black"),
                     breaks = c("LAB", "WT", "SIM"), name = NULL) +
  scale_y_continuous(limits = c(0, 0.02), breaks = c(0.00, 0.01, 0.02)) +
  coord_flip(xlim = c(0, 204.5992)) +
  facet_wrap(~behavior) +
  labs(x = "Height (cm)", y = "Density") +
  theme_bw(base_size = 30) +
  theme(legend.position = "right", legend.direction = "vertical",
        panel.spacing = unit(2, "lines")) 
  
  


## Figure 1 Final Assembly and Export
library(pdftools)
library(magick)

# read PPT-exported PDF
fig1_img <- image_read_pdf("Revised_Figures/Figure1_REVISED.pdf", density = 300)

# write out as .tiff
image_write(fig1_img, "Figure1_REVISED.tiff", format = "tiff", density = 300, compression = "lzw")



#####################################################
#####################################################

## GLMM used to investigate flight differences between LAB and WT

# prep data to be fitted (W treatment subset)
w_tracks <- subset(tracks, color == "W")
w_tracks$strain <- as.factor(w_tracks$strain)
w_tracks$unique.id <- as.factor(w_tracks$unique.id)

w_tracks$strain <- relevel(w_tracks$strain, 
                           ref = "LAB")
w_tracks$sex <- relevel(w_tracks$sex, 
                        ref = "Unfed")


# have to normalize the data
w_tracks <- w_tracks %>% 
  mutate(y_cm = y_cm / 204.5992)

# make y_cm 0s -> 7.239476e-06
w_tracks <- w_tracks %>% 
  mutate(y_cm = if_else(y_cm == 0, 4.740977e-06, y_cm))

# fit model
glmm.intrinsic.strain <- glmmTMB(
  y_cm ~ strain + sex + strain:sex +
    (1|unique.id), 
  data = w_tracks, 
  family = beta_family(link = "logit"), 
  na.action = "na.fail")

### Dredge 
flight.strain.dredge <- dredge(glmm.intrinsic.strain)

# put dredge results in table
flight.strain.table <- as.data.frame(flight.strain.dredge)

# extract additive model from dredge table and get result output
flight.add.model <- get.models(flight.strain.dredge, 2)[[1]]
summary(flight.add.model)

# Get coefficients and 95% CIs
confint(flight.add.model)

# track mean heights (UF, W)
track_uf_w_means <- w_tracks %>%
  group_by(strain) %>%
  summarize(mean_height = mean(y_cm * 204.5992, na.rm = TRUE))


## GLMM used to investigate landing differences between LAB and WT

# subset sticky dataframe to W only
sticky_w <- subset(sticky, color == "W")

sticky_w <- droplevels(sticky_w)

# re-level variables
sticky$strain <- relevel(sticky$strain, 
                         ref = "LAB")

sticky$sex <- relevel(sticky$sex, 
                      ref = "Unfed")


# let's fit the model
landing.glmm.intrinsic <- glmmTMB(
  resting.height ~ strain + sex + strain:sex +
    (1|replicate.id), 
  data = sticky, 
  family = beta_family(link = "logit"), 
  na.action = "na.fail")

# Dredge 
landing.glmm.dredge <- dredge(landing.glmm.intrinsic)

# put dredge results in table
landing.glmm.table <- as.data.frame(landing.glmm.dredge)

# extract additive model from dredge table and get result output
landing.add.model <- get.models(landing.glmm.dredge, 2)[[1]]
summary(landing.add.model)

confint(landing.add.model)

# landing mean heights (UF, W)
landing_uf_w_means <- sticky_w %>%
  group_by(strain) %>%
  summarize(mean_height = mean(resting.height * 192, na.rm = TRUE))



#####################################################
#####################################################

## Track GAM - Tensor (Temp and RH)
# subset track data according to strain and color (WT, W, Males and UFs)
wt.w.tracks <- subset(tracks, strain == "WT" & color == "W" & sex %in% c("Male", "Unfed")) %>% 
  droplevels()

# downsizing coordinates to keep every 10th for each track id
wt.w.tracks.downsize <- wt.w.tracks %>%
  group_by(id) %>%
  filter(row_number() %% 10 == 1) %>%  # Keep rows 1, 11, 21, 31, etc.
  ungroup()

# normalizing y_cm between 0 and 1 to fit beta 
wt.w.tracks.downsize <- wt.w.tracks.downsize %>% 
  mutate(y_cm = y_cm / 204.5992)

# GAM w/ Tensor (Temp & RH) and 3D visualization with improved aesthetics
flight.gam.wt.w <- gam(y_cm ~ te(temp, rh, k = 3) + sex + 
                         s(unique.id, bs = "re"), 
                         data = wt.w.tracks.downsize, 
                         family = betar(link = "logit"))

summary(flight.gam.wt.w)

# 3D raster surface
vis.gam(flight.gam.wt.w, theta = 215, phi = 30,
        view = c("temp", "rh"),
        too.far = 0.1,
        ticktype = "detailed",
        type = "response",
        plot.type = "persp",
        color = "heat",  # or try "heat", "terrain", "cm"
        xlab = "Temperature (°C)",
        ylab = "Relative Humidity (%)",
        zlab = "Flying Height (Normalized 0,1)",
        zlim = c(0,1),
        main = "Mosquito Flying Height Response Surface\n(Temperature × Humidity Interaction)",
        border = NA,        # Remove grid lines for cleaner look
        shade = 0.75,       # Add shading for depth perception
        expand = 0.6,       # Adjust z-axis scaling
        cex.lab = 2.5,      # Larger axis labels
        cex.main = 2,     # Main title size
        cex.axis = 2,     # Axis text size
        ltheta = 45,        # Light source angle
        lphi = 45,          # Light source elevation
        n.grid = 50)       # Higher resolution (default is usually 40)

# add legend
pred_range <- range(predict(flight.gam.wt.w, type = "response"))
legend("bottomright", inset = c(-0.6, -0.01), 
       legend = c("High", "Medium", "Low"),
       fill = rev(heat.colors(3)),
       title = "Flying Height",
       bty = "n", xpd = TRUE,
       cex = 2,
       title.cex = 2)


## Track GAM - RH smooth
model.wt.w.tracks.rh <- gam(y_cm ~ s(rh, k = 3) + sex +
                              s(unique.id, bs = "re"), 
                            data = wt.w.tracks.downsize, 
                            family = betar(link = "logit"))


summary(model.wt.w.tracks.rh)




## Track GAM - Temp smooth
model.wt.w.tracks.temp <- bam(y_cm ~ s(temp, k = 3) + sex +
                                s(unique.id, bs = "re"), 
                              data = wt.w.tracks.downsize, 
                              family = betar(link = "logit"))

summary(model.wt.w.tracks.temp)




## Track GAM - VPD smooth
# estimate VPD
# Calculate saturated vapor pressure (es) using Temp and the Tetens equation
wt.w.tracks.downsize <- wt.w.tracks.downsize %>% 
  mutate(es = 0.6108 * exp(17.27 * temp / (temp + 237.3)))

# Calculate actual vapor pressure (ea) using RH
wt.w.tracks.downsize <- wt.w.tracks.downsize %>% 
  mutate(ea = es * (rh / 100))

# Calculate VPD by doing (es - ea)
wt.w.tracks.downsize <- wt.w.tracks.downsize %>% 
  mutate(vpd = es - ea)

# Fit model
model.wt.w.tracks.vpd <- bam(y_cm ~ s(vpd, k = 3) + sex +
                               s(unique.id, bs = "re"), 
                             data = wt.w.tracks.downsize, 
                             family = betar(link = "logit"))


summary(model.wt.w.tracks.vpd)

# Extract and compare AIC across the four GAM structures
wt.w.track.gam.aic <- AIC(flight.gam.wt.w, model.wt.w.tracks.rh,
                          model.wt.w.tracks.temp, model.wt.w.tracks.vpd)


print(wt.w.track.gam.aic)



## With the flight GAMs complete, now do the same for landing GAMs
# Subset sticky data to just WT and W color (Males and UFs only)
wt.w.sticky <- subset(sticky, strain == "WT" & color == "W" & sex %in% c("Male", "Unfed"))


## Landing GAM - Tensor (mean_temp and mean_rh)
model.wt.w.sticky <- gam(resting.height ~ te(mean_temp, mean_rh, k = 3) + sex, 
                         data = wt.w.sticky, 
                         family = betar(link = "logit"))

summary(model.wt.w.sticky)

vis.gam(model.wt.w.sticky, theta = 215, phi = 30,
        view = c("mean_temp", "mean_rh"),
        too.far = 0.25,
        ticktype = "detailed",
        type = "response",
        plot.type = "persp",
        color = "heat",  # or try "heat", "terrain", "cm"
        contour.col = "white",
        xlab = "Temperature (°C)",
        ylab = "Relative Humidity (%)",
        zlab = "Landing Height (Normalized 0,1)",
        zlim = c(0,1),
        main = "Mosquito Landing Height Response Surface\n(Temperature × Humidity Interaction)",
        border = NA,        # Remove grid lines for cleaner look
        shade = 0.75,       # Add shading for depth perception
        expand = 0.6,       # Adjust z-axis scaling
        cex.lab = 2.5,      # Larger axis labels
        cex.main = 2,     # Main title size
        cex.axis = 2,     # Axis text size
        ltheta = 45,        # Light source angle
        lphi = 45,          # Light source elevation
        n.grid = 50)       # Higher resolution (default is usually 40)

pred_range <- range(predict(model.wt.w.sticky, type = "response"))
legend("bottomright", inset = c(-0.6, -0.01),
       legend = c("High", "Medium", "Low"),
       fill = rev(heat.colors(3)[c(1,2,3)]),
       title = "Landing Height",
       bty = "n", xpd = TRUE,
       cex = 2,
       title.cex = 2)



## Landing GAM - RH smooth only
model.wt.w.sticky.rh <- gam(resting.height ~ s(mean_rh, k = 3) + sex, 
                            data = wt.w.sticky, 
                            family = betar(link = "logit"))


summary(model.wt.w.sticky.rh)


## Landing GAM - Temp smooth only
model.wt.w.sticky.temp <- gam(resting.height ~ s(mean_temp, k = 3) + sex, 
                              data = wt.w.sticky, 
                              family = betar(link = "logit"))

summary(model.wt.w.sticky.temp)



## Landing GAM - VPD smooth only
# Estimating VPD
# Calculate saturated vapor pressure (es) using Temp and the Tetens equation
wt.w.sticky <- wt.w.sticky %>% 
  mutate(es = 0.6108 * exp(17.27 * mean_temp / (mean_temp + 237.3)))

# Calculate actual vapor pressure (ea) using RH
wt.w.sticky <- wt.w.sticky %>%
  mutate(ea = es * (mean_rh / 100))

# Calculate VPD by doing (es - ea)
wt.w.sticky <- wt.w.sticky %>%
  mutate(vpd = es - ea)


## now fit the GAM
model.wt.w.sticky.vpd <- gam(resting.height ~ s(vpd, k = 3) + sex, 
                             data = wt.w.sticky, 
                             family = betar(link = "logit"))


summary(model.wt.w.sticky.vpd)

## now get and compare AIC across the four model structures
wt.w.sticky.gam.aic <- AIC(model.wt.w.sticky, model.wt.w.sticky.rh,
                           model.wt.w.sticky.temp, model.wt.w.sticky.vpd)


print(wt.w.sticky.gam.aic)



### That is all for the GAMs; next is Firth GLM to assess relationship
### between VPD and the probability mosquitoes will land in the upper half of hut (vertically)

# un-normalize resting (landing) height
sticky.vpd <- sticky %>% 
  mutate(resting.height = resting.height * 192)


# create the binary response for landing (i.e. above vs below 1m)
sticky.vpd <- sticky.vpd %>% 
  mutate(upperhalf = ifelse(resting.height > 96, 1, 0))

# compute VPD for entire dataset
sticky.vpd <- sticky.vpd %>% 
  mutate(es = 0.6108 * exp(17.27 * mean_temp / (mean_temp + 237.3)))

sticky.vpd <- sticky.vpd %>%
  mutate(ea = es * (mean_rh / 100))

sticky.vpd <- sticky.vpd %>%
  mutate(vpd = es - ea)


# subset to WT, M and UF, B and W (uniform sticky traps)
sticky.vpd.wt.b.w <- subset(sticky.vpd, strain == "WT" & sex %in% c("Male", "Unfed") & color %in% c("B", "W"))

sticky.vpd.wt.b.w <- droplevels(sticky.vpd.wt.b.w)

# check for separation for the subset
with(sticky.vpd.wt.b.w, tapply(upperhalf, cut(vpd, 8), mean))

with(sticky.vpd.wt.b.w, table(cut(vpd, 8), upperhalf))


# Given separation issues across VPD bins (quasi-complete separation),
## A Firth GLM would be more appropriate 
library(brglm2)

# set reference levels for model fitting
sticky.vpd.wt.b.w$color <- relevel(factor(sticky.vpd.wt.b.w$color), ref = "W")
sticky.vpd.wt.b.w$sex <- relevel(factor(sticky.vpd.wt.b.w$sex), ref = "Unfed")

landing.firth.ab.1 <- glm(upperhalf ~ vpd + color + sex,
                          family = binomial,
                          method = "brglmFit",
                          type = "AS_mean",
                          data = sticky.vpd.wt.b.w)

summary(landing.firth.ab.1)

# Prediction grid for plotting (W color only)
pred_grid2 <- expand.grid(
  vpd = seq(min(sticky.vpd.wt.b.w$vpd), max(sticky.vpd.wt.b.w$vpd), length = 100),
  color = factor("W", levels = levels(sticky.vpd.wt.b.w$color)),
  sex = factor("Unfed", levels = levels(sticky.vpd.wt.b.w$sex)))

# Predictions on link scale with standard errors
pred_link <- predict(landing.firth.ab.1, newdata = pred_grid2, type = "link", se.fit = TRUE)

crit <- 1.96
pred_grid2 <- pred_grid2 %>%
  mutate(
    fit_link   = pred_link$fit,
    se_link    = pred_link$se.fit,
    lower_link = fit_link - crit * se_link,
    upper_link = fit_link + crit * se_link,
    fit_prob   = plogis(fit_link),
    lower_prob = plogis(lower_link),
    upper_prob = plogis(upper_link))

# Plot - W color only
ggplot(pred_grid2, aes(x = vpd, y = fit_prob)) +
  geom_line(linewidth = 1.5, color = "black") +
  geom_ribbon(aes(ymin = lower_prob, ymax = upper_prob), alpha = 0.2, color = NA, fill = "grey50") +
  scale_x_continuous(breaks = seq(0, ceiling(max(pred_grid2$vpd) * 2) / 2, by = 0.5)) +
  ylim(0,1) +
  labs( x = "VPD (kPa)",
        y = "Prob. of Landing in Upper Half",
        title = "White Traps (W)") +
  theme_bw(base_size = 45) +
  theme(plot.title = element_text(hjust = 0.5))


# Plot - B color only
# Prediction grid for plotting (B color only)
pred_grid3 <- expand.grid(
  vpd = seq(min(sticky.vpd.wt.b.w$vpd), max(sticky.vpd.wt.b.w$vpd), length = 100),
  color = factor("B", levels = levels(sticky.vpd.wt.b.w$color)),
  sex = factor("Unfed", levels = levels(sticky.vpd.wt.b.w$sex)))

# Predictions on link scale with standard errors
pred_link3 <- predict(landing.firth.ab.1, newdata = pred_grid3, type = "link", se.fit = TRUE)

crit <- 1.96
pred_grid3 <- pred_grid3 %>%
  mutate(
    fit_link   = pred_link3$fit,
    se_link    = pred_link3$se.fit,
    lower_link = fit_link - crit * se_link,
    upper_link = fit_link + crit * se_link,
    fit_prob   = plogis(fit_link),
    lower_prob = plogis(lower_link),
    upper_prob = plogis(upper_link))

# Plot - W color only
ggplot(pred_grid3, aes(x = vpd, y = fit_prob)) +
  geom_line(linewidth = 1.5, color = "black") +
  geom_ribbon(aes(ymin = lower_prob, ymax = upper_prob), alpha = 0.2, color = NA, fill = "grey50") +
  scale_x_continuous(breaks = seq(0, ceiling(max(pred_grid2$vpd) * 2) / 2, by = 0.5)) +
  ylim(0,1) +
  labs( x = "VPD (kPa)",
        y = "Prob. of Landing in Upper Half",
        title = "Black Traps (B)") +
  theme_bw(base_size = 45) +
  theme(plot.title = element_text(hjust = 0.5))


## Additional VPD summaries
# estimating VPD for all HOBO logger measurements
# Calculate saturated vapor pressure (es) using Temp and the Tetens equation
hobo <- hobo %>% 
  mutate(es = 0.6108 * exp(17.27 * Temp / (Temp + 237.3)))

# Calculate actual vapor pressure (ea) using RH
hobo <- hobo %>% 
  mutate(ea = es * (RH / 100))

# Calculate VPD by doing (es - ea)
hobo <- hobo %>% 
  mutate(vpd = es - ea)

# VPD and time of day collinearity
kruskal.test(vpd ~ Replicate, data = hobo)

# VPD by TOD - Direction and Magnitude
hobo %>% 
  group_by(HOBO, Replicate) %>%
  summarize(mean_vpd = mean(vpd),
            sd_vpd = sd(vpd),
            median_vpd = median(vpd),
            n = n(),
            .groups = "drop")



## Figure 2 Final Assembly and Export
# read PPT-exported PDF
fig2_img <- image_read_pdf("Revised_Figures/Figure2_REVISED.pdf", density = 300)

# write out as .tiff
image_write(fig2_img, "Figure2_REVISED.tiff", format = "tiff", density = 300, compression = "lzw")




#####################################################
#####################################################


## Color Attraction Strongly Modulates Low Height Preference

## Produce histograms for each color (for tracks and sticky)
## then extract frequencies of heights

## Tracks
# Subset each color by strain
# B
B_tracks_lab <- subset(tracks, strain == "LAB" & color == "B") %>%
  droplevels()
B_tracks_wt <- subset(tracks, strain == "WT" & color == "B") %>%
  droplevels()

# BW
BW_tracks_lab <- subset(tracks, strain == "LAB" & color == "BW") %>%
  droplevels()
BW_tracks_wt <- subset(tracks, strain == "WT" & color == "BW") %>%
  droplevels()

# WB
WB_tracks_lab <- subset(tracks, strain == "LAB" & color == "WB") %>%
  droplevels()
WB_tracks_wt <- subset(tracks, strain == "WT" & color == "WB") %>%
  droplevels()

# W
W_tracks_lab <- subset(tracks, strain == "LAB" & color == "W") %>%
  droplevels()
W_tracks_wt <- subset(tracks, strain == "WT" & color == "W") %>%
  droplevels()



# B histograms
# lab
B_lab_track_hist <- ggplot(B_tracks_lab, aes(x = y_cm, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the lab computed data
B_lab_track_freq <- ggplot_build(B_lab_track_hist)$data[[1]]

B_lab_track_freq <- B_lab_track_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "LAB",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "B") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)


# wt
B_wt_track_hist <- ggplot(B_tracks_wt, aes(x = y_cm, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the wt computed data
B_wt_track_freq <- ggplot_build(B_wt_track_hist)$data[[1]]

B_wt_track_freq <- B_wt_track_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "WT",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "B") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)



# BW histograms
# lab
BW_lab_track_hist <- ggplot(BW_tracks_lab, aes(x = y_cm, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the lab computed data
BW_lab_track_freq <- ggplot_build(BW_lab_track_hist)$data[[1]]

BW_lab_track_freq <- BW_lab_track_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "LAB",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "BW") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)


# wt
BW_wt_track_hist <- ggplot(BW_tracks_wt, aes(x = y_cm, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the wt computed data
BW_wt_track_freq <- ggplot_build(BW_wt_track_hist)$data[[1]]

BW_wt_track_freq <- BW_wt_track_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "WT",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "BW") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)



# WB histograms
# lab
WB_lab_track_hist <- ggplot(WB_tracks_lab, aes(x = y_cm, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the lab computed data
WB_lab_track_freq <- ggplot_build(WB_lab_track_hist)$data[[1]]

WB_lab_track_freq <- WB_lab_track_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "LAB",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "WB") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)


# wt
WB_wt_track_hist <- ggplot(WB_tracks_wt, aes(x = y_cm, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the wt computed data
WB_wt_track_freq <- ggplot_build(WB_wt_track_hist)$data[[1]]

WB_wt_track_freq <- WB_wt_track_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "WT",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "WB") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)




# W histograms
# lab
W_lab_track_hist <- ggplot(W_tracks_lab, aes(x = y_cm, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the lab computed data
W_lab_track_freq <- ggplot_build(W_lab_track_hist)$data[[1]]

W_lab_track_freq <- W_lab_track_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "LAB",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "W") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)


# wt
W_wt_track_hist <- ggplot(W_tracks_wt, aes(x = y_cm, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the wt computed data
W_wt_track_freq <- ggplot_build(W_wt_track_hist)$data[[1]]

W_wt_track_freq <- W_wt_track_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "WT",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "W") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)

# add BEHAVIOR column to all track freq dfs
# lab
B_lab_track_freq <- B_lab_track_freq %>% 
  mutate(behavior = "FLIGHT")

BW_lab_track_freq <- BW_lab_track_freq %>% 
  mutate(behavior = "FLIGHT")

WB_lab_track_freq <- WB_lab_track_freq %>% 
  mutate(behavior = "FLIGHT")

W_lab_track_freq <- W_lab_track_freq %>% 
  mutate(behavior = "FLIGHT")

# wt
B_wt_track_freq <- B_wt_track_freq %>% 
  mutate(behavior = "FLIGHT")

BW_wt_track_freq <- BW_wt_track_freq %>% 
  mutate(behavior = "FLIGHT")

WB_wt_track_freq <- WB_wt_track_freq %>% 
  mutate(behavior = "FLIGHT")

W_wt_track_freq <- W_wt_track_freq %>% 
  mutate(behavior = "FLIGHT")

##### Write out all!

write_csv(B_lab_track_freq,"Lab_Black_track_height_freq.csv") 
write_csv(BW_lab_track_freq,"Lab_BlackWhite_track_height_freq.csv") 
write_csv(WB_lab_track_freq,"Lab_WhiteBlack_track_height_freq.csv") 
write_csv(W_lab_track_freq,"Lab_White_track_height_freq.csv") 


write_csv(B_wt_track_freq,"WT_Black_track_height_freq.csv") 
write_csv(BW_wt_track_freq,"WT_BlackWhite_track_height_freq.csv") 
write_csv(WB_wt_track_freq,"WT_WhiteBlack_track_height_freq.csv") 
write_csv(W_wt_track_freq,"WT_White_track_height_freq.csv") 


## Repeat above for landing (sticky trap) data
# B
B_sticky_lab <- subset(sticky, strain == "LAB" & color == "B") %>%
  mutate(resting.height = resting.height * 192) %>%
  droplevels()
  
B_sticky_wt <- subset(sticky, strain == "WT" & color == "B") %>%
  mutate(resting.height = resting.height * 192) %>%
  droplevels()

# BW
BW_sticky_lab <- subset(sticky, strain == "LAB" & color == "BW") %>%
  mutate(resting.height = resting.height * 192) %>%
  droplevels()
  
BW_sticky_wt <- subset(sticky, strain == "WT" & color == "BW") %>% 
  mutate(resting.height = resting.height * 192) %>%
  droplevels()

# WB
WB_sticky_lab <- subset(sticky, strain == "LAB" & color == "WB") %>%
  mutate(resting.height = resting.height * 192) %>%
  droplevels()

WB_sticky_wt <- subset(sticky, strain == "WT" & color == "WB") %>%
  mutate(resting.height = resting.height * 192) %>%
  droplevels()

# W 
W_sticky_lab <- subset(sticky, strain == "LAB" & color == "W") %>%
  mutate(resting.height = resting.height * 192) %>%
  droplevels()

W_sticky_wt <- subset(sticky, strain == "WT" & color == "W") %>%
  mutate(resting.height = resting.height * 192) %>%
  droplevels()




# B histograms
# lab
B_lab_sticky_hist <- ggplot(B_sticky_lab, aes(x = resting.height, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the lab computed data
B_lab_sticky_freq <- ggplot_build(B_lab_sticky_hist)$data[[1]]

B_lab_sticky_freq <- B_lab_sticky_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "LAB",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "B") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)


# wt
B_wt_sticky_hist <- ggplot(B_sticky_wt, aes(x = resting.height, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the wt computed data
B_wt_sticky_freq <- ggplot_build(B_wt_sticky_hist)$data[[1]]

B_wt_sticky_freq <- B_wt_sticky_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "WT",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "B") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)



# BW histograms
# lab
BW_lab_sticky_hist <- ggplot(BW_sticky_lab, aes(x = resting.height, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the lab computed data
BW_lab_sticky_freq <- ggplot_build(BW_lab_sticky_hist)$data[[1]]

BW_lab_sticky_freq <- BW_lab_sticky_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "LAB",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "BW") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)


# wt
BW_wt_sticky_hist <- ggplot(BW_sticky_wt, aes(x = resting.height, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the wt computed data
BW_wt_sticky_freq <- ggplot_build(BW_wt_sticky_hist)$data[[1]]

BW_wt_sticky_freq <- BW_wt_sticky_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "WT",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "BW") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)



# WB histograms
# lab
WB_lab_sticky_hist <- ggplot(WB_sticky_lab, aes(x = resting.height, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the lab computed data
WB_lab_sticky_freq <- ggplot_build(WB_lab_sticky_hist)$data[[1]]

WB_lab_sticky_freq <- WB_lab_sticky_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "LAB",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "WB") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)


# wt
WB_wt_sticky_hist <- ggplot(WB_sticky_wt, aes(x = resting.height, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the wt computed data
WB_wt_sticky_freq <- ggplot_build(WB_wt_sticky_hist)$data[[1]]

WB_wt_sticky_freq <- WB_wt_sticky_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "WT",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "WB") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)




# W histograms
# lab
W_lab_sticky_hist <- ggplot(W_sticky_lab, aes(x = resting.height, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the lab computed data
W_lab_sticky_freq <- ggplot_build(W_lab_sticky_hist)$data[[1]]

W_lab_sticky_freq <- W_lab_sticky_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "LAB",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "W") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)


# wt
W_wt_sticky_hist <- ggplot(W_sticky_wt, aes(x = resting.height, fill = strain)) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), 
                 position = "identity", alpha = 0.7, bins = 20)

# Extract the wt computed data
W_wt_sticky_freq <- ggplot_build(W_wt_sticky_hist)$data[[1]]

W_wt_sticky_freq <- W_wt_sticky_freq %>%
  select(x, y, fill, xmin, xmax, count) %>%
  rename(strain = fill, 
         bin_center = x, 
         bin_start = xmin,
         bin_end = xmax, 
         relative_freq = y) %>% 
  mutate(strain = case_when(
    strain == "#F8766D" ~ "WT",
    TRUE ~ strain  # keeps any other colors unchanged
  ),
  color = "W") %>% 
  select(strain, color, bin_start, bin_center, bin_end, count, relative_freq)

# create a new column for BEHAVIOR in all landing dfs
# lab
B_lab_sticky_freq <- B_lab_sticky_freq %>% 
  mutate(behavior = "LANDING")

BW_lab_sticky_freq <- BW_lab_sticky_freq %>% 
  mutate(behavior = "LANDING")

WB_lab_sticky_freq <- WB_lab_sticky_freq %>% 
  mutate(behavior = "LANDING")

W_lab_sticky_freq <- W_lab_sticky_freq %>% 
  mutate(behavior = "LANDING")

# wt
B_wt_sticky_freq <- B_wt_sticky_freq %>% 
  mutate(behavior = "LANDING")

BW_wt_sticky_freq <- BW_wt_sticky_freq %>% 
  mutate(behavior = "LANDING")

WB_wt_sticky_freq <- WB_wt_sticky_freq %>% 
  mutate(behavior = "LANDING")

W_wt_sticky_freq <- W_wt_sticky_freq %>% 
  mutate(behavior = "LANDING")

# Write out all

write_csv(B_lab_sticky_freq,"Lab_Black_sticky_height_freq.csv") 
write_csv(BW_lab_sticky_freq,"Lab_BlackWhite_sticky_height_freq.csv") 
write_csv(WB_lab_sticky_freq,"Lab_WhiteBlack_sticky_height_freq.csv") 
write_csv(W_lab_sticky_freq,"Lab_White_sticky_height_freq.csv") 


write_csv(B_wt_sticky_freq,"WT_Black_sticky_height_freq.csv") 
write_csv(BW_wt_sticky_freq,"WT_BlackWhite_sticky_height_freq.csv") 
write_csv(WB_wt_sticky_freq,"WT_WhiteBlack_sticky_height_freq.csv") 
write_csv(W_wt_sticky_freq,"WT_White_sticky_height_freq.csv") 



## With height frequencies calculates, bind all dataframes
# to produce LOESS plots for Figure 3

# rbind all track and sticky freq dfs
df <- rbind(B_lab_track_freq, BW_lab_track_freq, WB_lab_track_freq, W_lab_track_freq,
            B_lab_sticky_freq, BW_lab_sticky_freq, WB_lab_sticky_freq, W_lab_sticky_freq,
            B_wt_track_freq, BW_wt_track_freq, WB_wt_track_freq, W_wt_track_freq,
            B_wt_sticky_freq, BW_wt_sticky_freq, WB_wt_sticky_freq, W_wt_sticky_freq)

# to stratify by strain and behavior, mutate new column
df <- df %>% 
  mutate(color = factor(color,
                        levels = c("W", "B", "WB", "BW"),
                        labels = c("White", "Black", "White:Black", "Black:White")))


# now plot using LOESS smoothing
ggplot(df, aes(x = bin_center, y = relative_freq, color = behavior, fill = behavior)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.3, alpha = 0.2, linewidth = 1.5) +
  facet_grid(strain ~ color) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "grey40", linewidth = 1) +
  scale_x_continuous(name = "Height (cm - bin center)") +
  scale_y_continuous(name = "Relative Frequency") +
  scale_color_manual(values = c("FLIGHT" = "#01665E","LANDING" = "#8E0152")) +
  scale_fill_manual(values = c("FLIGHT" = "#01665E", "LANDING" = "#8E0152")) +
  coord_flip(ylim = c(0, 0.3)) +
  theme_bw(base_size = 45) +
  theme(strip.background = element_rect(fill = "grey90", color = NA),
        panel.grid.major = element_line(color = alpha("grey80", 0.25)),  # 50% transparency
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legent.text = element_text(size = 30),
        axis.title = element_text(size = 45),
        axis.text = element_text(size = 30))



## Color GLMMs for both flight and landing 

# subset data (initially WT and Unfed)
wt_uf_tracks <- subset(tracks, strain == "WT" & sex == "Unfed") %>%
  droplevels()

# have to normalize the data
wt_uf_tracks <- wt_uf_tracks %>% 
  mutate(y_cm = y_cm / 204.5992)

# make y_cm 0s -> 8.088986e-06
wt_uf_tracks <- wt_uf_tracks %>% 
  mutate(y_cm = if_else(y_cm == 0, 8.088986e-06, y_cm))

# downsizing coordinates to keep every 10th for each track id
wt_uf_tracks.downsize <- wt_uf_tracks %>%
  group_by(id) %>%
  filter(row_number() %% 10 == 1) %>%  # Keep rows 1, 11, 21, 31, etc.
  ungroup()

# set W as reference color level
wt_uf_tracks.downsize$color <- relevel(wt_uf_tracks.downsize$color, 
                         ref = "W")


# fit model
glmm.extrinsic.wt <- glmmTMB(
  y_cm ~ color + (1 | unique.id), 
  data = wt_uf_tracks.downsize, 
  family = beta_family(link = "logit"), 
  na.action = "na.fail")


# Get model summary, coefficients and 95% CIs
summary(glmm.extrinsic.wt)
confint(glmm.extrinsic.wt)


## fitting this same GLMM but for landing
# subset data to UT WT
wt_uf_sticky <- subset(sticky, strain == "WT" & sex == "Unfed") %>%
  droplevels()

# set W as reference color level
wt_uf_sticky$color <- relevel(wt_uf_sticky$color, 
                              ref = "W")


# fit model
color.landing.glmm <- glmmTMB(
  resting.height ~ color + (1 | replicate.id), 
  data = wt_uf_sticky, 
  family = beta_family(link = "logit"), 
  na.action = "na.fail")

# Get model summary, coefficients and 95% CIs
summary(color.landing.glmm)
confint(color.landing.glmm)


#####################################################
#####################################################


### Considering Sex and Physiological Status Alongside Color Modulation

# Fit a new GLMM combining strain and sex/status with color as fixed effects

# Flight first
# downsize tracks data
tracks.slim <- tracks %>%
  group_by(id) %>%
  filter(row_number() %% 10 == 1) %>%  # Keep rows 1, 11, 21, 31, etc.
  ungroup()

# normalize y_cm to be between 0 and 1
tracks.slim <- tracks.slim %>% 
  mutate(y_cm = y_cm / 204.5992)

# make y_cm 0s -> 2.785935e-06
tracks.slim <- tracks.slim %>% 
  mutate(y_cm = if_else(y_cm == 0, 2.785935e-06, y_cm))

# set reference levels (W for color, LAB for strain, and Unfed for sex)
tracks.slim$color <- relevel(tracks.slim$color, 
                             ref = "W")

tracks.slim$strain <- relevel(tracks.slim$strain, 
                              ref = "LAB")

tracks.slim$sex <- relevel(tracks.slim$sex, 
                           ref = "Unfed")


## fitting a "global" GLMM
color.combined.glmm <- glmmTMB(
  y_cm ~ strain + sex + strain:sex + color + color:strain + color:sex +
    (1 | unique.id), 
  data = tracks.slim, 
  family = beta_family(link = "logit"), 
  na.action = "na.fail")

## dredge the model
color.combined.dredge <- dredge(color.combined.glmm)

# put dredge results in table
color.combined.table <- as.data.frame(color.combined.dredge)

# extract top model from dredge table and get result output
flight.combined.model <- get.models(color.combined.dredge, 6)[[1]]
summary(flight.combined.model)
confint(flight.combined.model)




### Now do the same for landing data
# set appropriate reference levels for color (W), strain (LAB), sex (UF)
sticky$color <- relevel(sticky$color,
                        ref = "W")

sticky$strain <- relevel(sticky$strain,
                         ref = "LAB")

sticky$sex <- relevel(sticky$sex,
                      ref = "Unfed")


## fitting a "global" GLMM
landing.combined.glmm <- glmmTMB(
  resting.height ~ strain + sex + strain:sex + color + strain:color + sex:color +
    (1 | replicate.id), 
  data = sticky, 
  family = beta_family(link = "logit"), 
  na.action = "na.fail")

## dredge the model
landing.combined.dredge <- dredge(landing.combined.glmm)

# put dredge results in table
landing.combined.table <- as.data.frame(landing.combined.dredge)

# extract top model from dredge table and get result output
landing.combined.model <- get.models(landing.combined.dredge, 1)[[1]]
summary(landing.combined.model)
confint(landing.combined.model)


#####################################################
#####################################################

### Progeny from Reciprocal Crosses Resemble Landing Height Preferences of WT Parents

## Comparisons for Cross A (LAB M x WT F): 
## 1) progeny (UF & M), 2) corresponding WT parent, & 3) corresponding LAB parent

## Comparisons for Cross B (LAB F x WT M: 
## 1) progeny (UF & M), 2) corresponding WT parent, & 3) corresponding LAB parent


## First Cross A (Lab Males ; WT Females)

# subset parent landing data (LAB M x WT F)
lab_males <- subset(sticky, strain == "LAB" & sex == "Male" & color == "B") %>%
  droplevels()

wt_females <- subset(sticky, strain == "WT" & sex == "Unfed" & color == "B") %>% 
  droplevels()

# keeping appropriate columns
lab_males <- lab_males[, c("strain", "color", "sex", "replicate", "resting.height")]
wt_females <- wt_females[, c("strain", "color", "sex", "replicate", "resting.height")]

# un-normalize landing height
lab_males <- lab_males %>% 
  mutate(resting.height = resting.height * 192)

wt_females <- wt_females %>% 
  mutate(resting.height = resting.height * 192)

# progeny
# read in data
crosses <- read_csv("crosses.2024.resting.updated.4.2025.csv")

# update column names
colnames(crosses) <- c("date", "time", "strain", "color", "sex", "replicate", "trap", "resting.height")

# convert landing height from m to cm
crosses <- crosses %>% 
  mutate(resting.height = resting.height * 100)

# set factors
crosses$strain <- as.factor(crosses$strain)
crosses$color <- as.factor(crosses$color)
crosses$sex <- as.factor(crosses$sex)
crosses$replicate <- as.factor(crosses$replicate)

# subset data for Cross A (LAB M x WT F)
cross_a <- subset(crosses, strain == "A") %>%
  droplevels()

# keeping necessary columns to match parent dfs
cross_a <- cross_a[, c("strain", "color", "sex", "replicate", "resting.height")]

# bind progeny and parent data
cross_a_parent_progeny <- rbind(lab_males, wt_females, cross_a)


# make labels more clear for parents and progeny
cross_a_parent_progeny <-cross_a_parent_progeny %>%
  mutate(strain = case_when(
    strain == "LAB" ~ "LAB_Male",
    strain == "WT" ~ "WT_Female",
    strain == "A" ~ "Progeny" # keeps all other values unchanged
  ))

# setting correct display order
cross_a_parent_progeny$strain <- factor(cross_a_parent_progeny$strain, 
                                        levels = c("LAB_Male", "WT_Female", "Progeny"))


## Fit GLMM 
# new column for replicate id
cross_a_parent_progeny <- cross_a_parent_progeny %>% 
  mutate(replicate.id = paste(strain, color, sex, replicate, sep = "_"))

cross_a_parent_progeny$replicate.id <- as.factor(cross_a_parent_progeny$replicate.id)

# normalize resting height between 0 and 1
cross_a_parent_progeny <- cross_a_parent_progeny %>% 
  mutate(resting.height = resting.height / 192)

# changing 0s to 0.002604167
cross_a_parent_progeny$resting.height[cross_a_parent_progeny$resting.height == 0.000000000] <- 0.002604167

# set cross as the reference
cross_a_parent_progeny$strain <- relevel(cross_a_parent_progeny$strain, 
                                         ref = "Progeny")



# okay now I can fit the model
cross.a.glmm <- glmmTMB(resting.height ~ strain + (1 | replicate.id), 
  data = cross_a_parent_progeny, 
  family = beta_family(link = "logit"), 
  na.action = "na.fail")

summary(cross.a.glmm)


## boxplot with significance markings
ggplot(cross_a_parent_progeny, aes(x = strain, y = resting.height *192, color = strain)) +
  geom_boxplot(size = 5, width = 0.75, fill = NA) +
  scale_color_manual(values = c("LAB_Male" = "skyblue", "WT_Female" = "darkorange", "Progeny" = "darkred")) +
  scale_x_discrete(labels = c("LAB_Male" = "LAB Male",
                              "WT_Female" = "WT Female",
                              "Progeny" = "Progeny")) +
  theme_bw(base_size = 45) +
  labs(title = "Cross 1: Lab Male x WT Female") +
  xlab("") + 
  ylab("Landing Height (cm)") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("Progeny", "LAB_Male")), 
              annotations = "*", textsize = 25, vjust = 0.65, 
              y_position = 210, size = 2.5, color = "red") +
  geom_signif(comparisons = list(c("Progeny", "WT_Female")), 
              annotations = "ns", textsize = 15, vjust = 0,  
              y_position = 195, size = 2.5, color = "black")
  

## Now Cross B (LAB Female x WT Male)

## load in datasets
# parent
lab_females <- subset(sticky, strain == "LAB" & sex == "Unfed" & color == "B") %>%
  droplevels()

wt_males <- subset(sticky, strain == "WT" & sex == "Male" & color == "B") %>%
  droplevels()

# keeping necessary columns
lab_females <- lab_females[, c("strain", "color", "sex", "replicate", "resting.height")]
wt_males <- wt_males[, c("strain", "color", "sex", "replicate", "resting.height")]

# un-normalizing landing height
lab_females <- lab_females %>% 
  mutate(resting.height = resting.height * 192)

wt_males <- wt_males %>% 
  mutate(resting.height = resting.height * 192)


# progeny
cross_b <- subset(crosses, strain == "B") %>%
  droplevels()

# keep same columns as parent dfs
cross_b <- cross_b[, c("strain", "color", "sex", "replicate", "resting.height")]

cross_b_parent_progeny <- rbind(lab_females, wt_males, cross_b)


# make labels more clear for parents and progeny
cross_b_parent_progeny <-cross_b_parent_progeny %>%
  mutate(strain = case_when(
    strain == "LAB" ~ "LAB_Female",
    strain == "WT" ~ "WT_Male",
    strain == "B" ~ "Progeny"))

# setting correct display order
cross_b_parent_progeny$strain <- factor(cross_b_parent_progeny$strain, 
                                        levels = c("LAB_Female", "WT_Male", "Progeny"))


## Fit GLMM 

# new column for replicate id
cross_b_parent_progeny <- cross_b_parent_progeny %>% 
  mutate(replicate.id = paste(strain, color, sex, replicate, sep = "_"))

cross_b_parent_progeny$replicate.id <- as.factor(cross_b_parent_progeny$replicate.id)

# normalize resting height between 0 and 1
cross_b_parent_progeny <- cross_b_parent_progeny %>% 
  mutate(resting.height = resting.height / 192)

# changing 0s to 0.002604167
cross_b_parent_progeny$resting.height[cross_b_parent_progeny$resting.height == 0.000000000] <- 0.002604167

# set cross as the reference
cross_b_parent_progeny$strain <- relevel(cross_b_parent_progeny$strain, 
                                         ref = "Progeny")


# okay now I can fit the model
cross.b.glmm <- glmmTMB(resting.height ~ strain + (1 | replicate.id), 
  data = cross_b_parent_progeny, 
  family = beta_family(link = "logit"), 
  na.action = "na.fail")

summary(cross.b.glmm)

## boxplot with significance markings
ggplot(cross_b_parent_progeny, aes(x = strain, y = resting.height * 192, color = strain)) +
  geom_boxplot(size = 5, width = 0.75, fill = NA) +
  scale_color_manual(values = c("LAB_Female" = "skyblue", "WT_Male" = "darkorange", "Progeny" = "darkred")) +
  scale_x_discrete(labels = c("LAB_Female" = "LAB Female",
                              "WT_Male" = "WT Male",
                              "Progeny" = "Progeny")) +
  theme_bw(base_size = 45) +
  labs(title = "Cross 2: Lab Female x WT Male") +
  xlab("") + 
  ylab("Landing Height (cm)") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("Progeny", "LAB_Female")), 
              annotations = "***", textsize = 25, vjust = 0.65, 
              y_position = 210, size = 2.5, color = "red") +
  geom_signif(comparisons = list(c("Progeny", "WT_Male")), 
              annotations = "ns", textsize = 15, vjust = 0, 
              y_position = 195, size = 2.5, color = "black")
  













