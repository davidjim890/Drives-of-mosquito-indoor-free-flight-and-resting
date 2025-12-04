Drivers of mosquito indoor free-flight and resting

This repository contains all datasets (raw and processed) and the full R script used to analyze and visualize experimental data for the manuscript Drivers of mosquito indoor free-flight and resting. All materials are provided to ensure full transparency and reproducibility.

Files within the repository:
organized_pfmd_track_datasets_raw.zip     # Raw PFMD mosquito trajectory CSVs
sticky_hobo_allcolors_raw.csv             # Raw sticky-trap resting behavior data
Processed_Data.zip                        # All processed datasets (described below)
Jimenez-Vallejo_et_al_2025_mosquito_behavior_script.R   # Full analysis script

###########################################

Raw Data
organized_pfmd_track_datasets_raw.zip 
    Contains all raw mosquito trajectory CSVs from the PFMD tracking system.
    Files are organized by experimental replicate.  
    Represents unaltered data as exported directly from the PFMD portal.

sticky_hobo_allcolors_raw.csv
    Contains all raw sticky-trap resting observations across all replicates.
    Includes mosquito captures by trap color and replicate information.


###########################################

Processed Data (inside “Processed_Data.zip”)

  Simulated Trajectory Data
    Files with “SIM” in the filename represent simulated mosquito trajectories.
    Generated using step-length and turning-angle distributions from “White-only” sticky trap treatments (unbiased behavior).
    Produced using a random-walk model (see Methods for a detailed description of model).


  Processed PFMD Trajectory Data
  "cleaned_track_data_1.21.2025.csv"
          All PFMD trajectory CSVs merged into one file.
          Includes all relevant design factors.
          Datetimes were manually corrected to match true replicate timing.
          
  "track_hobo_merged_8.19.2025.csv"
          PFMD trajectory data merged with HOBO logger environmental data.
          Millisecond-precision timestamps enable accurate matching of temperature and RH.
          
  "track_vpd_data.csv"
          Same as above but includes Vapor Pressure Deficit (VPD) calculated for each trajectory coordinate (see Methods for VPD equation used).


  HOBO Logger Data
  "hobo_all_data.csv"
          Raw HOBO logger temperature and RH measurements for all replicates.

  "hobo_data_clean.csv"
          Cleaned and annotated version of the HOBO dataset.
          Includes design-factor variables for merging with PFMD and sticky-trap data.


  Sticky Trap + HOBO Merged Data
  sticky_hobo_model_ready_8_19_2024.csv
          Sticky-trap data merged with HOBO logger data.
          Each captured mosquito is assigned estimated temperature and RH based on height and time.

  sticky_vpd_data.csv
          Same as above, but with VPD added.


############################################

Analysis Script
Jimenez-Vallejo_et_al_2025_mosquito_behavior_script.R
  This R script contains all analyses used in the manuscript, including:
      Data cleaning and preprocessing
      Random-walk simulations
      Statistical modeling
      Integration of environmental covariates
      Figure generation for manuscript and supplementary materials
      The script includes extensive comments to support reproducibility.


###########################################

License
Creative Commons Attributions 4.0 International



############################################

Citation

If you use these data or code, please cite:

Jiménez-Vallejo et al. (2025). Drivers of mosquito indoor free-flight and resting.
Zenodo. DOI: https://doi.org/10.5281/zenodo.17817358













          
