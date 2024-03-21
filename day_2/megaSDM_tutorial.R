#### MegaSDM Demonstration ####
# This package can efficiently create and project species distribution models
# using the MaxEnt framework and parallel processing. It can find and download
# occurrence data for a list of species on GBIF (Global Biodiversity
# Information Facility), environmentally subsample the occurrences to mitigate
# spatial bias, generate background (pseudo-absence) points, train the model and
# project it to different times (incorporating dispersal rate of each species and
# intermediate range fluctuations), and create species richness maps for each time
# period and taxon.

#Install megaSDM
#install.packages("devtools")
library(devtools)

devtools::install_github("brshipley/megaSDM", build_vignettes = TRUE)
library(megaSDM)
library(terra)

#If you haven't used a java GUI before, ensure that maxent.jar can open:
#Navigate to the file path shown here, and click on the file to open. If
#a dialog box open without any issues, you're good to go!
MaxEntFile <- list.files(system.file("extdata", package = "megaSDM"),
                          pattern = ".jar$",
                          full.names = TRUE, recursive = TRUE)

# Set up Working Directory------------------------------------------------------

WD <- "/Users/tjs/R_packages/R_projects/spatialSampling"
setwd(WD)

## Define Output for Data/Models------------------------------------------------
envoutput <- "TestRun"

## Load and Manage Climate Rasters----------------------------------------------
#megaSDM uses a "training area", where the model is trained on, and a
#"study area", where the model is projected to. These can be different, but
#they don't have to be, so in this example, we are training and projection on
#the same area.

#To train the model, we'll use the current/modern data to pair with the modern
#occurrences.
input_TA <- list.files("ClimateData/Current_0k",
                       pattern = ".grd$",
                       full.names = TRUE)

#We also need to direct R towards the historical/paleoclimate layers:
LGM <- list.files("ClimateData/LGM_21k",
                  pattern = ".grd$",
                  full.names = TRUE)

BA <- list.files("ClimateData/BA_15k",
                 pattern = ".grd$",
                 full.names = TRUE)

EH <- list.files("ClimateData/EH_11k",
                 pattern = ".grd$",
                 full.names = TRUE)

LH <- list.files("ClimateData/LH_4k",
                 pattern = ".grd$",
                 full.names = TRUE)


# Here we define the extent of the training and study regions in
# c(xmin, xmax, ymin, ymax) form and format the environmental layrs to make sure
#they are integrable.

TSEnv <- TrainStudyEnv(input_TA = input_TA,
                       output = envoutput,
                       clipTrain = c(-91.5, -75, 25.5, 36),
                       clipStudy = c(-91.5, -75, 25.5, 36))

climlayers <- list(LH, EH, BA, LGM)

# The "time_periods" argument must contain the current time (the time of the training
# and study rasters) first and then the time periods for the forecast/hindcast.

PredictEnv(studylayers = TSEnv$study,
           futurelayers = climlayers,
           time_periods = c(2000, -4000, -11000, -15000, -21000),
           output = envoutput,
           scenario_name = "PaleoClim")

## Occurrence Collection and Management-----------------------------------------
# megaSDM can download occurrences from GBIF if desired:
# This function only takes occurrences from the described trainingarea extent.
# The defined extent should be the same (or similar to) as the extent of the training area.
# Given in latitude/longitude coordinates:
extent_occ <- c(-91.5, -75, 25.5, 36)

#Dfeine the species we want ot look at!
spplist <- c("Podomys floridanus")

# Define the file folder where the occurrences will be written, within the working directory
# (if this folder doesn't already exist, megaSDM will make it)
occ_output <- "occurrences"

#Collect the occurrences:
Occurrences <- OccurrenceCollection(spplist = spplist,
                                    output = occ_output,
                                    trainingarea = extent_occ)


# NOTE: when running this using R Markdown, you may get "incomplete final line..."
#    warnings. However, they do not appear to affect the total number or identity
#    of the occurrence points and when the code is run off of the console, the
#    warnings do not appear.

# Because one species was renamed, rename species list to reflect taxonomy changes
spplist <- Occurrences$Scientific.Name

# The OccurrenceManagement function will ensure that the occurrence data are
# properly formatted and extract the enviromental data for each occurrence point.
# In addition, if requested the function environmentally subsampled the occurrence
#data to reduce environmental bias.

# First, get the list of the occurrence files
occlist <- list.files(occ_output, pattern = ".csv", full.names = TRUE)

OccurrenceManagement(occlist = occlist,
                     output = occ_output,
                     envextract = TRUE,
                     envsample = TRUE,
                     nbins = 25,
                     envdata = TSEnv$training)

## Background Point Creation----------------------------------------------------
# Although the best method of generating background points is still up for debate,
# spatially-constrained background points can be more effective than randomly-
# generated background points for SDMs. Therefore, many people choose to
#randomly sample background points within a buffer of the occurrences.

#To do thisin megaSDM, we need to get the list of occurrence files again, even
#if they were written out in the same folder as before. This ensures that the
#occurrence files are properly formatted.
occlist <- list.files(occ_output, pattern = ".csv", full.names = TRUE)

# The location to print out the background buffers (.shp) (will be created if it doesn't exist)
buff_output <- "TestRun/buffers"

# Generates buffers for each species. buff_distance argument tells how wide the
#buffer should be (if set to NA, the distances between the occurrence points
#themselves govern the buffer width).
BackgroundBuffers(occlist = occlist,
                  envdata = TSEnv$training,
                  buff_output,
                  ncores = 2,
                  buff_distance = NA)

# Set the parameters for the background point generation
# (how many points, and how spatially-constrained)

# How many background points should be generated per species?
nbg <- 1000

# What proportion of the background points should be sampled from within the buffers?
spatial_weights <- 0.5

# Should the background points be environmentally subsampled (Varela) or
# randomly distributed (random)?
sampleMethod <- "Varela"

# Because we want a partial spatial constraint (50% of points within the buffer), we must make a
# list of the buffer files to use in the creation of the background points. In the example,
# these files are created from the BackgroundBuffers function, but they can also be generated
# outside of megaSDM and brought in here.

bufflist <- list.files(buff_output, pattern = ".shp$", full.names = TRUE)

# Define the location where the background points will be printed out to (.shp)
# (This directory will be created if it doesn't already exist)
bg_output <- "TestRun/backgrounds"

BackgroundPoints(spplist = spplist,
                 envdata = TSEnv$training,
                 output = bg_output,
                 nbg = nbg,
                 spatial_weights = spatial_weights,
                 buffers = bufflist,
                 method = sampleMethod,
                 ncores = 1)

## Environmental Variable Selection---------------------------------------------
# Define a list of the environmental variables to keep for each species
# In this example, we simply want all of the species to have the same environmental variables.
envvar <- rep("bio_1,bio_12,bio_14,bio_6,bio_9", length = length(occlist))

# Define a list of the background point files
# (either created in the BackgroundPoints function or generated separately)
bglist <- list.files(bg_output, pattern = ".csv", full.names = TRUE)

# In this example, megaSDM overwrites the occurrence and background points,
# but they could be placed in a different folder if requested.
VariableEnv(occlist = occlist,
            bglist = bglist,
            env_vars = envvar,
            occ_output = occ_output,
            bg_output = bg_output)

## MaxEnt Modelling and Projection----------------------------------------------
# First, define a list of all background and occurrence point files
occlist <- list.files(occ_output, pattern = ".csv", full.names = TRUE)
bglist <- list.files(bg_output, pattern = ".csv", full.names = TRUE)

# Define where the results of the MaxEnt model runs will be printed out to (as .lambdas files)
model_output <- "TestRun/models"

# "nrep" is set to 4, meaning that the MaxEnt algorithm will run 4 times with different
# subsets of occurrence points for a better representation of the habitat suitability.
MaxEntModel(occlist = occlist,
            bglist = bglist,
            model_output = model_output,
            ncores = 2,
            nrep = 4,
            alloutputs = FALSE)

# First, create a list of the time periods and climate scenarios used in the analysis
# (starting with the year the model is trained on)
time_periods <- c(2000, -4000, -11000, -15000, -21000)
scenarios <- c("PaleoClim")

# Define the directory where the current study area rasters are located
#    (generated from the TrainStudyEnv function or brought in from a separate location)
study_dir <- "TestRun/studyarea"

# Define the directories where the future study area rasters are location
#    (generated from the PredictEnv function or brought in from a separate location)


# Define a list of directories for the projected climate layers,
# separated into the different climate scenarios and years:
#    list(c(Scenario1Year1, Scenario1Year2),
#         c(Scenario2Year1, Scenario2Year2))

predictdir <- list(c("TestRun/PaleoClim/-4000",
                     "TestRun/PaleoClim/-11000",
                     "TestRun/PaleoClim/-15000",
                     "TestRun/PaleoClim/-21000"))

# Define Where the results will be printed out.
# For this example, We'll define a new folder within
# the working directory that is specifically for the model
# projecions and analysis.

result_dir <- "Results"

# Other options are also available (check the documentation page)
MaxEntProj(input = model_output,
           time_periods = time_periods,
           scenarios = scenarios,
           study_dir = study_dir,
           predict_dirs = predictdir,
           output = result_dir,
           aucval = 0.7,
           ncores = 1)

## Create Maps of Range Shifts--------------------------------------------------
# The time maps will be written out to the directory supplied in "result_dir"
 result_dir <- "Results"

 createTimeMaps(result_dir = result_dir,
                time_periods = time_periods,
                scenarios = scenarios,
                dispersal = FALSE,
                ncores = 1)

## Create Plots of range size change--------------------------------------------
 additionalStats(result_dir = result_dir,
                 time_periods = time_periods,
                 scenarios = scenarios,
                 dispersal = FALSE,
                 ncores = 1)

## Dispersal Data---------------------------------------------------------------
 # Add in dispersal data (normally you would read a .csv file with two columns, but in this example
 # the data is just added in by hand here).
 dispersaldata <- data.frame(Species = spplist, Rate = c(0.035))

 dispersalRate(result_dir = result_dir,
               dispersaldata = dispersaldata,
               time_periods = time_periods,
               scenarios = scenarios,
               hindcast = TRUE,
               startpoint = c("PaleoClim", "-21000"),
               ncores = 1)

 # Repeat the time map and additional stats steps for the dispersal constrained data.
 # Set dispersal = TRUE this time.
 createTimeMaps(result_dir,
                time_periods,
                scenarios,
                dispersal = TRUE,
                dispersaldata = dispersaldata,
                ncores = 1)

 # The additional stats function will compare the species ranges between the
 # dispersal-constrained and the regular data.
 additionalStats(result_dir,
                 time_periods,
                 scenarios,
                 dispersal = TRUE,
                 dispersaldata = dispersaldata,
                 ncores = 1)
