library(shiny)
library(bslib)
library(bsplus)
library(reactable)
library(waiter)

library(httr)
library(htmltools)
library(XML)

library(RColorBrewer)
library(viridis)
library(ggiraph)
library(png)

library(tidyverse)

# load supporting file (mostly in-house functions)
source("shiny_multiomic_web_source.R")

# Load Necessary Data -----
# result of metabolomic estimation result from each tissue (CellMet 2025)
dat_metab <- readRDS("data/CellMet2025_metabolites_adjustedLFC_2025-11-25.rds")
# compound IDs annotations
dat_compound_id <- read.csv("data/compound_mapping_res_2025-11-25.csv")

# load KEGG pathway IDs
kegg_ids <- readRDS("data/KEGG_pathway_ID_list_2025_12_17_exported.rds")


# Define Objects Required for Dynamic UI Change -----
# change tissue list in response to datasource (RNA/Protein)
choices_tissue <- list(
  rna = c("BAT", "Bone", "Brain", "GAT", "Heart", "Kidney", "Limb_Muscle", "Liver", "Lung", "Marrow", "MAT",
          "Pancreas", "SCAT", "Skin", "Small_Intestine", "Spleen", "WBC"),
  protein = c("kidney", "hippocampus", "heart", "fat", "cerebellum", "striatum", "spleen", "skeletalmuscle", "lung", "liver"),
  not_show = NULL
)

# change analytical design based on the data source (RNA/Protein)
choices_design <- list(
  rna = c("linear model using all age groups" = "linear",
          "two age groups comparison" = "two_group"),
  protein = c("two age groups comparison" = "two_group")
)

# Age group can be chosen in bulkRNA dataset
# Web-ver only: only 3-month vs 21-month comparison is available

