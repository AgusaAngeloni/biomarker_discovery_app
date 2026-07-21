############################################################
# Script name: generate_pb_data.R
# Description: Generate a genome-wide DNA methylation 
# beta-value matrix from PB IDAT files using the 
# sesame preprocessing framework.
#
# Raw IDAT files were downloaded from GEO and correspond
# to peripheral blood PB samples profiled using
# the Illumina Infinium MethylationEPIC array.
#
# GEO Series:
# - GSE247869
# - GSE270856
#
# Prerequisites:
#  Download and extract all IDAT files from
#  GEO https://www.ncbi.nlm.nih.gov/geo/
#  accessions GSE247869 and GSE270856
#  into:
#      data/geo/IDATS/
#  The directory should contain paired
#  		*_Red.idat and *_Grn.idat files for all samples.
#
# Processing workflow:
# 1. Read all IDAT files from the input directory.
# 2. Perform quality control and preprocessing using
# sesame ("QCDPB" pipeline).
# 3. Extract beta values for all CpG sites.
# 4. Build a methylation matrix with CpG probes as rows
# and samples as columns.
# 5. Export the resulting matrix as a parquet file.
#
# Input:
# data/geo/IDATS/
#
# Output:
# data/methy_pb.parquet
#
# Requirements:
# - R >= 4.3
# - sesame
# - dplyr
# - arrow
#
############################################################

# ---- Install and load required packages ----

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sesame")

library(sesame)
library(dplyr)
library(arrow)

# ---- Input directory containing IDAT files ----

idat_pb <- "data/geo/IDATS"

# ---- Read and preprocess methylation data ----

betasAllpb <- openSesame(
  idat_pb,
  prep = "QCDPB"
)

# ---- Convert to data frame ----

betasAllpbdf <- as.data.frame(
  betasAllpb
)

# Add probe identifiers

betasAllpbdf$site_id <- rownames(
  betasAllpb
)

# ---- Export parquet file ----

write_parquet(
  betasAllpbf,
  "data/raw/methy_pb.parquet"
)
