# Disentangling bias for non-destructive insect metabarcoding

This repository contains the code required for reproducing the analyses presented in our manuscript: _Martoni, F., Piper, A. M., Rodoni, B. C., & Blacket, M. J. (2022). Disentangling bias for non-destructive insect metabarcoding. PeerJ, 10, e12981._

## Data

The input sequencing data are not included in the repository for size reasons, and are instead available from the SRA under accession: [PRJNA767112](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA767112)

However RDS files holding intermediate data objects such as the OTU and taxonomy tables suitable for performing the analyses are contained within output/rds and these can be used to reproduce the statistical analyses.

## Analysis 

The bioinformatic workflow used to analyse the sequence data can be found in the bioinformatics.rmd file in the root directory or rendered: [Here](https://alexpiper.github.io/disentangling_metabarcoding_bias_MS/bioinformatics.html)

The statistical analyses can be found in the statistics.rmd file in the root directory or rendered: [Here](https://alexpiper.github.io/disentangling_metabarcoding_bias_MS/statistics.html)
