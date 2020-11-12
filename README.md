# LiLab_fitBMandOUmodels
# PARTI codes to fit BM models 
## 1. System requirements
- PBS (Protable Batch System)
- centOS 7
- python 2.7
- R 3.5
- R pacakage phytools 0.7-70

## 2. Installation guide

No installation needed. Directly run runall_BM.py

## 3. Demo/BMmodel

### Input:
- cat.tre
- forselection_Subcutaneous_adipose_mean.csv
- forselection_Subcutaneous_adipose_sd.csv

### Output:
- Subcutaneous_adipose-cat.BMresults.txt

## 4. Instructions for use

Change the hard-coded paths in the code to your own paths where you place the input and output files

# PARTII codes to fit OU model
## 1. System requirements
- PBS (Protable Batch System)
- centOS 7
- python 2.7
- R 3.5
- R pacakage phytools 0.7-70, geiger 2.0.7

## 2. Installation guide
No installation needed
Directly run runall_OU.py

## 3. Demo/OUmodel
### Input:
- cat.tre
- forselection_Subcutaneous_adipose_mean_*.csv
- forselection_Subcutaneous_adipose_sd_*.csv

### Output:
- Subcutaneous_adipose-cat-*.OUresults.txt

### 4. Instructions for use
Change the hard-coded paths in the code to your own paths where you place the input and output files

