1.Overview of CRSCNV
CRSCNV is a cross-model based statistical approach to detect CNVs in individual samples
from next-generation sequencing data.In order to verify the validity of CRSCNV,we use it
to detect simulation data and real data compared with several existing methods,which 
proves to be a very effective and reliable tool.
2.Usage
2.1 Operating environment
Our software is running on windows system.It is developed using R and Python language.Users 
need to install the R and Python build environment before running the software.We need to 
install the R packages(readr and PythonInR) and Python packages(numpy and numba)
2.2 Input
The user opens the input configuration file(RD_filename), enters the RD file in the first 
column, and enters the installation path of Python in the second column.
2.3 Run command
library(PythonInR)
library(readr)
source("CRSCNV.R")
2.4 Output
The output is saved in the result file
3 Application
We provide an application example for the user to further use the software.