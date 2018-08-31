#!/usr/bin/Rscript
# This script can be used to analyze a batch of sequences without having to use
# GUI or R directly.
# It checks to see if the required tools are installed (devtools, mlplasmids)

if(!"devtools" %in% rownames(installed.packages())) {
  install.packages("devtools", repos='http://cran.us.r-project.org')
}

if(!"Biostrings" %in% rownames(installed.packages())) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("Biostrings")
}
if(!"mlplasmids" %in% rownames(installed.packages())) {
  devtools::install_git("https://gitlab.com/sirarredondo/mlplasmids",
                        repos='http://cran.us.r-project.org')
}
library(mlplasmids)

print("USAGE: Rscript run_mlplasmids.R ./path/to/assembly.fasta")
# ENABLE command line arguments
args <- commandArgs(TRUE)


# Change the following object with the system location of your input file
my_path <- (args[1])
if(is.na(args[2])
stop()
if (!is.args[2])
# example_prediction <- plasmid_classification(path_input_file = my_path, full_output = TRUE, prob_threshold=.8, species = "Escherichia coli")
example_prediction <- plasmid_classification(path_input_file = my_path,  prob_threshold=.8, species = "Escherichia coli")
write.table(x=example_prediction, file=args[2], row.names=F, sep="\t")
