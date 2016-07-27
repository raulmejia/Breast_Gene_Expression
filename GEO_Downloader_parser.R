#install GEOquery
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description="Collapse the expression matrix by gene symbols with B statistical")
parser$add_argument("--matrix", required=TRUE, help="normalized expression matrix")
args <- parser$parse_args()
print(args$matrix)
edata <- read.table(file=args$matrix, header=FALSE, sep="\t")
print(edata[,1])

############################################################################  
# A) If you have the GSM vector try to Download the .CEL files through:  ###
############################################################################

# You need your own GSM vector, or for this time you can use the following:
############### Example GSM vector #################
#MyGSMList = c("GSM272923","GSM272924","GSM272925")##
####################################################

sapply(as.character(edata[,1]),function(x){ getGEOSuppFiles(x)})


