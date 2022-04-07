
##### ANNOTATE (by overlap) a BED File ######
# Annotate and write any bed file that overlap with :
# - CGI Annotated with the closest promoter
# - Promoter (-1500 +1000), you can chage the ranges
# - Genebody
#
# Genes are given as Symbol, ENSEMBL ID, ENtrez ID and one refseq


##################################################################################################################################

# require R > 4.0

# Set your working directory with the following files:
# -your bed file
# -Annot_Oli.RData
# -Annot_Oli.R

# run the following commands
# Load the RData
load("Annot_Oli.RData")
library(GenomicRanges)


# load your bed file (a tab delimited file with a least 3 columns "chr" "start" "end")
# if you have an header
mybed <- read.delim("Yourdatatoannotate.bed", header=TRUE)

# if you don't have an header
mybed <- read.delim("Yourdatatoannotate.bed", header=FALSE)
colnames(mybed)=c("chrom","Start","End") # you can complete with all col names



# set a GRanges object
mybed <- as(mybed,"GRanges")

# Set promoter ranges
#Upstream
UP <- 1500
#Downstream
DWN <- 1000


# Set your output basename, No space please!
OUTPUT <- "my_ChIPseq_Peaks"

# Annotate and write outputs, run with this commands
AnnotCGI(mybed)
AnnotProm(mybed)
AnnotGB(mybed)





#########
#Bed for test. run this to test the package, you will get your outputs
# remove them after from your dir
AnnotCGI(mybedtest)
AnnotProm(mybedtest)
AnnotGB(mybedtest)
#########


##################################################################################################################################
# Function used in the RData
# Don't touch this
##################################################################################################################################

AnnotCGI <- function(mybed){                                                          # 3 args
  annot <- mm10CGIprom.gr
  overlap <- subsetByOverlaps(mybed, annot)                                     # Subset Diffmeth by overlap with annot
  dist2annot <- distanceToNearest(overlap, annot)                                 # Calculate distance between features
  exportdiffannot <- data.frame(overlap, Overlap2_ = annot[dist2annot@to])        # Export annoted Diff meth
  write.table(exportdiffannot, 
              paste0(OUTPUT, "_overlap_CGI.txt"),
              sep='\t',
              dec=",",
              col.names=T,
              row.names=F,
              quote=F) 
}



AnnotProm <- function(mybed){                                                          # 3 args
  prom.gr = GenomicRanges::promoters(mm10Genes.gr,upstream=UP, downstream=DWN, use.names=F)
  prom.gr
  annot <- prom.gr
  overlap <- subsetByOverlaps(mybed, annot)                                     # Subset Diffmeth by overlap with annot
  dist2annot <- distanceToNearest(overlap, annot)                                 # Calculate distance between features
  exportdiffannot <- data.frame(overlap, Overlap2_ = annot[dist2annot@to])        # Export annoted Diff meth
  write.table(exportdiffannot, 
              paste0(OUTPUT, "_overlap_Prom.txt"),
              sep='\t',
              dec=",",
              col.names=T,
              row.names=F,
              quote=F) 
}


AnnotGB <- function(mybed){                                                          # 3 args
  annot <- mm10Genes.gr
  overlap <- subsetByOverlaps(mybed, annot)                                     # Subset Diffmeth by overlap with annot
  dist2annot <- distanceToNearest(overlap, annot)                                 # Calculate distance between features
  exportdiffannot <- data.frame(overlap, Overlap2_ = annot[dist2annot@to])        # Export annoted Diff meth
  write.table(exportdiffannot, 
              paste0(OUTPUT, "_overlap_Genebody.txt"),
              sep='\t',
              dec=",",
              col.names=T,
              row.names=F,
              quote=F) 
}



##################################################################################################################################



