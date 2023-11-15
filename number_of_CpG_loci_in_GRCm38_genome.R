# How to find CpG sites in mouse genome GRCm38?
# important link for this pipeline: https://stackoverflow.com/questions/33648447/mapping-cpg-coordinate-using-a-bam-file

# Start

# Loading related packages

library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)

#Create the chromossome set

# Find chromosome names in mouse genome

chr_all = paste("chr", c(1:19, "X", "Y", "M"), sep = "")

#CpG Mapping on Genome
CpGpos<-lapply(chr_all,function(x) { 
  RES<-matchPattern("CG",Mmusculus[[x]])
  RES<-cbind(rep(x,length(RES)),RES@ranges@start)
  return(RES) 
})

CpGpos<-do.call("rbind",CpGpos)
CpGstartend<-cbind(CpGpos,as.numeric(CpGpos[,2])+1)
CpGstartend<-as.data.frame(CpGstartend,stringsAsFactors=F)
CpGstartend$V1 <- mgsub::mgsub(CpGstartend$V1, "chr", "")
write.table(CpGstartend, "/path_to_directory/CPGstartend_GRCm38_mouse_genome.txt", sep="\t", col.names=F, row.names=T, quote = F)
str(CpGstartend)

# Another pipeline: https://support.bioconductor.org/p/95239/

names(Mmusculus)[1:22]
cgs <- lapply(chr_all, function(x) start(matchPattern("CG", Mmusculus[[x]])))
cpgr <- do.call(c, lapply(1:22, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs[[x]], width = 2))))
str(cpgr)

# Both pipelines give the same results in terms of number of CpGs and coordinates of CpGs in mouse genome.

# End
