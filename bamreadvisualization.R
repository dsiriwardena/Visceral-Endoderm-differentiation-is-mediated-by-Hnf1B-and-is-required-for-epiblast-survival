if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("karyoploteR")
BiocManager::install("Rsamtools")
BiocManager::install("BSgenome")
BiocManager::install("Rsubread")

#load library
library(Rsamtools)
library(karyoploteR)
library(BSgenome)
library(GenomicAlignments)
library(rtracklayer)
library(Rsubread)
setwd("C:/Users/Dylan/Documents/cambridge/forDylan")
samples<-read.table("visceralendoderm_cells.csv",sep=",",header=T)
for (i in 1:nrow(samples)){
  temp_name<-toString(samples[i,1])
  temp_shortn<-toString(samples[i,2])
  temp_emb<-toString(samples[i,3])
  bam <- scanBam(temp_name)
  bai<-indexBam(temp_name)
  jpeg( paste("C:/",temp_emb,"/",temp_shortn,"exon1.jpg",sep = ""), width = 900, height = 900)
    kp <- plotKaryotype(genome = "mm10", chromosomes = "chr11",zoom=toGRanges("chr11:83851049-83851371"))
    kpAddBaseNumbers(kp, tick.dist = 5000, add.units = TRUE)
    
    kp <- kpPlotBAMCoverage(kp, data=temp_name)
    kp <- kpPlotBAMCoverage(kp, data=temp_name, max.valid.region.size = 2e6)
    kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage)
    #https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotBAMCoverage/PlotBAMCoverage.html
  dev.off()
  
  jpeg( paste("C:/exon2/",temp_emb,"/",temp_shortn,"exon2.jpg",sep = ""), width = 900, height = 900)
  kp <- plotKaryotype(genome = "mm10", chromosomes = "chr11",zoom=toGRanges("chr11:83855904-83856103"))
  kpAddBaseNumbers(kp, tick.dist = 50, add.units = TRUE)
  
  kp <- kpPlotBAMCoverage(kp, data=temp_name)
  kp <- kpPlotBAMCoverage(kp, data=temp_name, max.valid.region.size = 2e6)
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage)
  #https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotBAMCoverage/PlotBAMCoverage.html
  dev.off()
}

