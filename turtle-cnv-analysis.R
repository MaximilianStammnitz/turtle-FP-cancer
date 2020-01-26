####################################################
#        Sea Turtle FP Tumour CNV Analysis         #
#                                                  #
# 26.01.2020                                       #
# TCG Cambridge, Department of Veterinary Medicine #
# mrs72@cam.ac.uk                                  #
####################################################

# Libraries and Dependencies
library(cn.mops)
library(DNAcopy)
library(GenomicRanges)

setwd("/Users/ms37/Desktop/Labwork/Side Projects/Turtles/src")


# 1. Input 5000 bp read counts from three FP tumours
####################################################

counts.Lil.5000bp <- read.table('../liRRF4dnavliCSVS2dna_readcounts_5000_BP.txt', header = T, strip.white = T, check.names = F, sep = '\t')
counts.27L.5000bp <- read.table('../27L1Fdnav27L5Hdna_readcounts_5000_BP.txt', header = T, strip.white = T, check.names = F, sep = '\t')
counts.27K.5000bp <- read.table('../27K2Fdnav27K4Hdna_readcounts_5000_BP.txt', header = T, strip.white = T, check.names = F, sep = '\t')
colnames(counts.Lil.5000bp)[4:5] <- colnames(counts.27L.5000bp)[4:5] <- colnames(counts.27K.5000bp)[4:5] <- c('TUMOUR', 'HOST')


# 2. Concatenate contigs into one big chromosome and convert counts into GR object
##################################################################################

countMatrix.to.countRange <- function(x, win.size){

  # a. 'Uniformize' chromosome locations
  x <- x[c(x[,3] - x[,2]) == c(win.size - 1),]
  x.contigs <- unique(as.character(x[,1]))
  for (i in 2:length(x.contigs)){
    
    x.contigs.tmp <- x.contigs[i]
    x.contigs.tmp.ind <- which(x[,1] == x.contigs.tmp)
    x[x.contigs.tmp.ind,2:3] <- x[x.contigs.tmp.ind,2:3] + x[tail(which(x[,1] == x.contigs[i-1]), 1),3]
    
  }

  # b. GRanges reconversion
  range_1 <- '1'
  range_2 <- IRanges::IRanges(start = x[,2], end = x[,3])
  range_counts <- as.matrix(x[,-c(1:3)])
  rownames(range_counts) <- NULL; mode(range_counts) <- "integer"
  GR <- GenomicRanges::GRanges(seqnames = range_1, ranges = range_2)
  IRanges::values(GR) <- range_counts
  GR <- sortSeqlevels(GR)
  
  # c. Output
  return(GR)
}
counts.Lil.5000bp.GR <- countMatrix.to.countRange(counts.Lil.5000bp, 5000)
counts.27L.5000bp.GR <- countMatrix.to.countRange(counts.27L.5000bp, 5000)
counts.27K.5000bp.GR <- countMatrix.to.countRange(counts.27K.5000bp, 5000)


# 3. Normalise coverages
########################

counts.Lil.5000bp.GR.norm <- normalizeGenome(counts.Lil.5000bp.GR, normType = "mode")
counts.27L.5000bp.GR.norm <- normalizeGenome(counts.27L.5000bp.GR, normType = "mode")
counts.27K.5000bp.GR.norm <- normalizeGenome(counts.27K.5000bp.GR, normType = "mode")


# 4. Segment readcounts of the three tumours vs. their matched host tissue
##########################################################################

counts.Lil.5000bp.segmented <- referencecn.mops(cases = counts.Lil.5000bp.GR.norm[,'TUMOUR'],
                                                controls = counts.Lil.5000bp.GR.norm[,'HOST'],
                                                norm = 0,
                                                minReadCount = 3,
                                                segAlgorithm = 'fast',
                                                I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 8, 16, 32, 64), 
                                                classes = c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6", 
                                                            "CN7","CN8","CN16","CN32","CN64","CN128"))
counts.Lil.5000bp.segmented <- calcIntegerCopyNumbers(counts.Lil.5000bp.segmented)

counts.27L.5000bp.segmented <- referencecn.mops(cases = counts.27L.5000bp.GR.norm[,'TUMOUR'],
                                                controls = counts.27L.5000bp.GR.norm[,'HOST'],
                                                norm = 0,
                                                minReadCount = 3,
                                                segAlgorithm = 'fast',
                                                I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 8, 16, 32, 64), 
                                                classes = c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6", 
                                                            "CN7","CN8","CN16","CN32","CN64","CN128"),
                                                returnPosterior = T)
counts.27L.5000bp.segmented <- calcIntegerCopyNumbers(counts.27L.5000bp.segmented)

counts.27K.5000bp.segmented <- referencecn.mops(cases = counts.27K.5000bp.GR.norm[,'TUMOUR'],
                                                controls = counts.27K.5000bp.GR.norm[,'HOST'],
                                                norm = 0,
                                                minReadCount = 3,
                                                segAlgorithm = 'fast',
                                                I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 8, 16, 32, 64), 
                                                classes = c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6", 
                                                            "CN7","CN8","CN16","CN32","CN64","CN128"))
counts.27K.5000bp.segmented <- calcIntegerCopyNumbers(counts.27K.5000bp.segmented)


# 5. Segmentation plots
#######################

nice.seg.plot <- function(x.segmented, x.normalised.counts, title){
  
  ratio <- elementMetadata(x.normalised.counts)[,'TUMOUR']/elementMetadata(x.normalised.counts)[,'HOST']
  ratio <- ratio*2
  
  pdf(paste0('../results/', title, '.pdf'),
      width = 14, height = 7)
  mar.default <- c(5, 4, 4, 2) + 0.1
  par(mar = mar.default + c(0, 3, 0, 0))
  
  # initialise plot with (normalised) read-count ratios: tumour/host across window
  plot(x = as.matrix(ranges(x.normalised.counts))[,1], y = ratio, 
       type = 'p', xlab = 'Genome Position', ylab = 'Copy Number', pch = 16, cex = 0.1,
       col = 'grey', ylim = c(0, 4), main = title, cex.lab = 2, cex.main = 3)
  
  # add segmentation lines in red
  segs <- as.data.frame(segmentation(x.segmented))
  segs[,'CN'] <- as.numeric(str_split_fixed(segs[,'CN'], 'N', 2)[,2])
  for (i in 1:nrow(segs)){
    
    lines(x = c(segs[i,'start'], segs[i,'end']), y = rep(segs[i,'CN'],2), col = 'red', lwd = 2)
    
  }
  
  dev.off()
  
}

nice.seg.plot(x.segmented = counts.Lil.5000bp.segmented,
              x.normalised.counts = counts.Lil.5000bp.GR.norm,
              title = 'Lilac FP - segmented CNVs')

nice.seg.plot(x.segmented = counts.27L.5000bp.segmented,
              x.normalised.counts = counts.27L.5000bp.GR.norm,
              title = '27L - segmented CNVs')

nice.seg.plot(x.segmented = counts.27K.5000bp.segmented,
              x.normalised.counts = counts.27K.5000bp.GR.norm,
              title = '27K - segmented CNVs')

# 6. Re-annotate the old contig names
#####################################

counts.Lil.5000bp <- counts.Lil.5000bp[c(counts.Lil.5000bp[,3] - counts.Lil.5000bp[,2]) == c(5000 - 1),]
counts.27L.5000bp <- counts.27L.5000bp[c(counts.27L.5000bp[,3] - counts.27L.5000bp[,2]) == c(5000 - 1),]
counts.27K.5000bp <- counts.27K.5000bp[c(counts.27K.5000bp[,3] - counts.27K.5000bp[,2]) == c(5000 - 1),]

segs.counts.Lil.5000bp <- as.data.frame(integerCopyNumber(counts.Lil.5000bp.segmented))
segs.counts.Lil.5000bp <- cbind(segs.counts.Lil.5000bp[,1:3], counts.Lil.5000bp[,1:3],
                                elementMetadata(counts.Lil.5000bp.GR.norm), segs.counts.Lil.5000bp[,6])

segs.counts.27L.5000bp <- as.data.frame(integerCopyNumber(counts.27L.5000bp.segmented))
segs.counts.27L.5000bp <- cbind(segs.counts.27L.5000bp[,1:3], counts.27L.5000bp[,1:3],
                                elementMetadata(counts.27L.5000bp.GR.norm), segs.counts.27L.5000bp[,6])

segs.counts.27K.5000bp <- as.data.frame(integerCopyNumber(counts.27K.5000bp.segmented))
segs.counts.27K.5000bp <- cbind(segs.counts.27K.5000bp[,1:3], counts.27K.5000bp[,1:3],
                                elementMetadata(counts.27K.5000bp.GR.norm), segs.counts.27K.5000bp[,6])

colnames(segs.counts.27K.5000bp) <- colnames(segs.counts.27L.5000bp) <- colnames(segs.counts.Lil.5000bp) <- 
  c('CHR', 'CHR_START', 'CHR_END', 'SCAFFOLD', 'SCAFFOLD_START', 'SCAFFOLD_END', 'TUMOUR_COUNTS_norm', 
    'HOST_COUNTS_norm', 'COPY_NUMBER_pred')