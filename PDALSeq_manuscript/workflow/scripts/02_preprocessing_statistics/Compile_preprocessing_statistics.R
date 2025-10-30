
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

calcReads = function(x = "data/02_preprocessing_statistics/Merged/mapping_quality/{mapq}/{data}.txt"){
    df = read.delim(x, header = F, sep = "") 
    return(df)
}

calcFLAG = function(x = "data/02_preprocessing_statistics/Merged/bam_flags/{mapq}/{data}.txt"){
    df = read.delim(x, header = F, sep = "") %>% arrange(V2)
    unmapped = df$V1[which(df$V2 %in% c(4))]
    if (length(unmapped) == 0){unmapped = 0}
    mapped = sum(df$V1[which(df$V2 %in% c(0,16))])
    supplementary = sum(df$V1[which(df$V2 %in% c(2048,2064))])
    return(data.frame(unmapped,mapped,supplementary ))
}

df.index = read.delim(args[1], sep = "\t")

Raw.fastq.million.reads = c()
Trimmed.fastq.million.reads = c()
Mapped.bam.primary.million.reads = c()
Mapped.bam.supplementary.million.reads = c()
Mapped.bam.unmapped.million.reads = c()
Mapped.deduplicated.bam.primary.million.reads = c()
Mapped.deduplicated.bam.supplementary.million.reads = c()
Mapped.deduplicated.bam.unmapped.million.reads = c()
Mapped.deduplicated.high.MapQ.bam.primary.million.reads = c()
Mapped.deduplicated.high.MapQ.bam.supplementary.million.reads = c()
Mapped.deduplicated.high.MapQ.bam.unmapped.million.reads = c()

for (i in 1:nrow(df.index)){
    x = df.index$Sequencing_sample[i]

    df.raw.fastq = calcReads(paste("data/02_preprocessing_statistics/Raw_fastq/", x, ".txt", sep = ""))
    Raw.fastq.million.reads[i] = df.raw.fastq$V1[1]
    
    df.Trimmed.fastq = calcReads(paste("data/02_preprocessing_statistics/Trimmed_fastq/", x, ".txt", sep = ""))
    Trimmed.fastq.million.reads[i] = df.Trimmed.fastq$V1[1]
    
    df.mapped.bam = calcFLAG(paste("data/02_preprocessing_statistics/Mapped_bam/", x, ".txt", sep = ""))
    Mapped.bam.primary.million.reads[i] = df.mapped.bam$mapped[1]
    Mapped.bam.supplementary.million.reads[i] = df.mapped.bam$supplementary[1]
    Mapped.bam.unmapped.million.reads[i] = df.mapped.bam$unmapped[1]
    
    df.mapped.deduplicated.bam = calcFLAG(paste("data/02_preprocessing_statistics/Mapped_deduplicated_bam/", x, ".txt", sep = ""))
    Mapped.deduplicated.bam.primary.million.reads[i] = df.mapped.deduplicated.bam$mapped[1]
    Mapped.deduplicated.bam.supplementary.million.reads[i] = df.mapped.deduplicated.bam$supplementary[1]
    Mapped.deduplicated.bam.unmapped.million.reads[i] = df.mapped.deduplicated.bam$unmapped[1]
    
    df.deduplicated.high.MapQ.bam = calcFLAG(paste("data/02_preprocessing_statistics/Mapped_deduplicated_high_MapQ_bam/", x, ".txt", sep = ""))
    Mapped.deduplicated.high.MapQ.bam.primary.million.reads[i] = df.deduplicated.high.MapQ.bam$mapped[1]
    Mapped.deduplicated.high.MapQ.bam.supplementary.million.reads[i] = df.deduplicated.high.MapQ.bam$supplementary[1]
    Mapped.deduplicated.high.MapQ.bam.unmapped.million.reads[i] = df.deduplicated.high.MapQ.bam$unmapped[1]
}

df = df.index[,1:5]
df$Raw.fastq.million.reads = Raw.fastq.million.reads
df$Trimmed.fastq.million.reads = Trimmed.fastq.million.reads
df$Mapped.bam.primary.million.reads = Mapped.bam.primary.million.reads
df$Mapped.bam.supplementary.million.reads = Mapped.bam.supplementary.million.reads
df$Mapped.bam.unmapped.million.reads = Mapped.bam.unmapped.million.reads
df$Mapped.deduplicated.bam.primary.million.reads = Mapped.deduplicated.bam.primary.million.reads
df$Mapped.deduplicated.bam.supplementary.million.reads =  Mapped.deduplicated.bam.supplementary.million.reads
df$Mapped.deduplicated.bam.unmapped.million.reads = Mapped.deduplicated.bam.unmapped.million.reads
df$Mapped.deduplicated.high.MapQ.bam.primary.million.reads = Mapped.deduplicated.high.MapQ.bam.primary.million.reads
df$Mapped.deduplicated.high.MapQ.bam.supplementary.million.reads = Mapped.deduplicated.high.MapQ.bam.supplementary.million.reads
df$Mapped.deduplicated.high.MapQ.bam.unmapped.million.reads = Mapped.deduplicated.high.MapQ.bam.unmapped.million.reads

write.csv(df, args[2], row.names = F)
