
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

MAPQ = c("All_reads", "High_MapQ")

calcQMAP = function(x = "data/02_preprocessing_statistics/Merged/mapping_quality/{mapq}/{data}.txt"){
    df = read.delim(x, header = F, sep = "") %>% arrange(V2)
    TOTAL = sum(df$V1)
    less.than.equal.to.Q10 = sum(df$V1[which(df$V2 <= 10)])
    greater.than.Q10 = sum(df$V1[which(df$V2 >= 10)])
    greater.than.Q20 = sum(df$V1[which(df$V2 >= 20)])
    greater.than.Q30 = sum(df$V1[which(df$V2 >= 30)])
    return(data.frame(less.than.equal.to.Q10, greater.than.Q10, greater.than.Q20, greater.than.Q30))
}

calcFLAG = function(x = "data/02_preprocessing_statistics/Merged/bam_flags/{mapq}/{data}.txt"){
    df = read.delim(x, header = F, sep = "") %>% arrange(V2)
    TOTAL = sum(df$V1[which(df$V2 %in% c(0,4,16))])
    if(length(which(df$V2 %in% c(4))) > 0){unmapped = df$V1[which(df$V2 %in% c(4))]}else{unmapped = 0} 
    mapped = sum(df$V1[which(df$V2 %in% c(0,16))])
    secondary = sum(df$V1[which(df$V2 %in% c(2048,2064))])
    return(data.frame(TOTAL,unmapped,mapped,secondary))
}

calcUMR = function(x = "data/02_preprocessing_statistics/Merged/primary_versus_supplemental/{mapq}/{data}.txt"){
    df = read.delim(x, header = T,sep = "")
    TOTAL = sum(df$Primary.score.minus.suplemental.score)
    HighMQ = which(df$Primary.score.minus.suplemental.score >= 10)
    if (length(HighMQ) > 0){
    	uniquely.mapped = sum(df$count[HighMQ])
    }else{
    	uniquely.mapped = sum(df$count)
    }
    return(data.frame(uniquely.mapped))
}


df.index = read.delim(args[1], sep = "\t")

df.i = data.frame(expand.grid(unique(df.index$Experiment), MAPQ))
colnames(df.i) = c("Experiment", "MapQ")

l.df.out = {}

for (i in 1:nrow(df.i)){
    print(i)
    df.FLAG = calcFLAG(paste("data/02_preprocessing_statistics/Merged/bam_flags/", df.i$MapQ[i], "/", df.i$Experiment[i], ".txt", sep = ""))
    df.Quality = calcQMAP(paste("data/02_preprocessing_statistics/Merged/mapping_quality/", df.i$MapQ[i], "/", df.i$Experiment[i], ".txt", sep = ""))
    df.UMR = calcUMR(paste("data/02_preprocessing_statistics/Merged/primary_versus_supplemental/", df.i$MapQ[i], "/", df.i$Experiment[i], ".txt", sep = ""))
    l.df.out[[i]] = cbind(df.UMR, df.FLAG, df.Quality)
    l.df.out[[i]]$MapQ = df.i$MapQ[i]
    l.df.out[[i]]$Experiment = df.i$Experiment[i]
}

df.out = bind_rows(l.df.out)

write.csv(df.out, args[2], row.names = F)
