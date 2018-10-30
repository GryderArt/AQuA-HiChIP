### APA plots and comparisons

### 1.Set up project and samples. Configure variables to adopt code into your own file system.
    setwd("K:/projects/ChIP_seq/projects/RMS_Epigenetics/RH4/")
    project.folder = "VASE/"
    sample.list = list.dirs(path = project.folder, full.names = F, recursive = F)
    sample.list = sample.list[grep("RH4",sample.list)]
    
### 2. Extract APA matrices
    APAmatrix.scale = 10000
    APA.all = as.data.frame(read.table(paste(project.folder,sample.list[1],"/",APAmatrix.scale,"/gw/APA.txt",sep=""), sep=",", header=F))
    APA.all$V1 = as.numeric(gsub("\\[|\\]", "", APA.all$V1)) ;APA.all$V21 = as.numeric(gsub("\\[|\\]", "", APA.all$V21))
    APA.all$comparison.name = sample.list[1]
    APA.all = APA.all[0,]
    
    lapply(sample.list, function(x) {
      ##load and merge sample data
      APA = as.data.frame(read.table(paste(project.folder,x,"/",APAmatrix.scale,"/gw/APA.txt",sep=""), sep=",", header=F))
      APA$V1 = as.numeric(gsub("\\[|\\]", "", APA$V1)) ;APA$V21 = as.numeric(gsub("\\[|\\]", "", APA$V21))
      APA$comparison.name = x
      APA.all <<- rbind(APA.all,APA)
    })
    APA.DMSO = APA.all[grep("D6", APA.all$comparison.name),]
    APA.Ent6 = APA.all[grep("Ent6", APA.all$comparison.name),]
    
          library(pheatmap)
          pheatmap(APA.all[,1:21], cluster_rows = F, cluster_cols = F)
          pheatmap(APA.DMSO[,1:21], cluster_rows = F, cluster_cols = F)
          pheatmap(APA.Ent6[,1:21], cluster_rows = F, cluster_cols = F)
    
### 3. Convert to contacts per million, then AQuA normalize
    
    sample1 = "Sample_RH4_D6_H3K27ac_HiChIP_HKJ22BGX7";sample2 = "Sample_RH4_Ent6_H3K27ac_HiChIP_HKJ22BGX7"
        removable.string = "Sample_" ; name1 = gsub(removable.string,"",sample1) ; name2 = gsub(removable.string,"",sample2) 
        mm10_name1 = gsub(removable.string,"mm10_",sample1) ; mm10_name2 = gsub(removable.string,"mm10_",sample2)
    mergestat.all = read.table("VASE/mergestat.HiChIP.all.txt", sep="\t", header = T)    #generated from plotVirtual4C.R
    
    APA.DMSO.CPM = APA.DMSO[,1:21]*1000000/(mergestat.all[2,name1]+mergestat.all[2,mm10_name1]) ; APA.Ent6.CPM = APA.Ent6[,1:21]*1000000/(mergestat.all[2,name2]+mergestat.all[2,mm10_name2])
    AQuAfactor1 = mergestat.all[2,name1]/mergestat.all[2,mm10_name1] ; AQuAfactor2 = mergestat.all[2,name2]/mergestat.all[2,mm10_name2]
    #AQuAfactor1 = 1; AQuAfactor2 = 1  #set AQUA to 1, just to compare
    APA.DMSO.RCPM = APA.DMSO.CPM*AQuAfactor1; APA.Ent6.RCPM = APA.Ent6.CPM*AQuAfactor2
    APA.DMSO.RCPM$comparison.name = APA.DMSO$comparison.name; APA.Ent6.RCPM$comparison.name = APA.Ent6$comparison.name
    APA.EvD.RCPM = APA.Ent6.RCPM[,1:21]-APA.DMSO.RCPM[,1:21]; APA.EvD.RCPM$comparison.name = APA.DMSO$comparison.name

          breakList = seq(0, 200000, by = 2000)
          pheatmap(subset(APA.DMSO.RCPM[,1:21], APA.DMSO.RCPM$comparison.name=="RH4_D6_H3K27ac_intraAPA"), cluster_rows = F, cluster_cols = F, breaks = breakList)
          pheatmap(subset(APA.Ent6.RCPM[,1:21], APA.DMSO.RCPM$comparison.name=="RH4_D6_H3K27ac_intraAPA"), cluster_rows = F, cluster_cols = F, breaks = breakList)
          pheatmap(subset(APA.EvD.RCPM[,1:21], APA.DMSO.RCPM$comparison.name=="RH4_D6_H3K27ac_intraAPA"), cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white", "mediumvioletred"))(33))
    
    ## plots that compare DMSO and Ent with the same scale set by APA MAX
          
    comparison.pick = "CRC_RH4_D6_H3K27ac_intraAPA"
    comparison.pick = "CRC_RH4_D6_H3K27ac_outsideVASE_APA"
    
    quant_cut = 0.95  #caps the contact map plot values at a given percentile
    APA.matrix.control = as.matrix(subset(APA.DMSO.RCPM[,1:21], APA.DMSO.RCPM$comparison.name==comparison.pick)); APA.matrix.treated = as.matrix(subset(APA.Ent6.RCPM[,1:21], APA.DMSO.RCPM$comparison.name==comparison.pick))
    quantile(APA.matrix.control, probs = c(quant_cut))
    
      APA.max = max(quantile(APA.matrix.control, probs = c(quant_cut)),quantile(APA.matrix.treated, probs = c(quant_cut)))
      #APA.max = 1200
      APA.max.control = APA.matrix.control; APA.max.control[APA.max.control>APA.max] <- APA.max
      APA.max.treated = APA.matrix.treated; APA.max.treated[APA.max.treated>APA.max] <- APA.max
      
      breakList = seq(0, APA.max, by = APA.max/100)
      breakList.delta = seq(-APA.max, APA.max, by = APA.max/50)
      pheatmap(APA.max.control, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white", "red"))(100), breaks = breakList)
      pheatmap(APA.max.treated, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white", "red"))(100), breaks = breakList)
      #pheatmap(log2(APA.max.treated)-log2(APA.max.control), cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("dodgerblue", "white", "mediumvioletred"))(60), breaks = seq(0,2.2,0.04))
      pheatmap(APA.max.treated-APA.max.control, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("dodgerblue", "white", "mediumvioletred"))(100),breaks = breakList.delta)
    
    ## cut rows to show an APA slice
      library(tidyr)
      
      APA.DMSO.RCPM.xslice= APA.DMSO.RCPM[seq(11, nrow(APA.DMSO.RCPM), 21), ]; APA.DMSO.RCPM.xlong = gather(APA.DMSO.RCPM.xslice, distance, APAscore, 1:(ncol(APA.DMSO.RCPM.xslice)-1))
          APA.DMSO.RCPM.xlong$distance = (as.numeric(gsub("V", "", APA.DMSO.RCPM.xlong$distance))-11)*APAmatrix.scale
      APA.Ent6.RCPM.xslice= APA.Ent6.RCPM[seq(11, nrow(APA.Ent6.RCPM), 21), ]; APA.Ent6.RCPM.xlong = gather(APA.Ent6.RCPM.xslice, distance, APAscore, 1:(ncol(APA.Ent6.RCPM.xslice)-1))
          APA.Ent6.RCPM.xlong$distance = (as.numeric(gsub("V", "", APA.Ent6.RCPM.xlong$distance))-11)*APAmatrix.scale
      APA.xlong = rbind(APA.DMSO.RCPM.xlong,APA.Ent6.RCPM.xlong)
      
      library(ggplot2)
        ggplot(APA.xlong[grep("CRC_R",APA.xlong$comparison.name),], aes(x=distance, y=APAscore,color=comparison.name))+theme_bw()+stat_smooth(method = lm, formula = y ~ poly(x, 10), se = FALSE)
        ggplot(APA.xlong[grep("CRC_TES",APA.xlong$comparison.name),], aes(x=distance, y=APAscore,color=comparison.name))+theme_bw()+stat_smooth(method = lm, formula = y ~ poly(x, 10), se = FALSE)
      
      
      
    