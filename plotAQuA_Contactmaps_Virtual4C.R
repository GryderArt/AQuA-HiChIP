###HiC analysis, setup for compatibility with HiC-pro and Juicebox/juicer

### 1.Set up project and samples. Configure variables to adopt code into your own file system.

      setwd("K:/projects/HiC/projects/")
      project.title = "RH4_H3K27ac_HiChIP"
      project.folder = paste(project.title,"/HiCpro_OUTPUT/hic_results/data/",sep="")
      sample.list = list.dirs(path = project.folder, full.names = F, recursive = F)
      project.data.folder = paste(project.title,"/HiCpro_OUTPUT/hic_results/data/",sep="")

      project.title.mm10 = paste(project.title,"_mm10",sep="")
      project.folder.mm10 = paste(project.title.mm10,"/HiCpro_OUTPUT/hic_results/data/",sep="")  #folder changed to "stats" in newer HiC Pro pipeline for mergestat 
      sample.list.mm10 = list.dirs(path = project.folder.mm10, full.names = F, recursive = F)
      project.data.folder.mm10 = paste(project.title.mm10,"/HiCpro_OUTPUT/hic_results/data/",sep="")
      
### 2. Obtain spike in read counts
  
      mergestat.all = as.data.frame(read.table(paste(project.folder,"/",sample.list[1],"/",sample.list[1],"_allValidPairs.mergestat",sep=""), sep="\t", header=F)); mergestat.all = as.data.frame(mergestat.all$V1)
      
      lapply(sample.list, function(x) {
        ##load and merge sample data
        mergestat <- read.table(paste(project.folder,"/",x,"/",x,"_allValidPairs.mergestat",sep=""), sep="\t", header=F)
        mergestat.sample = as.data.frame(mergestat[,2])
        removable.string = "Sample_" ; sample.name = gsub(removable.string,"",x) ; colnames(mergestat.sample) = c(sample.name)
        mergestat.all <<- cbind(mergestat.all,mergestat.sample)
        })
      
      lapply(sample.list.mm10, function(x) {
        ##load and merge sample data
        mergestat <- read.table(paste(project.folder.mm10,"/",x,"/",x,"_allValidPairs.mergestat",sep=""), sep="\t", header=F)
        mergestat.sample = as.data.frame(mergestat[,2])
        removable.string = "Sample_" ; sample.name = gsub(removable.string,"mm10_",x) ; colnames(mergestat.sample) = c(sample.name)
        mergestat.all <<- cbind(mergestat.all,mergestat.sample)
        })
      
      lapply(sample.list, function(x) {
        ##get percentage hg19 and mm10
        removable.string = "Sample_" ; hg19 = gsub(removable.string,"",x) ; mm10 = gsub(removable.string,"mm10_",x)
        mergestat.prct = as.data.frame(mergestat.all[,c(hg19,mm10)])
        mergestat.prct$prct_hg19 = mergestat.prct[,1]/(mergestat.prct[,1]+mergestat.prct[,2])
        mergestat.prct$prct_mm10 = mergestat.prct[,2]/(mergestat.prct[,1]+mergestat.prct[,2])
        colnames(mergestat.prct) = paste(colnames(mergestat.prct),x,sep="_")
        mergestat.all <<- cbind(mergestat.all,mergestat.prct[,3:4])
      })
  
      write.table(mergestat.all, paste("stats/",project.title,".mergestat.txt", sep=""),sep="\t", col.names = T, row.names = F)
      
### 3. Load and transform sparse matrix (from juicer: java -jar /usr/local/apps/juicer/juicer-1.5.6/scripts/juicer_tools.jar dump)
      #if HiCcompare isn't installed, remove comments from the next 2 lines of code and run them.
      #source("https://bioconductor.org/biocLite.R")
      #biocLite("HiCcompare")
      library(HiCcompare)
      sample1 = "Sample_RH4_D6_H3K27ac_HiChIP_HKJ22BGX7";sample2 = "Sample_RH4_Ent6_H3K27ac_HiChIP_HKJ22BGX7"
          removable.string = "Sample_" ; name1 = gsub(removable.string,"",sample1) ; name2 = gsub(removable.string,"",sample2) 
          mm10_name1 = gsub(removable.string,"mm10_",sample1) ; mm10_name2 = gsub(removable.string,"mm10_",sample2)
          matrix.suffix = ".chr11.176MB-178MB.5KB.matrix.txt"
      sparMat1 = read.table(paste(project.data.folder,sample1,"/",sample1,matrix.suffix,sep=""), header=F, sep="\t"); sparMat2 = read.table(paste(project.data.folder,sample2,"/",sample2,matrix.suffix,sep=""), header=F, sep="\t")
      denMat1 <- sparse2full(sparMat1, hic.table = FALSE, column.name = NA); denMat2 <- sparse2full(sparMat2, hic.table = FALSE, column.name = NA)
      denMat1.CPM = denMat1*1000000/(mergestat.all[2,name1]+mergestat.all[2,mm10_name1]) ; denMat2.CPM = denMat2*1000000/(mergestat.all[2,name2]+mergestat.all[2,mm10_name2])
      AQuAfactor1 = mergestat.all[2,name1]/mergestat.all[2,mm10_name1] ; AQuAfactor2 = mergestat.all[2,name2]/mergestat.all[2,mm10_name2]
      denAQuA1 = denMat1.CPM*AQuAfactor1; denAQuA2 = denMat2.CPM*AQuAfactor2

### 4. Build matrices with AQuA maximums, plot HEATMAPS
            library(pheatmap)
            library(ggplot2)
            library(gridExtra)
            quant_cut = 0.98 #caps the contact map plot values at a given percentile
  
            AQuAmax = max(quantile(denAQuA1, probs = c(quant_cut)),quantile(denAQuA2, probs = c(quant_cut)))
            denAQuA1max = denAQuA1; denAQuA1max[denAQuA1max>AQuAmax] <- AQuAmax
            denAQuA2max = denAQuA2; denAQuA2max[denAQuA2max>AQuAmax] <- AQuAmax

            denDeltaAQuA = denAQuA2 - denAQuA1; 
               DeltaAQuA.max = AQuAmax*0.5
            denDeltaAQuA.max = denDeltaAQuA; denDeltaAQuA.max[denDeltaAQuA.max>DeltaAQuA.max] <- DeltaAQuA.max
                                             denDeltaAQuA.max[denDeltaAQuA.max<(-DeltaAQuA.max)] <- (-DeltaAQuA.max)
            
            breakList = seq(-(DeltaAQuA.max*1.05), (DeltaAQuA.max*1.05), by = AQuAmax/100)

            #combine the plots!
            p1 = pheatmap((denAQuA1max), cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white", "red"))(50))
            p2 = pheatmap((denAQuA2max), cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white", "red"))(50))
            p3 = pheatmap(denDeltaAQuA.max , cluster_rows = F, cluster_cols = F, breaks = breakList, color = colorRampPalette(c("dodgerblue", "white", "mediumvioletred"))(length(breakList)))
            grid.arrange(arrangeGrob(p1[[4]],p2[[4]],p3[[4]],nrow=1,top = paste("AQuA: ",name1,", ",name2,", and delta from *",matrix.suffix,sep="")))
            
            #save it!!
            g <- arrangeGrob(p1[[4]],p2[[4]],p3[[4]],nrow=1,top = paste("AQuA: ",name1,", ",name2,", and delta from *",matrix.suffix,sep=""))
            ggsave(file=paste(project.data.folder,sample2,"/",sample2,matrix.suffix,".3AQuAplots.pdf",sep=""), g, width = 14, height = 4.5) #saves g
            
      #repeat with no AQuA
           
            CPMmax = max(quantile(denMat1.CPM, probs = c(quant_cut)),quantile(denMat2.CPM, probs = c(quant_cut)))
            denMat1.CPMmax = denMat1.CPM; denMat1.CPMmax[denMat1.CPMmax>CPMmax] <- CPMmax
            denMat2.CPMmax = denMat2.CPM; denMat2.CPMmax[denMat2.CPMmax>CPMmax] <- CPMmax
            
            denDelta.CPM = denMat2.CPM - denMat1.CPM
               Delta.CPM.max = CPMmax*0.8
               denDelta.CPM.max = denDelta.CPM; denDelta.CPM.max[denDelta.CPM.max>Delta.CPM.max] <- Delta.CPM.max
                                                denDelta.CPM.max[denDelta.CPM.max<(-Delta.CPM.max)] <- (-Delta.CPM.max)
               
            breakList.CPM = seq(-(Delta.CPM.max*1.05), (Delta.CPM.max*1.05), by = CPMmax/100)

            #combine the plots!
            p1 = pheatmap((denMat1.CPMmax), cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white", "red"))(50))
            p2 = pheatmap((denMat2.CPMmax), cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white", "red"))(50))
            p3 = pheatmap(denDelta.CPM.max, cluster_rows = F, cluster_cols = F, breaks = breakList.CPM, color = colorRampPalette(c("dodgerblue", "white", "mediumvioletred"))(length(breakList.CPM)))
            grid.arrange(arrangeGrob(p1[[4]],p2[[4]],p3[[4]],nrow=1,top = paste("HiChIP: ",name1,", ",name2,", and delta from *",matrix.suffix,sep="")))
            
            #save it!!
            g <- arrangeGrob(p1[[4]],p2[[4]],p3[[4]],nrow=1,top = paste("HiChIP: ",name1,", ",name2,", and delta from *",matrix.suffix,sep=""))
            ggsave(file=paste(project.data.folder,sample2,"/",sample2,matrix.suffix,".3plots.CPM.pdf",sep=""), g, width = 14, height = 4.5) #saves g
            
      
### 5. Extract a virtual 4C viewpoint, make bedgraphs, plot VIRTUAL4C
      virt4C_viewpoint_chr = 'chr16'
      virt4C_viewpoint = '960000'
      df_AQuA1 = as.data.frame(denAQuA1, row.names = row.names(denAQuA1), col.names = col.names(denAQuA1))
      df_AQuA2 = as.data.frame(denAQuA2, row.names = row.names(denAQuA2), col.names = col.names(denAQuA2))
      
      virt4C.bedgraph1 = as.data.frame(df_AQuA1[,virt4C_viewpoint])
        colnames(virt4C.bedgraph1) = "AQuA_contact_freq"
        virt4C.bedgraph1$chr = virt4C_viewpoint_chr
        virt4C.bedgraph1$start = as.numeric(rownames(df_AQuA1))
        virt4C_binsize = virt4C.bedgraph1[2,c('start')] - virt4C.bedgraph1[1,c('start')]
        virt4C.bedgraph1$stop = virt4C.bedgraph1$start + virt4C_binsize
        virt4C.bedgraph1 = virt4C.bedgraph1[,c(2,3,4,1)]
      
      virt4C.bedgraph2 = as.data.frame(df_AQuA2[,virt4C_viewpoint])
        colnames(virt4C.bedgraph2) = "AQuA_contact_freq"
        virt4C.bedgraph2$chr = virt4C_viewpoint_chr
        virt4C.bedgraph2$start = as.numeric(rownames(df_AQuA2))
        virt4C_binsize = virt4C.bedgraph2[2,c('start')] - virt4C.bedgraph2[1,c('start')]
        virt4C.bedgraph2$stop = virt4C.bedgraph2$start + virt4C_binsize
        virt4C.bedgraph2 = virt4C.bedgraph2[,c(2,3,4,1)]
        
      window_coordinates = paste(virt4C_viewpoint_chr,":",virt4C.bedgraph1[1,c('start')],"-",max(virt4C.bedgraph1$stop),sep="")
    ## spline to smooth
          spline1 <- as.data.frame(spline(virt4C.bedgraph1$start, virt4C.bedgraph1$AQuA_contact_freq))
          spline2 <- as.data.frame(spline(virt4C.bedgraph2$start, virt4C.bedgraph2$AQuA_contact_freq))
          spline_delta = spline1; spline_delta$y = (spline2$y-spline1$y)
          spline_fold = spline1; spline_fold$y = (spline2$y+0.1)/(spline1$y+0.1) ; spline_fold = spline_fold[4:(nrow(spline_fold)-3),]
              
          spline1.bedgraph = spline1; spline1.bedgraph[,1]=round(x = spline1.bedgraph[,1], digits = 0)+virt4C_binsize/2; spline1.bedgraph$chr = virt4C_viewpoint_chr ; spline1.bedgraph = spline1.bedgraph[,c(3,1,1,2)]
          spline2.bedgraph = spline2; spline2.bedgraph[,1]=round(x = spline2.bedgraph[,1], digits = 0)+virt4C_binsize/2; spline2.bedgraph$chr = virt4C_viewpoint_chr ; spline2.bedgraph = spline2.bedgraph[,c(3,1,1,2)]
          spline_delta.bedgraph = spline_delta; spline_delta.bedgraph[,1]=round(x = spline_delta.bedgraph[,1], digits = 0)+virt4C_binsize/2; spline_delta.bedgraph$chr = virt4C_viewpoint_chr ; spline_delta.bedgraph = spline_delta.bedgraph[,c(3,1,1,2)]
          spline_fold.bedgraph = spline_fold; spline_fold.bedgraph[,1]=round(x = spline_fold.bedgraph[,1], digits = 0)+virt4C_binsize/2; spline_fold.bedgraph$chr = virt4C_viewpoint_chr ; spline_fold.bedgraph = spline_fold.bedgraph[,c(3,1,1,2)]
    ## write out bedgraphs  
      write.table(virt4C.bedgraph1, file=paste(project.data.folder,sample1,virt4C_viewpoint_chr,virt4C_viewpoint,".bedgraph",sep=""), sep="\t", row.names=FALSE, col.names=F, quote=FALSE)
      write.table(virt4C.bedgraph2, file=paste(project.data.folder,sample2,virt4C_viewpoint_chr,virt4C_viewpoint,".bedgraph",sep=""), sep="\t", row.names=FALSE, col.names=F, quote=FALSE)
   
      write.table(spline1.bedgraph, file=paste(project.data.folder,sample1,virt4C_viewpoint_chr,virt4C_viewpoint,".spline.bedgraph",sep=""), sep="\t", row.names=FALSE, col.names=F, quote=FALSE)
      write.table(spline2.bedgraph, file=paste(project.data.folder,sample2,virt4C_viewpoint_chr,virt4C_viewpoint,".spline.bedgraph",sep=""), sep="\t", row.names=FALSE, col.names=F, quote=FALSE)
      write.table(spline_delta.bedgraph, file=paste(project.data.folder,sample1,"_",sample2,virt4C_viewpoint_chr,virt4C_viewpoint,".spline_delta.bedgraph",sep=""), sep="\t", row.names=FALSE, col.names=F, quote=FALSE)
      write.table(spline_fold.bedgraph, file=paste(project.data.folder,sample1,"_",sample2,virt4C_viewpoint_chr,virt4C_viewpoint,".spline_fold.bedgraph",sep=""), sep="\t", row.names=FALSE, col.names=F, quote=FALSE)
      
    ## make some plots
      
      library(ggplot2)
      ggplot(virt4C.bedgraph1, aes(x=(start+virt4C_binsize/2),y=AQuA_contact_freq))+geom_line()+geom_line(data=virt4C.bedgraph2, color = 'red')+
        theme_bw()+geom_vline(xintercept = as.numeric(virt4C_viewpoint),color = "blue")+geom_vline(xintercept = (as.numeric(virt4C_viewpoint)+virt4C_binsize),color = "blue")
      
      ggplot(spline1, aes(x=(x+virt4C_binsize/2),y=y))+geom_line()+geom_line(data=spline2, color = 'red')+
        theme_bw()+geom_vline(xintercept = as.numeric(virt4C_viewpoint),color = "blue")+geom_vline(xintercept = (as.numeric(virt4C_viewpoint)+virt4C_binsize),color = "blue")
      
      ggplot(spline1, aes(x=(x+virt4C_binsize/2),y=y))+geom_line(data=spline1, color = 'red')+geom_line(data=spline_delta, color = 'purple')+
        theme_bw()+geom_vline(xintercept = as.numeric(virt4C_viewpoint),color = "blue")+geom_vline(xintercept = (as.numeric(virt4C_viewpoint)+virt4C_binsize),color = "blue")
      
      #check spline coordinates by overlapping with raw bedgraph
      ggplot(virt4C.bedgraph1, aes(x=(start+virt4C_binsize/2),y=AQuA_contact_freq))+geom_line()+geom_line(data=spline1, aes(x=x,y=y), color = 'red')+
        theme_bw()+geom_vline(xintercept = as.numeric(virt4C_viewpoint),color = "blue")+geom_vline(xintercept = (as.numeric(virt4C_viewpoint)+virt4C_binsize),color = "blue")
      
      
      
      
      
      
    
    