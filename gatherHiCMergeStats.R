###########################################
#  Collect mergestat info for AQuA-HiChIP #
#  Berkley Gryder (gryderart#gmail.com)   #
#  June 2019                              #
###########################################

## 1. Set project global parameters
      root.dir = "K:/projects/HiC/"
      setwd(root.dir)
      output.dir = "projects/stats/"

## 2. Gather projects and sample lists
      HiC.allfolders = data.frame(list.dirs(path="projects", full.names=T, recursive=T))
      colnames(HiC.allfolders) = "folders"
      library(stringr)
      HiC.allfolders$project = str_split_fixed(HiC.allfolders$folders, "/", 3)[,2]
      
      HiCsamples = HiC.allfolders[grep("HiCpro_OUTPUT/hic_results/stats/", HiC.allfolders$folders), ]
      HiCsamples = HiCsamples[grep("/tmp", HiCsamples$folders,invert = T), ]
      HiCsamples$sample = str_split_fixed(HiCsamples$folders, "/", 6)[,6]
      
      HiCsamples$mergestat.path = paste(root.dir,"projects/",HiCsamples$project,"/HiCpro_OUTPUT/hic_results/stats/",HiCsamples$sample,"/",HiCsamples$sample,"_allValidPairs.mergestat",sep="")
      HiCsamples$mergestat.exists = file.exists(HiCsamples$mergestat.path); HiCsamples = subset(HiCsamples, HiCsamples$mergestat.exists == TRUE)
      
      HiCsamples.mm10 = HiCsamples[grep("_mm10", HiCsamples$project,invert = F), ]; HiCsamples.mm10$genome = "mm10"
      HiCsamples.hg19 = HiCsamples[grep("_mm10", HiCsamples$project,invert = T), ]; HiCsamples.hg19$genome = "hg19"
      HiCsamples = rbind(HiCsamples.hg19, HiCsamples.mm10)
      
      HiCsamples.AQuApaired = subset(HiCsamples.hg19, HiCsamples.hg19$sample %in% HiCsamples.mm10$sample)
        AQuA.mm10.reftable = HiCsamples.mm10[,c("sample","mergestat.path")]; colnames(AQuA.mm10.reftable)=c("sample","mergestat.path.mm10")
        library(plyr)
        HiCsamples.AQuApaired =  join(HiCsamples.AQuApaired, AQuA.mm10.reftable, by="sample")
        
      sample.list = HiCsamples.AQuApaired$sample

## 3. Obtain spike in read counts
      
      mergestat.all = as.data.frame(read.table(HiCsamples[1,c("mergestat.path")], sep="\t", header=F))
      mergestat.all = as.data.frame(mergestat.all$V1); colnames(mergestat.all) = "statistic"
      
      lapply(sample.list, function(x) {
        ##load and merge sample data
        sample.df = subset(HiCsamples.AQuApaired, HiCsamples.AQuApaired$sample == x)
        mergestat <- read.table(sample.df$mergestat.path, sep="\t", header=F)
        mergestat.sample = as.data.frame(mergestat[,2])
        removable.string = "Sample_" ; sample.name = gsub(removable.string,"",x) ; colnames(mergestat.sample) = c(sample.name)
        
        mergestat.all <<- cbind(mergestat.all,mergestat.sample)
      })
      
      mergestat.all.mm10 = as.data.frame(read.table(HiCsamples[1,c("mergestat.path")], sep="\t", header=F))
      mergestat.all.mm10 = as.data.frame(mergestat.all.mm10$V1); colnames(mergestat.all.mm10) = "statistic"
      mergestat.all.mm10$statistic = as.factor(paste(mergestat.all.mm10$statistic,"mm10",sep="."))
      
      lapply(sample.list, function(x) {
        ##load and merge sample data
        sample.df = subset(HiCsamples.AQuApaired, HiCsamples.AQuApaired$sample == x)
        mergestat <- read.table(sample.df$mergestat.path.mm10, sep="\t", header=F)
        mergestat.sample = as.data.frame(mergestat[,2])
        removable.string = "Sample_" ; sample.name = gsub(removable.string,"",x) ; colnames(mergestat.sample) = c(sample.name)
        
        mergestat.all.mm10 <<- cbind(mergestat.all.mm10,mergestat.sample)
      })
      
      mergestat.AQuA <- rbind(mergestat.all,mergestat.all.mm10)
      
      write.table(mergestat.AQuA, file = "projects/stats/All.samples.mergestat.AQuA.txt",sep="\t", col.names = T, row.names = F)
      
      #transform the table so AQuA stats can be calculated for each sample!
      library(reshape2)
      
      mergestat.AQuA.t = as.data.frame(t(mergestat.AQuA[,2:ncol(mergestat.AQuA)]))
        colnames(mergestat.AQuA.t) = mergestat.AQuA$statistic
        mergestat.AQuA.t$prct_hg19 = mergestat.AQuA.t$valid_interaction_rmdup/(mergestat.AQuA.t$valid_interaction_rmdup+mergestat.AQuA.t$valid_interaction_rmdup.mm10)
        mergestat.AQuA.t$prct_mm10 = mergestat.AQuA.t$valid_interaction_rmdup.mm10/(mergestat.AQuA.t$valid_interaction_rmdup+mergestat.AQuA.t$valid_interaction_rmdup.mm10)
          
        mergestat.AQuA.t$prct_dup_hg19 = 1-mergestat.AQuA.t$valid_interaction_rmdup/mergestat.AQuA.t$valid_interaction
        mergestat.AQuA.t$prct_dup_mm10 = 1-mergestat.AQuA.t$valid_interaction_rmdup.mm10/mergestat.AQuA.t$valid_interaction.mm10
        
        mergestat.AQuA.t$AQuA.factor = mergestat.AQuA.t$valid_interaction_rmdup/mergestat.AQuA.t$valid_interaction_rmdup.mm10
      
      write.table(mergestat.AQuA.t, file = "projects/stats/All.samples.mergestat.AQuA.t.txt",sep="\t", col.names = T, row.names = T)
      
      
      
      
      
      
      
