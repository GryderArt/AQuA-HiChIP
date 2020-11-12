#########################################
#  batch convert HiCpro to Juicebox     #
#  (.validpairs to .hic)                #
#  Berkley Gryder (gryderart#gmail.com) #
#  June 2019                            #
#########################################

## 1. Set project global parameters
      setwd("K:/projects/HiC/")
      output.dir = "manage_samples/"

## 2. Gather validpair and .hic file lists
      HiC.allfolders = data.frame(list.dirs(path="projects", full.names=T, recursive=T))
      colnames(HiC.allfolders) = "folders"
      library(stringr)
      HiC.allfolders$project = str_split_fixed(HiC.allfolders$folders, "/", 3)[,2]
      
      HiCsamples = HiC.allfolders[grep("HiCpro_OUTPUT/hic_results/data/", HiC.allfolders$folders), ]
      HiCsamples = HiCsamples[grep("/tmp", HiCsamples$folders,invert = T), ]
      HiCsamples$sample = str_split_fixed(HiCsamples$folders, "/", 6)[,6]
      
      HiCsamples$validpairs.path = paste("K:/projects/HiC/projects/",HiCsamples$project,"/HiCpro_OUTPUT/hic_results/data/",HiCsamples$sample,"/",HiCsamples$sample,".allValidPairs",sep="")
      HiCsamples$validpair.exists = file.exists(HiCsamples$validpairs.path)
      
      HiCsamples$hic.path = paste("K:/projects/HiC/projects/",HiCsamples$project,"/HiCpro_OUTPUT/hic_results/data/",HiCsamples$sample,"/",HiCsamples$sample,".allValidPairs.hic",sep="")
      HiCsamples$hic.exists = file.exists(HiCsamples$hic.path)
      
## 3. Define unconverted validpair files and genome build
      
      HiCsamples.convert = subset(HiCsamples, HiCsamples$validpair.exists == TRUE); HiCsamples.convert = subset(HiCsamples.convert, HiCsamples.convert$hic.exists == FALSE)

      HiCsamples.convert.mm10 = HiCsamples.convert[grep("_mm10", HiCsamples.convert$project,invert = F), ]; HiCsamples.convert.mm10$genome = "mm10"
      HiCsamples.convert.hg19 = HiCsamples.convert[grep("_mm10", HiCsamples.convert$project,invert = T), ]; HiCsamples.convert.hg19$genome = "hg19"
      HiCsamples.convert = rbind(HiCsamples.convert.hg19, HiCsamples.convert.mm10)
      
## 4. Build sbatch commands for convertHiCPro2Juicebox.sh

      sbatch.config = "sbatch -J HiC2Juice --partition=ccr --time=96:00:00 --mem=121g --cpus-per-task=4  --gres=lscratch:200 "
      script = " /data/khanlab/projects/HiC/scripts/convertHiCPro2Juicebox.sh"
      
      HiCsamples.convert$sbatch.command = paste(sbatch.config,
                                                "--export=HIC_PROJECT=",HiCsamples.convert$project,
                                                ",GENOME=",HiCsamples.convert$genome,
                                                ",HIC_SAMPLE=",HiCsamples.convert$sample,
                                                script,sep = "")
      
      #write to file
      write.table(HiCsamples.convert$sbatch.command,file = paste(output.dir,"multisbatch_HiCPro2Juicebox.sh",sep = ""),col.names = F, row.names = F, quote = F)
      
## 5. Instructions for command line execution of the shell script we just made:
      #log in to biowulf cluster, then navigate to the output folder where the shell script (.sh file) is:
      #>>>command line>>>$ cd /data/khanlab/projects/HiC/manage_samples/
      #then run it with:
      #>>>command line>>>$ sh multisbatch_HiCPro2Juicebox.sh
      
      