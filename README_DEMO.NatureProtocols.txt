#############################################
#   README for AQuA-HiChIP code, 2018		#
#	by Berkley Gryder, gryderart@gmail.com	#
#############################################

### System requirements:
	a. Software dependencies to install: 
		HiC-Pro 2.10.0
		bowtie  2-2.3.4
		R 3.5.0
		openmpi 3.0.0 
		GSL 2.4
		gcc  7.2.0
		samtools 1.8
		juicer  1.5.6
		python 2.7
	b. Harware suggestion:
		HiC-Pro works well on a cluster node with 4 CPUs, 32 Gb of memory, and 200 Gb scratch disk space
		***not tested for operation on a "normal" desktop computer***
	c. Installation:
		Follow instructions for installation at the following websites:
			https://github.com/nservant/HiC-Pro
			https://github.com/aidenlab/juicer/wiki/Juicer-Tools-Quick-Start
			https://www.rstudio.com/products/rstudio/download/#download
			http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2

### DEMO with instructions for use

	***CODE IS SENSITIVE TO FOLDER PATHS, WHICH MUST BE UPDATED FROM THE DEMO TO MATCH YOUR LOCAL MACHINE***

	d. Download FASTQ files from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120770
	e. Place FASTQ files into directories (2 fastq files per directory) with the following structure:
		projects/
		├── RH4_H3K27ac_HiChIP_hg19
		│   └── DATA
		│       ├── Sample_RH4_D6_H3K27ac_HiChIP_HKJ22BGX7
		│       └── Sample_RH4_Ent6_H3K27ac_HiChIP_HKJ22BGX7
		└── RH4_H3K27ac_HiChIP_mm10
			└── DATA
				├── Sample_RH4_D6_H3K27ac_HiChIP_HKJ22BGX7
				└── Sample_RH4_Ent6_H3K27ac_HiChIP_HKJ22BGX7
			
				***note: this uses the same FASTQ files for hg19 and mm10 in parallel***
	f. Prepare configuration file for human and mouse genomes.  
			DEMO configuration files and chrom*.sizes files	here:
				https://github.com/GryderArt/AQuA-HiChIP/tree/master/reference_files
	g. Prepare bowtie 2 indexes for hg19 and mm19 (http://bowtie-bio.sourceforge.net/tutorial.shtml#newi)
			(or, locate local instances of these already in use, and adjust configuration file accordingly)
	h. Prepare in silico digested genomes
		example code for digesting hg19: /usr/local/Anaconda/envs_app/hicpro/2.10.0/HiC-Pro_2.10.0/bin/utils/digest_genome.py -r dpnii  -o dpnii.ucsc.hg19.bed /data/khanlab/projects/HiC/reference_files/bowtie2_index/ucsc.hg19.fasta
	i. Run HiC-Pro
		-example code for mm10: /usr/local/Anaconda/envs_app/hicpro/2.10.0/HiC-Pro_2.10.0/bin/HiC-Pro -i /data/khanlab/projects/HiC/projects/RH4_H3K27ac_HiChIP_mm10/DATA/ -o /data/khanlab/projects/HiC/projects/RH4_H3K27ac_HiChIP_mm10/HiCpro_OUTPUT/ -c /data/khanlab/projects/HiC/reference_files/config_khanlab_mm10.txt -p
		-this will create a shell script, so move folders: cd /data/khanlab/projects/HiC/projects/RH4_H3K27ac_HiChIP_mm10/HiCpro_OUTPUT/
		-then, run the shell script with: sbatch -J HiCstep1 --time=24:00:00 --mem=121g --cpus-per-task=4  --gres=lscratch:200 HiCPro_step1_mm10_HiCpro.sh
		-once step 1 finishes, then run: sbatch -J HiCstep2 --time=24:00:00 --mem=121g --cpus-per-task=4  --gres=lscratch:200 HiCPro_step2_mm10_HiCpro.sh
	j. Convert to juicer compatible .hic format
		example code for hg19: /usr/local/Anaconda/envs_app/hicpro/2.10.0/HiC-Pro_2.10.0/bin/utils/hicpro2juicebox.sh -i /data/khanlab/projects/HiC/projects/RH4_H3K27ac_HiChIP_hg19/HiCpro_OUTPUT/hic_results/data/Sample_RH4_D6_H3K27ac_HiChIP_HKJ22BGX7/Sample_RH4_D6_H3K27ac_HiChIP_HKJ22BGX7_allValidPairs -g hg19 -j /usr/local/apps/juicer/juicer-1.5.6/scripts/juicer_tools.jar
	k. Extraction of DEMO matrix using the juicer tool "dump"
		example code: java -jar /usr/local/apps/juicer/juicer-1.5.6/scripts/juicer_tools.jar dump observed NONE /data/khanlab/projects/HiC/projects/RH4_H3K27ac_HiChIP/HiCpro_OUTPUT/hic_results/data/Sample_RH4_Ent6_H3K27ac_HiChIP_HKJ22BGX7/Sample_RH4_Ent6_H3K27ac_HiChIP_HKJ22BGX7_allValidPairs.hic 11:17600000:17800000 11:17600000:17800000 BP 5000 Sample_RH4_Ent6_H3K27ac_HiChIP_HKJ22BGX7.chr11.176MB-178MB.5KB.matrix.txt
	l. Open R-studio to visualize heatmap
		-download here: https://github.com/GryderArt/AQuA-HiChIP/blob/master/plotAQuA_Contactmaps_Virtual4C_DEMO.R
		-to run a demo starting with pre-made matrices, create folders to mimic HiC-pro structure:
				./projects/HiC/projects/RH4_H3K27ac_HiChIP/HiCpro_OUTPUT/hic_results/data/
				./projects/HiC/projects/RH4_H3K27ac_HiChIP_mm10/HiCpro_OUTPUT/hic_results/data/
		-into these folders, place sample DEMO matrix for control and treated cells here: https://github.com/GryderArt/AQuA-HiChIP/tree/master/DEMO_data
		-to get mm10 and hg19 sample stats in DEMO without running HiC-pro, download file: https://github.com/GryderArt/AQuA-HiChIP/blob/master/DEMO_data/mergestat.HiChIP.all.txt
			place mergestat.HiChIP.all.txt in this folder:  ./projects/HiC/projects/
		-run R-code: Step 1, 3 and 4 to generate heatmaps with AQuA normalization
		-run R-code: Step 5 to demonstrate the Virtual 4C code
		-plots should resemble those inline in the Protocol text
	



			
			

	

			
		
		