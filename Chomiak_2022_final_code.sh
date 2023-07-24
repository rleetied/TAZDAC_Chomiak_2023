#Software versions
STAR_2.6.1a_08-27
trim_galore-0.5.0
python-2.7.15
squire-0.9.9.92
bowtie2-2.3.5
samblaster-0.1.24
samtools-1.9
siQ-ChIP-Feb2021
siQ-ChIP-June2022
UCSC_Browser_Tools-2017-03-15
bedtools-v2.25.0
R-4.2.2
biscuit-0.3.16

############ Processing Total RNA-seq reads ################
#Trim PE reads
trim_galore --quality 25 --fastqc --output_dir ../trimmed_fastq --paired $SAMPLE_L000_R1_001.fastq.gz $SAMPLE_L000_R2_001.fastq.gz
#Align totalRNA-seq PE reads with STAR GRCh38.primary_assembly.genome.fa gencode.v29.primary_assembly.annotation.gtf
STAR --genomeDir ../STAR_index_20181015 \                                                                                                                                   
--readFilesIn $SAMPLE_L000_R1_001_val_1.fq.gz $SAMPLE_L000_R2_001_val_2.fq.gz \                                                                                                                          
--twopassMode Basic \                                                                                                                                                                                              
--outReadsUnmapped None \                                                                                                                                                                                          
--runThreadN 8 \                                                                                                                                                                                                   
--readFilesCommand zcat \                                                                                                                                                                                          
--outFileNamePrefix $SAMPLE \                                                                                                                                                                                 
--sjdbOverhang 50 \                                                                                                                                                                                                
--quantMode GeneCounts                                                                                                                                                                                       
--outSAMtype SortedByCoordinate     
#*ReadsPerGene.out.tab files were further processed using edgeR (see Chomiak_Analysis_code.r)

#Quantifying transposable/repeat elements from total RNA-seq
#Alignment
squire Map -1 "$SAMPLE_L000_R1_001_val_1.fq.gz" -2 "$SAMPLE_L000_R2_001_val_2.fq.gz" -f ../alignment_indexes/hg38_index/SQuIRE_index/squire_fetch -o ../aligned/squire/squire_map/bp50/$SAMPLE -b hg38 -r 50 -p 28 -n $SAMPLE
#Quantification                                                                                                                                                                                                                   
squire Count -m ../aligned/squire/squire_map/bp50/$SAMPLE -n $SAMPLE -c ../alignment_indexes/hg38_index/SQuIRE_index/squire_clean -o ../TE_analysis/aligned/squire_count -t ../TE_analysis/aligned/squire_count/tmp -f ../alignment_indexes/hg38_index/SQuIRE_index/squire_fetch -b hg38 -p 28 -r 50          
#Differential expression analysis
squire Call --group1 $Treatment_rep_1,$Treatment_rep_2,$Treatment_rep_3 --group2 $Control_rep_1,$Control_rep_2,$Control_rep_3 --condition1 $Treatment --condition2 $Control --projectname $ControlvsTreatment --pthreads 8 --output_format pdf -i ../TE_analysis/aligned/squire_count --call_folder squire_call/$ControlvsTreatment --verbosity

########## Processing ChIP-sequencing reads ################
#Trim PE reads
trim_galore --quality 25 --fastqc --output_dir ../trimmed_fastq --paired $SAMPLE_L000_R1_001.fastq.gz $SAMPLE_L000_R2_001.fastq.gz
#Align ChIP-seq PE reads with bowtie2
bowtie2 -I 0 -X 700 --end-to-end --sensitive -x ../hg38_fasta_SNAP/bowtie_index/hg38.SNAP -1 $SAMPLE_L000_R1_001_val_1.fq.gz -2 $SAMPLE_L000_R2_001_val_2.fq.gz -S $SAMPLE_PE.bt2.sam
#Mark and remove duplicate reads with samblaster
samblaster -i $SAMPLE_PE.bt2.sam -o $SAMPLE_PE.rmdup.bt2.sam -r
#Sort and index BAM files with samtools
samtools view -S -b $SAMPLE_PE.rmdup.bt2.sam > $SAMPLE_PE.rmdup.bt2.bam   
samtools sort $SAMPLE_PE.rmdup.bt2.bam -o $SAMPLE_PE.rmdup.bt2.sort.bam
samtools index $SAMPLE_PE.rmdup.bt2.sort.bam

########## Preparation of ChIP-seq aligned reads for quantitative siQ-ChIP analysis ###########
#Generate bed file with chr start stop fragment_length from aligned SAM file (for siQ-ChIP-Feb2021 and siQ-ChIP-June2022)
awk -v MAQ=20 '$5>=MAQ && $2==99 || $5>=MAQ && $2==163 {print $3"\t"$4"\t"$4+$9-1}' $SAMPLE_PE.rmdup.bt2.sam | awk '$2<=$3 {print $1"\t"$2"\t"$3"\t"$3-2"}' | sort -k1,1 -k2,2n > $SAMPLE_fragments.bed

#bigwigs used for deeptools analysis and browser shots were generated from the February 2021 version of siQ-ChIP and deposited under GEO:GSE236897
./runCrunch.sh $SAMPLE.IP.bed $SAMPLE.INPUT.bed params.in $SAMPLE
#conversion of bed to bigwig
bedGraphToBigWig SIQ_$SAMPLE.bed hg38.chrom.sizes.sort $SAMPLE.100.bw

#Change in response analysis and determination of enriched regions was conducted using the June 2022 version of siQ-ChIP
#Prepare params.in and EXPlayout files as described on https://github.com/BradleyDickson/siQ-ChIP
#Execute siQ-ChIP on HPC
nohup ./getsiq.sh > Errors.out &
#Determination of conserved regions of change between biological duplicates was conducted with a mixture of bedtools and R scripting
#Isolate repsonse column ($7) from siQ-ChIP output for each biological replicate and write out new peak file for each comparison to control
#Using DAC30 vs Veh comparisons as an example:
awk '{print $2"\t"$3"\t"$4"\t"$7"\t""DAC301vsVeh1"}' DAC301vsVeh1 > DAC301vsVeh1_peaks.bed
awk '{print $2"\t"$3"\t"$4"\t"$7"\t""DAC301vsVeh2"}' DAC301vsVeh2 > DAC301vsVeh2_peaks.bed
#Intersect peaks from treatment comparison to Veh replicates
bedtools intersect -a DAC301vsVeh1_peaks.bed -b DAC301vsVeh2_peaks.bed > DAC301_conserved_peaks_all.bed

#In R, read in the data file and calculate the average log2FC for an individual biological replicate versus the control replicates
library(data.table)
library(dplyr)

DAC301 <- fread("DAC301_conserved_peaks_all.bed")
DAC301 <- mutate(DAC301, avg_FC = (V4 + V9)/2)
DAC301 <- mutate(DAC301, Log2FC = log2(avg_FC))
DAC301 <- DAC301[,c(1:3,12)]
write.table(DAC301, file = "DAC301_Log2FC_all.bed", quote = F, row.names = F, col.names = F, sep = "\t")

#Using bedtools, determine conserved peaks between the treatment biological replicates using 'closest' from bedtools
bedtools closest -a DAC301_log2FC_all.bed -b DAC302_Log2FC_all.bed -d > DAC301.closest.DAC302.Log2FC.bed

#In R, read in the data file and make final cut-offs for determining conserved peaks
DAC30 <- fread("DAC301.closest.DAC302.Log2FC.bed")
DAC30_conserved <- subset(DAC30, V9 <= 200) #peaks between biological replicates need to be within 200 bp of each other
DAC30_conserved <- mutate(DAC30_conserved, avg_log2FC = (V4+V8)/2)
DAC30_conserved <- DAC30_conserved[,c(1:3,10)]

write.table(DAC30_conserved, file = "DAC30_conserved_peaks_Log2FC_all_regions.bed", quote = F, sep = "\t", row.names = F, col.names = F)

#ChromHMM/repliseq enrichment analysis
java -mx4000M -jar ChromHMM.jar OverlapEnrichment -labels HCT116_15_segments_renamed.bed $SAMPLE.directory_with_regions_of_interest $SAMPLE_chromHMM

#R script for generating heatmaps for chromHMM/repliseq overlap enrichment
library(pheatmap)

H3K27me3_change <- fread("H3K27me3_change_enrichment.txt")
colnames(H3K27me3_change) <- c("State", "Genome", "DAC300_decrease", "DAC300_increase","DAC30_decrease","DAC30_increase")
x <- H3K27me3_change[,c(1,5,3,6,4)]
x <- data.frame(x[,-1], row.names = x$State)
x2 <- data.matrix(x)

col <- colorRampPalette(c("white","blue"))(10) #color for decrease
col <- colorRampPalette(c("white","orangered"))(10) #color for increase

pheatmap(x2, col = col, cluster_cols = F, cluster_rows = F, scale = c("none"), breaks = seq(0,10, length.out = 11), 
         fontsize = 12, angle_col = 90, display_numbers = T)

################## Processing Enzymatic-Methyl (EM)-seq reads #######################
#Trim EM-seq PE reads
trim_galore --quality 25 --fastqc --output_dir ../trimmed_fastq --paired $SAMPLE_L000_R1_001.fastq.gz $SAMPLE_L000_R2_001.fastq.gz
#these are directional libs                                                                                                                                                                                        
#Launch biscuit alignment                                     
/biscuit align -@ 30 -R '@RG\tLB:hg38\tID:WGBS_$SAMPLE\tPL:Illumina\tPU:novoseq6000\tSM:$SAMPLE' \                                           
../alignment_indexes/hg38_index/biscuit_index/hg38.final.fa \                                                                                                    
$SAMPLE_L000_R1_001_val_1.fq.gz $SAMPLE_L000_R2_001_val_2.fq.gz | \                                                                                                                                      
samblaster --addMateTags | parallel --tmpdir ${tmpdir} --pipe --tee {} ::: 'samblaster -a -e -u ${SAMPLE[i]}.clipped.fastq -d $SAMPLE.disc.hg38.sam -s $SAMPLE.split.hg38.sam -o /dev/null' \            
'samtools view -hb | samtools sort -@ 8 -m 5G -o ${SAMPLE[i]}.sorted.markdup.withdisc_split_clip.hg38.bam -O BAM -'                                                                                                
                                                                                                                                                                                                                   
#you can collect the split and discordant bams out of smoove too                                                                                                                                                   
#sort and convert to bams                                                                                                                                                                                          
samtools sort -@ 20 -o $SAMPLE.disc.hg38.bam -O BAM $SAMPLE.disc.hg38.sam                                                                                                                                
                                                                                                                                                                                                                   
#index                                                                                                                                                                                                             
samtools index $SAMPLE.disc.hg38.bam                                                                                                                                                                          
                                                                                                                                                                                                                   
#sort and convert to bams                                                                                                                                                                                          
samtools sort -o $SAMPLE.split.hg38.bam -O BAM $SAMPLE.split.hg38.sam                                                                                                                                    
                                                                                                                                                                                                                   
#index                                                                                                                                                                                                             
samtools index $SAMPLE.split.hg38.bam                                                                                                                                                                         
                                                                                                                                                                                                                   
#compress                                                                                                                                                                                                          
pigz -p 40 $SAMPLE.clipped.fastq                                                                                                                                                                              
                                                                                                                                                                                                                   
#index the full bam                                                                                                                                                                                                
samtools index $SAMPLE.sorted.markdup.withdisc_split_clip.hg38.bam                                                                                                                                            
                                                                                                                                                                                                                   
#clean up the sam files                                                                                                                                                                                              
if [[ -f $SAMPLE.disc.hg38.bam ]] && [[ -f $SAMPLE.split.hg38.bam ]]; then rm $SAMPLE*.sam; fi                                                                                                      

#combine sample output and call CpG methylation
biscuit pileup -@ 30 -o 20201002_ChIP_BS_CpG.vcf ../alignment_indexes/hg38_index/biscuit_index/hg38.final.fa Veh1.sorted.markdup.withdisc_split_clip.hg38.bam Veh2.sorted.markdup.withdisc_split_clip.hg38.bam DAC301.sorted.markdup.withdisc_split_clip.hg38.bam DAC302.sorted.markdup.withdisc_split_clip.hg38.bam DAC3001.sorted.markdup.withdisc_split_clip.hg38.bam DAC3002.sorted.markdup.withdisc_split_clip.hg38.bam

bgzip 20201002_ChIP_BS_CpG.vcf

tabix -p vcf 20201002_ChIP_BS_CpG.vcf.gz

biscuit vcf2bed 20201002_ChIP_BS_CpG.vcf.gz > 20201002_ChIP_BS_CpG_all_samples.bed

#merge CpG calls from both strands
biscuit mergecg ../alignment_indexes/hg38_index/biscuit_index/hg38.final.fa 20201002_ChIP_BS_CpG_all_samples.bed > 20220131_mergeCG_output.bed



