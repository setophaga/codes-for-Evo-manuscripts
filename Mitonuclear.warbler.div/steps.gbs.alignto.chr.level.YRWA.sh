#step 1: demultiplexing
perl ../tools/GBS_demultiplexer_30base.pl plate1_barcode.txt HI.4196.002.towa_hewa_hz_plate1_R1.fastq HI.4196.002.towa_hewa_hz_plate1_R2.fastq clean_data/silu_plate1
perl ../tools/GBS_demultiplexer_30base.pl plate2_barcode.txt HI.4411.005.towa_hewa_hz_plate2_R1.fastq HI.4411.005.towa_hewa_hz_plate2_R2.fastq clean_data/silu_plate2
perl ../tools/GBS_demultiplexer_30base.pl plate3_barcode.txt HI.4411.006.towa_hewa_hz_plate3_R1.fastq HI.4411.006.towa_hewa_hz_plate3_R2.fastq clean_data/silu_plate3
perl ../tools/GBS_demultiplexer_30base.pl plate4_barcode.txt HI.4524.003.towa_hewa_hz_plate4_R1.fastq HI.4524.003.towa_hewa_hz_plate4_R2.fastq clean_data/silu_plate4
perl GBS_demultiplexer_30base.pl maddie.p1.barcode.txt HI.4840.005.TOWA_I_R1.fastq HI.4840.005.TOWA_I_R2.fastq clean_data/maddie_plate1
perl GBS_demultiplexer_30base.pl maddie.p2.barcode.txt HI.4840.007.TOWA_II_R1.fastq HI.4840.007.TOWA_II_R2.fastq clean_data/maddie_plate2


#step 2: trim, 2.1: make prefix list
awk '{print "silu_plate1_"$1}' plate1_barcode.txt > prefix.silu.p1.list
awk '{print "silu_plate2_"$1}' plate2_barcode.txt > prefix.silu.p2.list
awk '{print "silu_plate3_"$1}' plate3_barcode.txt > prefix.silu.p3.list
awk '{print "silu_plate4_"$1}' plate4_barcode.txt > prefix.silu.p4.list
awk '{print "maddie_plate1_"$1}' barcode.maddie.plate1 > prefix.maddie.p1.list
awk '{print "maddie_plate2_"$1}' barcode.maddie.plate2 > prefix.maddie.p2.list

#step2.2: trim
while read prefix
do
java -jar tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads 1 clean_data/"$prefix"_R1.fastq clean_data/"$prefix"_R2.fastq clean_data_trim/"$prefix"_R1.fastq clean_data_trim/"$prefix"_R1_unpaired.fastq clean_data_trim/"$prefix"_R2.fastq clean_data_trim/"$prefix"_R2_unpaired.fastq TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:30
done < prefix.silu.p1.list

while read prefix
do
java -jar tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads 1 clean_data/"$prefix"_R1.fastq clean_data/"$prefix"_R2.fastq clean_data_trim/"$prefix"_R1.fastq clean_data_trim/"$prefix"_R1_unpaired.fastq clean_data_trim/"$prefix"_R2.fastq clean_data_trim/"$prefix"_R2_unpaired.fastq TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:30
done < prefix.maddie.p1.list
#step 3 algin

lane="silu3"
runbarcode="C1JDVACXX"

while read prefix
do
## run bwa
#ref=/data/chend8/GBS/gbs.2020.ref.dec.2020/ref/mywagenomev2.1.fa
bwa mem -M ref/mywagenomev2.1.fa  clean_data_trim/"$prefix"_R1.fastq clean_data_trim/"$prefix"_R2.fastq > sam/"$prefix".sam

## add read group headers, convert to bam, sort and index
picard AddOrReplaceReadGroups I=sam/"$prefix".sam O=bam/"$prefix".bam RGID=$lane RGPL=ILLUMINA RGLB=LIB."$prefix" RGSM="$prefix" RGPU=$runbarcode SORT_ORDER=coordinate CREATE_INDEX=TRUE

## merge se and pe bam files with samtools and index
samtools index bam/"$prefix".bam
samtools sort bam/"$prefix".bam -o bam/"$prefix".sorted.bam

done < prefix.silu.p3.list

#step 4 SNP calling
gatk='GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar'
while read prefix
do
java -jar $gatk \
-nct 26 \
-T HaplotypeCaller \
-R $ref \
-I bam/"$prefix".sorted.bam \
-ERC GVCF \
-o gvcf/"$prefix".gvcf
done < prefix.silu.p1.list


#step5: combo.vcf
gatk='GBS/tools/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar'
log="GBS/gbs.2020.ref.dec.2020/log"
gvcf="GBS/gbs.2020.ref.dec.2020/gvcf"
ref="GBS/gbs.2020.ref.dec.2020/ref/mywagenomev2.1.fa"

##run this with headers in a combo.gvcftovcf.sh
while read prefix
 do
 java -jar $gatk \
 -nct 26 \
 -T HaplotypeCaller \
 -R $ref \
 -I bam/"$prefix".sorted.bam \
 -log log/"$prefix".SNPcalling.log \
 --emitRefConfidence GVCF \
 --fix_misencoded_quality_scores \
 -variant_index_type LINEAR \
 -variant_index_parameter 128000 \ 
 -o gvcf/"$prefix".gvcf \
 done < prefix.list


#ste6: genptype to GVCF 
gatk='/data/chend8/GBS/tools/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar'
ref="/data/chend8/GBS/gbs.2020.ref.dec.2020/ref/mywagenomev2.1.fa"
java -jar $gatk \
-T GenotypeGVCFs \
-nt 10 \
-R $ref \
-V gvcf.silu.maddie.list \
--disable_auto_index_creation_and_locking_when_reading_rods \
-o silu.maddie.vcf \
-log gvcf.intoVCF.log



#filtering VCF
gatk='/data/chend8/GBS/tools/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar'
ref="/data/chend8/GBS/gbs.2020.ref.dec.2020/ref/mywagenomev2.1.fa"
vcftools --vcf silu.maddie.vcf --remove-indels --recode --recode-INFO-all --out silu.maddie.indels
vcftools --vcf silu.maddie.indels.recode.vcf --max-alleles 4 --recode --recode-INFO-all --out silu.maddie.indel.maxA4
vcftools --vcf silu.maddie.indel.maxA4.recode.vcf --minGQ 20 --recode --recode-INFO-all --out silu.maddie.indel.maxA4.gq
vcftools --vcf silu.maddie.indel.maxA4.recode.vcf --maf 0.05 --recode --recode-INFO-all --out  silu.maddie.indel.maxA4.gp.maf
vcftools --vcf silu.maddie.indel.maxA4.gp.maf.recode.vcf --max-missing 0.7 --recode --recode-INFO-all --out silu.maddie.indel.maxA4.gp.maf.maxmiss0.7
vcftools --vcf silu.maddie.vcf --minGQ 20 --recode --recode-INFO-all --out silu.maddie.minGQ20
#calculate FST
#1, between inland TOWA and HEWA
vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --weir-fst-pop ../lists/HEWA.list --weir-fst-pop ../lists/inlandTOWA.list --out hewa.inlandtowa

#2, between inland TOWA and coastalTOWA.nonOutlier.list
 vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --weir-fst-pop ../lists/ocoastalTOWA.list --weir-fst-pop ../lists/inlandTOWA.list --out ocoastalTOWA.inlandtowa

#3, between inland TOWA and valdez
 vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --weir-fst-pop ../lists/valdez.list --weir-fst-pop ../lists/inlandTOWA.list --out valdez.inlandtowa

#4, between inland TOWA and haida gwaii
 vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --weir-fst-pop ../lists/haidagwaii.list --weir-fst-pop ../lists/inlandTOWA.list --out haidagwaii.inlandtowa
#5, between HEWA and coastalTOWA.nonOutlier.list
 vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --weir-fst-pop ../lists/ocoastalTOWA.list --weir-fst-pop ../lists/HEWA.list --out ocoastalTOWA.hewa
#6, between HEWA and valdez
 vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --weir-fst-pop ../lists/valdez.list --weir-fst-pop ../lists/HEWA.list --out valdez.hewa
#7, between HEWA and haida gwaii
 vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --weir-fst-pop ../lists/haidagwaii.list --weir-fst-pop ../lists/HEWA.list --out haidagwaii.hewa
#8, between valdez and haida gwaii
 vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --weir-fst-pop ../lists/haidagwaii.list --weir-fst-pop ../lists/valdez.list --out haidagwaii.valdez
#9, between haida gwaii and oc
 vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --weir-fst-pop ../lists/haidagwaii.list --weir-fst-pop ../lists/ocoastalTOWA.list --out haidagwaii.ocoastalTOWA
#10, between valdez and oc
 vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --weir-fst-pop ../lists/valdez.list  --weir-fst-pop ../lists/ocoastalTOWA.list --out valdez.ocoastalTOWA

#make vcf file with only these individuals 
 vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --keep ../lists/hewa.citowa.197.sorted.list --recode --recode-INFO-all --out  silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.197indv

#calculate allele frequency for HEWA
vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --keep ../lists/HEWA.list --freq  --out hewa
#calculate allele frequency for inland TOWA
vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --keep ../lists/inlandTOWA.list --freq --out inlandtowa
#calculate allele frequency for haida gwaii
vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --keep ../lists/haidagwaii.list --freq --out haidagwaii
##calculate allele frequency for  valdez
vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.recode.vcf --keep ../lists/valdez.list --freq --out valdez


#make input for GENABLE
#turn VCF into plink
vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.197indv.recode.vcf  --plink-tped --out silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.197indv

#extract positions along chr5, 1a, 4a, z
grep -F -w -f pattern.file mywagenomev2.1.all.noseq.gff > chr5.genes #extract chromosome 5 genes
grep -F -w -f pattern.file.chr5.exon.chr4a mywagenomev2.1.all.noseq.gff > chr5.exon #extract chromosome 5 genes

grep -F -w -f pattern.file.chrz mywagenomev2.1.all.noseq.gff > chrz.genes #extract chromosome z genes
grep -F -w -f pattern.file.chr1a mywagenomev2.1.all.noseq.gff > chr1a.genes #extract chromosome 1a genes
grep -F -w -f pattern.file.chr4a mywagenomev2.1.all.noseq.gff > chr4a.genes #extract chromosome 4a genes
##extract sequences for these candidate genes
#chr5
samtools faidx mywagenomev2.1.fa chr5:3869384-5010857 > chr5.peak.g1
samtools faidx mywagenomev2.1.fa chr5:38813528-38819339 > chr5.peak2.g1
samtools faidx mywagenomev2.1.fa chr5:38824164-38840412 > chr5.peak2.g2
#chrz mt association
samtools faidx mywagenomev2.1.fa chrz:11387594-11463763 > chrz.peak1.g1
samtools faidx mywagenomev2.1.fa chrz:11791989-12057618 > chrz.peak1.g2
samtools faidx mywagenomev2.1.fa chrz:16796558-16824247 > chrz.peak2.g1
samtools faidx mywagenomev2.1.fa chrz:16852943-16854937 > chrz.peak2.g2
samtools faidx mywagenomev2.1.fa chrz:16833380-16839300 > chrz.peak2.g3
samtools faidx mywagenomev2.1.fa chrz:63745467-63774490 > chrz.peak3.g1
samtools faidx mywagenomev2.1.fa chrz:63779401-63788384 > chrz.peak3.g2
$chrz fst peaks
samtools faidx mywagenomev2.1.fa chrz:50100432-50145287 > chrz.fst1.g1
samtools faidx mywagenomev2.1.fa chrz:50146181-50187539 > chrz.fst1.g2
samtools faidx mywagenomev2.1.fa chrz:57089559-57091548 > chrz.fst2.g1
samtools faidx mywagenomev2.1.fa chrz:57216097-57218574 > chrz.fst2.g2

#######Nature Comm reviews--> calculate diversity and dxy 
#PI 
vcftools --vcf silu.maddie.indels.recode.vcf --max-alleles 2 --minGQ 20 --max-missing 0.7  --recode --recode-INFO-all --out silu.maddie.indels.maxA2.minQ20.maxmiss0.7
#1, pi SOCC
vcftools --vcf silu.maddie.indels.maxA2.minQ20.maxmiss0.7.recode.vcf --keep ../lists/HEWA.list --site-pi  --out socc.site
#2, pi inland STOW
vcftools --vcf silu.maddie.indels.maxA2.minQ20.maxmiss0.7.recode.vcf --keep ../lists/inlandTOWA.list --site-pi  --out stow.site
#3, pi valdez
vcftools --vcf silu.maddie.indels.maxA2.minQ20.maxmiss0.7.recode.vcf --keep ../lists/valdez.list --site-pi  --out valdez.site
#4, pi haidagwaii
vcftools --vcf silu.maddie.indels.maxA2.minQ20.maxmiss0.7.recode.vcf --keep ../lists/haidagwaii.list --site-pi  --out haidagwaii.site
#5, pi o-stow
vcftools --vcf silu.maddie.indels.maxA2.minQ20.maxmiss0.7.recode.vcf --keep ../lists/ocoastalTOWA.list --site-pi  --out oc.stow.site

#PI window PI size=100bp
#1, pi SOCC
vcftools --vcf silu.maddie.indels.maxA2.minQ20.maxmiss0.7.recode.vcf --keep ../lists/HEWA.list  --window-pi 100  --out socc #--chr chr5
#2, pi inland STOW
vcftools --vcf silu.maddie.indels.maxA2.minQ20.maxmiss0.7.recode.vcf --keep ../lists/inlandTOWA.list --window-pi 100   --out stow
#3, pi valdez
vcftools --vcf silu.maddie.indels.maxA2.minQ20.maxmiss0.7.recode.vcf --keep ../lists/valdez.list  --window-pi 100   --out valdez
#4, pi haidagwaii
vcftools --vcf silu.maddie.indels.maxA2.minQ20.maxmiss0.7.recode.vcf --keep ../lists/haidagwaii.list --window-pi 100   --out haidagwaii
#5, pi oc-stow
vcftools --vcf silu.maddie.indels.maxA2.minQ20.maxmiss0.7.recode.vcf --keep ../lists/ocoastalTOWA.list --window-pi 100   --out oc.stow
#6, pi for all coastal stow
cat haidagwaii.list valdez.list ocoastalTOWA.list > coastalSTOW.list
vcftools --vcf silu.maddie.indels.maxA2.minQ20.maxmiss0.7.recode.vcf --keep ../lists/coastalSTOW.list --window-pi 100   --out cstow


#LD --ld-window-bp <integer>
#1, LD valdez
vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.197indv.recode.vcf --keep ../lists/valdez.list --chr chr5 --from-bp 2500000 --to-bp 7000000 --geno-r2 --out valdez.ch5.r2
#2, LD haidagwaii
vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.197indv.recode.vcf --keep ../lists/haidagwaii.list --chr chr5   --from-bp 2500000 --to-bp 7000000 --geno-r2  --out haidagwaii.ld
#3, LD o-stow
vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.197indv.recode.vcf --keep ../lists/ocoastalTOWA.list  --chr chr5  --from-bp 2500000 --to-bp 7000000 --geno-r2  --out oc.stow.ld
#4, all coastal STOW
vcftools --vcf silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.197indv.recode.vcf --keep ../lists/coastalSTOW.list  --chr chr5  --from-bp 2500000 --to-bp 7000000 --geno-r2  --out c.stow.ld


#prepare filtration VCF that only include coastal stow and run PIXY
#https://pixy.readthedocs.io/en/1.0.0.beta1/guide/pixy_guide.html#site-level-filtration
#step1.1: keep invariant sites vcf
#make a vcf file including variants
#combine gvcf 
gatk='GBS/tools/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar'
ref="GBS/gbs.2020.ref.dec.2020/ref/mywagenomev2.1.fa"
java -jar $gatk \
-T GenotypeGVCFs \
-nt 10 \
-R $ref \
-V gvcf.silu.maddie.list \
--include-non-variant-sites \
--disable_auto_index_creation_and_locking_when_reading_rods \
-o silu.maddie.invariant.vcf \
-log gvcf.intoVCF.log
#step1.2: filter vcf
#filter
vcftools --vcf silu.maddie.invariant.recode.vcf --remove-indels --max-alleles 2 --minGQ 20 --max-missing 0.7  --recode --recode-INFO-all --out silu.maddie.indels.maxA2.minQ20.maxmiss0.7
vcftools --vcf silu.maddie.indels.maxA2.minQ20.maxmiss0.7.recode.vcf --keep ../lists/coastalSTOW.list --recode --recode-INFO-all --out silu.maddie.indels.maxA2.minQ20.maxmiss0.7.all.cSTOW
#here is the output:    silu.maddie.indels.maxA2.minQ20.maxmiss0.7.all.cSTOW.recode.vcf
#step2: make population file
n=read.table("coastalSTOW.list",header=F)
p=c(rep("haidagwaii", 10),rep("valdez", 15),rep("ocstow", 33))
sink("pop.cstow.sub.list"); for(i in 1:length(n[,1])){cat(n[i,1], "\t", p[i], "\n", sep="")}; sink()
sink("pop.cstow.list"); for(i in 1:length(n[,1])){cat(n[i,1], "\t", "cstow", "\n", sep="")}; sink()

t=read.table("socc.stow.cstow.list",header=F)
p=c(rep("socc",62),rep("istow",77),rep("haidagwaii", 10),rep("valdez", 15),rep("ocstow", 33))
sink("pop.socc.istow.cstow.sub.196.list"); for(i in 1:length(t[,1])){cat(t[i,1], "\t", p[i], "\n", sep="")}; sink()
pc=c(rep("socc",62),rep("istow",77),rep("cstow",58))
sink("pop.socc.istow.cstow.all.196.list"); for(i in 1:length(t[,1])){cat(t[i,1], "\t", pc[i], "\n", sep="")}; sink()

#this one include subpop: pop.cstow.sub.list
#step3: pixy
#step3.1: prepare input for pixy
#step3.1.1: GenomicsDBImport
gatk='/GBS/tools/gatk-4.0.12.0/gatk-package-4.0.12.0-local.jar'
ref="/GBS/gbs.2020.ref.dec.2020/ref/mywagenomev2.1.fa"
java -jar $gatk \
       GenomicsDBImport \
       --genomicsdb-workspace-path my_database_chr5 \
       --batch-size 50 \
       -L chr5 \
       --sample-name-map silu.maddie.map \
       --reader-threads 10

#step3.1.2: GenotypeGVCFs:
#pipeline over here --> https://pixy.readthedocs.io/en/1.0.0.beta1/generating_invar/generating_invar.html
gatk='/GBS/tools/gatk-4.0.12.0/gatk-package-4.0.12.0-local.jar'
ref="/GBS/gbs.2020.ref.dec.2020/ref/mywagenomev2.1.fa"
java -jar $gatk GenotypeGVCFs \
-R $ref -V gendb://my_database_chr5 \
-all-sites -L chr5 -O pixy.chr5.vcf.gz

#step3.1.3: filtering:
vcftools --gzvcf pixy.chr5.vcf.gz \
--remove-indels \
--max-missing 0.7 \
--minQ 20 \
--recode --stdout | gzip -c > pixy.chr5.filtered.vcf.gz

#move pixy.chr5.vcf.gz to working directory
#step3.2 run pixy 
conda create --name pixy python=3.8 numpy=1.20.2
conda activate pixy
conda install --yes -c conda-forge pixy
conda install --yes -c bioconda htslib
bgzip pixy.chr5.filtered.vcf  
tabix  pixy.chr5.filtered.vcf.gz
#or 
bcftools index  pixy.chr5.filtered.vcf.gz

pixy --stats pi \
--vcf pixy.chr5.filtered.vcf.gz \
--populations pop.cstow.list \
--window_size 1000 \
--n_cores 5 \
--chromosomes 'chr5'

#run with sub populations
pixy --stats pi fst dxy \
--vcf pixy.chr5.filtered.vcf.gz \
--populations pop.cstow.sub.list \
--window_size 1000 \
--n_cores 5 \
--chromosomes 'chr5'

conda create --prefix /GBS/pixy python=3.8 numpy=1.20.2
