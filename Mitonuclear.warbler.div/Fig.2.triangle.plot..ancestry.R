library(GenABEL)
library(SNPRelate)
library(gdsfmt)
library(ResourceSelection)
#source("https://bioconductor.org/biocLite.R")
#biocLite("GWASTools")
library("GWASTools")
library(plotrix)


#d=read.csv("/Users/siluwang/Desktop/1_PhD.projects/1.0.CoastalTOWA/1.Oct.coastal.towa.2019.MS.version/Ref.yrwa.Dec.2020/fst.peaks.hewa.itowa.csv")
d=read.csv("fst.0.3peaks.hewa.itowa.csv")
g=snpgdsOpen("silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.197indv.recode.gds",  readonly = TRUE)


allele <- read.gdsn(index.gdsn(g, "snp.allele"))
genome <- read.gdsn(index.gdsn(g, "genotype"))
id=read.gdsn(index.gdsn(g, "sample.id"))

chr1.cand=d$POS[which(d$CHROM=="chr1")]
chr1a.cand=d$POS[which(d$CHROM=="chr1a")]
chr4a.cand=d$POS[which(d$CHROM=="chr4a")]
chr5.cand=d$POS[which(d$CHROM=="chr5")]
chr20.cand=d$POS[which(d$CHROM=="chr20")]
chrz.cand=d$POS[which(d$CHROM=="chrz")]

chr10.cand=d$POS[which(d$CHROM=="chr10")]
chr11.cand=d$POS[which(d$CHROM=="chr11")]
chr12.cand=d$POS[which(d$CHROM=="chr12")]
chr13.cand=d$POS[which(d$CHROM=="chr13")]
chr14.cand=d$POS[which(d$CHROM=="chr14")]
chr15.cand=d$POS[which(d$CHROM=="chr15")]
chr17.cand=d$POS[which(d$CHROM=="chr17")]
chr18.cand=d$POS[which(d$CHROM=="chr18")]
chr19.cand=d$POS[which(d$CHROM=="chr19")]
chr.cand=d$POS[which(d$CHROM=="chr10")]


#function to find genotype for each canddiate snp
genotype.chr.pos=function(chr, pos)
{snps={} #find columns of chr5 candidate SNPs
for(i in 1:length(pos))
{col=genome[, intersect(which(snpgdsSNPList(g)$chrom==chr), which(snpgdsSNPList(g)$position==pos[i]))]
	snps=cbind(snps, col)}
colnames(snps)=as.character(paste(chr, pos, sep="."))
return(snps)}

chr1d=genotype.chr.pos(chr=1, pos=chr1.cand)
chr1ad=genotype.chr.pos(chr="1a", pos=chr1a.cand)
chr4ad=genotype.chr.pos(chr="4a", pos=chr4a.cand)
chr5d=genotype.chr.pos(chr=5, pos=chr5.cand)
chr20d=genotype.chr.pos(chr=20, pos=chr20.cand)
chrzd=genotype.chr.pos(chr="z", pos=chrz.cand)


###match with saved
saved=read.csv("silu.maddie.197.genomicEV.mt.csv")
#for each SNP, find the HEWA and TOWA version, and change the scorign so that HEWA=0, TOWA=1, and heterozygotes =0.5
dd=data.frame(chr1d, chr1ad, chr4ad, chr5d, chr20d, chrzd, id)
table(saved$pop.ids) 
for(j in 1:(ncol(dd)-1))
	{dd[which(dd[,j]==3), j]=NA #change the missing values to NA
	dd[which(dd[,j]==1), j]=0.5 #heterozygotes (1) replaced by 0.5
	hewa=dd[which(saved$pop.ids=="HEWA"),j] #take out HEWA rows
	itowa=dd[which(saved$pop.ids=="inlandTOWA"),j] #take out i-TOWA rows
	 cat(" snp", j, "\n")
	cat("hewa:", mean(hewa, na.rm=T), "\n"); 
	cat("itowa:", mean(itowa, na.rm=T), "\n");
	if(mean(hewa, na.rm=T)>1.2 && mean(itowa, na.rm=T) < 0.9) #in situation hewa has the 2 genotype and towa has 0
		{	
		rowsto0=which(dd[,j]==2);  #these rows currently is 2 shoud be 0
		rowsto1=which(dd[,j]==0);  #thrsese rows currently is 0 should be 1
		dd[rowsto0, j]=0 #assign
		dd[rowsto1, j]=1
			} else if(mean(hewa, na.rm=T)<0.9 && mean(itowa, na.rm=T) > 1.2) #case when hewa is 0 and towa is 2
				{rowsto1=which(dd[,j]==2); dd[rowsto1, j]=1;  #simply change towa to 1
				} else {print("ERROR: CHANGE CUTOFF~~~~~~~~~~");}}

#check if HEWA is mostly 0 
for(j in 1:(ncol(dd)-1))	
{	hewa=dd[which(saved$pop.ids=="HEWA"),j]
	itowa=dd[which(saved$pop.ids=="inlandTOWA"),j]
	print(j)
	print(hewa); cat("\n")}

ddfull=data.frame(dd, saved)
ddfull=ddfull[order(ddfull$pop.SP, ddfull$mt.hap),] #careful!!! changed orders for samples by population
##check  that the ddfull matches between genotype info and background info
sum(ddfull$id==ddfull$pca.sample.id)==nrow(ddfull)
#assign sub-pops
mtd=ddfull[which(is.na(ddfull$mt.hap)==0),]
haida=mtd[which(mtd$pop.SP=="hadiagwaii"),]
val=mtd[which(mtd$pop.SP=="valdez"),]
octowa=mtd[intersect(which(mtd$pop.SP=="coastalTOWA"),which(is.na(rowMeans(mtd[,1:50], na.rm=T))==0)),]; 
itowa=mtd[which(mtd$pop.SP=="inlandTOWA"),]
hewa=mtd[c(which(mtd$pop.SP=="inlandHEWA"),which(mtd$pop.SP=="coastalHEWA")),]

#funciton to plot heatmap for populations
heatmap.popd=function(popd){
	pplot=popd[order(popd$mt.hap),]; plot=pplot[,c(which(colnames(pplot)=="mt.hap"),1:(ncol(dd)-1))]
	rownames(plot)=pplot$sample.id
	 color2D.matplot(plot,  extremes = c("turquoise", "plum"))
	 abline(v=c(1, (1+length(chr1.cand)),(1+length(chr1.cand)+length(chr1a.cand)),(1+length(chr1.cand)+length(chr1a.cand)+length(chr4a.cand)),(1+length(chr1.cand)+length(chr1a.cand)+length(chr4a.cand)+length(chr5.cand)),(1+length(chr1.cand)+length(chr1a.cand)+length(chr4a.cand)+length(chr5.cand)+length(chr20.cand)), (1+length(chr1.cand)+length(chr1a.cand)+length(chr4a.cand)+length(chr5.cand)+length(chr20.cand)+length(chrz.cand))), lwd=2, lty=2)	}
heatmap.popd(val)
heatmap.popd(haida)
heatmap.popd(octowa)
heatmap.popd(itowa)
heatmap.popd(hewa)
#figure 2 HEWA | c-TOWA | i-towa
set.seed(26);hewa.srow=sample(1:nrow(hewa), 10)
set.seed(26);itowa.srow=sample(1:nrow(itowa), 10)
hewa.sp=hewa[hewa.srow,]; itowa.sp=itowa[itowa.srow,]
octowa.par=rbind(hewa.sp, octowa, itowa.sp)
heatmap.popd(octowa.par)


##########FIGURRE 2B #make HEATMAP SOCC, I-STOW, C-STOW
hewa.rows=c(which(mtd$pop.SP=="inlandHEWA"),which(mtd$pop.SP=="coastalHEWA")); itowa.rows=which(mtd$pop.SP=="inlandTOWA")
set.seed(23);hewa.srow=sample(hewa.rows, 10)
set.seed(23);itowa.srow=sample(itowa.rows, 10)
ctowa=mtd[c(hewa.srow,which(mtd$pop.SP=="hadiagwaii"),intersect(which(mtd$pop.SP=="coastalTOWA"),which(is.na(rowMeans(mtd[,1:50], na.rm=T))==0)),which(mtd$pop.SP=="valdez"),itowa.srow),]
ctowa.plot=ctowa[,c(which(colnames(ctowa)=="mt.hap"),1:(ncol(dd)-1))]
 color2D.matplot(ctowa.plot,  extremes = c("turquoise", "plum"))
 abline(v=c(1, (1+length(chr1.cand)),(1+length(chr1.cand)+length(chr1a.cand)),(1+length(chr1.cand)+length(chr1a.cand)+length(chr4a.cand)),(1+length(chr1.cand)+length(chr1a.cand)+length(chr4a.cand)+length(chr5.cand)),(1+length(chr1.cand)+length(chr1a.cand)+length(chr4a.cand)+length(chr5.cand)+length(chr20.cand)), (1+length(chr1.cand)+length(chr1a.cand)+length(chr4a.cand)+length(chr5.cand)+length(chr20.cand)+length(chrz.cand))), lwd=3, lty=1)	
abline(h=(nrow(ctowa)-10), lwd=3)
abline(h=(nrow(ctowa)-nrow(haida)-10), lwd=3)
abline(h=(nrow(ctowa)-nrow(haida)-nrow(octowa)-10), lwd=3 )
abline(h=(nrow(ctowa)-nrow(haida)-nrow(octowa)-nrow(val)-10), lwd=3 )
###############


ctowa.n=ctowa[,1:50]
ctowa.h=ctowa[,1:50]; #make heterozygosity matrix, first assign ancestry scores for it, then update the columns
ctowa$pop.SP=droplevels(ctowa$pop.SP); levels(ctowa$pop.SP)
for(j in 1:ncol(ctowa.h))
{to0=c(which(ctowa.n[,j]==1),which(ctowa.n[,j]==0)) #homozygous sites to 0
to1=which(ctowa.n[,j]==0.5)#homozygous sites to 1
ctowa.h[to0,j]=0;ctowa.h[to1,j]=1;}
unique(ctowa.h, na.rm=T)#should only have 0 and 1
hi.50snps=rowMeans(ctowa.n, na.rm=T)
het.50snps=rowMeans(ctowa.h, na.rm=T)
ctowa.col=c("gold","seagreen","royalblue")#sequence accroding to levels(ctowa$pop.SP)
rgb.convert=function(col.vect,fr)
{cols={}
for(i in col.vect)
	{cl=col2rgb(i)/255; cols=c(cols, rgb(cl[1,1],cl[2,1],cl[3,1], fr))}; return(cols)}
plot(hi.50snps, het.50snps, xlim=c(0,1), ylim=c(0,1), pch=21, bg=rgb.convert(ctowa.col[as.numeric(ctowa$pop.SP)],0.5), xlab="Admixture Proportion", ylab="Heterozygosity")
abline(0,2);abline(2,-2)

##calculate HI and HET for parentals
par=rbind(hewa, itowa)
par.n=par[,1:50]
par.h=par[,1:50]; #make heterozygosity matrix, first assign ancestry scores for it, then update the columns
par$pop.ids=droplevels(par$pop.ids); levels(par$pop.ids)
for(j in 1:ncol(par.h))
{to0=c(which(par.n[,j]==1),which(par.n[,j]==0)) #homozygous sites to 0
to1=which(par.n[,j]==0.5)#homozygous sites to 1
par.h[to0,j]=0;par.h[to1,j]=1;}
unique(par.h, na.rm=T)#should only have 0 and 1
hi.50snps.par=rowMeans(par.n, na.rm=T)
het.50snps.par=rowMeans(par.h, na.rm=T)

##########FIGURRE 2 A: TRIANGLE PLOT 
par.col=c("turquoise", "plum")
plot(hi.50snps, het.50snps, xlim=c(0,1), ylim=c(0,1), pch=21,bg=rgb.convert(ctowa.col[as.numeric(ctowa$pop.SP)],0.8),  xlab="Admixture Proportion", ylab="Heterozygosity", cex=1.3)
points(hi.50snps.par, het.50snps.par, pch=16, col=rgb.convert(par.col[as.numeric(par$pop.ids)],0.7))
points(hi.50snps, het.50snps, xlim=c(0,1), ylim=c(0,1), pch=21,bg=rgb.convert(ctowa.col[as.numeric(ctowa$pop.SP)],0.8), cex=1.3)
abline(0,2);abline(2,-2)
legend("topleft", c("SOCC", "i-STOW"), pch=16, col=par.col , cex=0.8, pt.cex=1)
legend("topright", c( "Haida Gwaii", "Valdez", "oc-STOW"), pt.cex=1.3 ,pch=21,pt.bg=rgb.convert(c("seagreen","royalblue", "gold" ), 0.8), cex=0.7 )

table(par$pop.ids)
table(ctowa$pop.SP)
