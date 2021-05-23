library(GenABEL)
library(SNPRelate)
library(gdsfmt)
library(ResourceSelection)
library("GWASTools")
########**********************************************
#ph=read.table("plate1-4_withpopinfo_distr_7plumage_age.cor.transform.txt", header=T)
#STEP1: make 
saved=read.csv("silu.maddie.197.genomicEV.mt.csv") #mt=read.csv("~/Desktop/1_PhD.projects/1.0.CoastalTOWA/mt.pop.226indv.towa.hewa.burke.Jan2019.csv")
saved$sex=1
d=saved[,c(1, which(colnames(saved)=="sex"),2:ncol(saved))]
colnames(d)[1]="id"
#write.table(d, "silu.maddie.197.genomicEV.mt.jan.2021.txt", row.names=F)
########**********************************************
#convert.snp.tped(tpedfile="silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.197indv.tped",tfamfile="silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.197indv.tfam", "silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.197indv.plink.raw",strand = "+")
##
p=read.table("silu.maddie.197.genomicEV.mt.txt", header=1)
#tohe=load.gwaa.data(phe="plate1-4_withpopinfo_distr_7plumage_age.cor.transform.txt",gen="tohe.p1234.0.7miss.plink.raw",force=TRUE) 
tohe=load.gwaa.data(phe="silu.maddie.197.genomicEV.mt.jan.2021.txt",gen="silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.197indv.plink.raw",force=TRUE) 
#tohe=load.gwaa.data(phe="plate1-4_withpopinfo_distr_8plumage_agg.transform.cheekHZ.txt",gen="tohe.p1234.0.7miss.plink.raw",force=TRUE)
#tohe=load.gwaa.data(phe="plate1-4_withpopinfo_distr_8plumage_agg.transform.txt",gen="TOHE.plate1234.rmOutGL.qual.indels.gq.maf.0.7miss.plink.raw",force=TRUE)#sep, 2019
class(tohe)
chromosome(tohe)
nids(tohe) #Number of individuals
nsnps(tohe) #Number of SNPS
#coding(tohe) # The reference and alternative allele at each SNP.
#strand(tohe) # Strand for each SNP: set to "+" for all SNPs
descriptives.marker(tohe)
summary(tohe)# Summary of data: see page 53 of tutorial for an explination of the columns
head(summary(tohe)) 
#summary(tohe)$Position
descriptives.trait(tohe) # Statistics on phenotypic traits
names(phdata(tohe))

tohe.gkin=ibs(tohe,weight="freq")  
replace=hom(tohe)$Var
diag(tohe.gkin)=replace
library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")
image(tohe.gkin,col=pal) 

#RUN MT ASSOCAITION TEST
mt.cor=egscore(mt.hap,data=tohe,kin=tohe.gkin, clambda =F) 
summary(mt.cor)
lambda(mt.cor)
estlambda(mt.cor[, "P1df"], plot=TRUE, pch=16 ) 
plot(mt.cor, df="Pc1df",cex=1.3,col=c("coral ", "gold"),pch=16, main="Association P-values mtDNA haplotype")
abline(h=6, col="red")
descriptives.scan(mt.cor, sort="Pc1df") 
result_mt=results(mt.cor); 
result_mt.outlier=subset(result_mt,Pc1df<=1e-6)
result_mt.outlier
write.csv(result_mt.outlier, "silu.maddie.197.mt.10e-6.pop.cor.csv")

#######extract genotypes for candidate SNPs 
g=snpgdsOpen("silu.maddie.indel.maxA4.gp.maf.maxmiss0.7.197indv.recode.gds",  readonly = TRUE)

allele <- read.gdsn(index.gdsn(g, "snp.allele"))
genome <- read.gdsn(index.gdsn(g, "genotype"))
id=read.gdsn(index.gdsn(g, "sample.id"))
#all FST> 0.6 peaks on chr5: 3880486 3969384 4007220 4364934 4397756 4400734 4468249 4519381 4621259 4675788 4729930 4806290 4856848 4898166 5010857
chr5.cand=c(3880486,3969384, 4007220, 4364934, 4397756, 4400734, 4468249, 4519381, 4621259, 4675788, 4729930, 4806290 ,4856848 ,4898166, 5010857)
chr5.snp=4364934
chr5.snps={} #find columns of chr5 candidate SNPs
for(i in 1:length(chr5.cand))
{col=genome[, intersect(which(snpgdsSNPList(g)$chrom=="5"), which(snpgdsSNPList(g)$position==chr5.cand[i]))]
	chr5.snps=cbind(chr5.snps, col)}
dim(chr5.snps)
colnames(chr5.snps)=paste("chr5.snp", 1:length(chr5.cand), sep="")

###match with saved
saved=read.csv("silu.maddie.197.genomicEV.mt.csv")
#for each SNP, find the HEWA and TOWA version, and change the scorign so that HEWA=0, TOWA=1, and heterozygotes =0.5
dd=data.frame(chr5.snps,chrz.snps, id)
table(saved$pop.ids) 
for(j in 1:length(c(chr5.cand, chrz.cand)))
	{dd[which(dd[,j]==3), j]=NA
	dd[which(dd[,j]==1), j]=0.5
	hewa=dd[which(saved$pop.ids=="HEWA"),j]
	itowa=dd[which(saved$pop.ids=="inlandTOWA"),j]
	 cat(" snp", j, "\n")
	cat("hewa:", mean(hewa, na.rm=T), "\n"); 
	cat("itowa:", mean(itowa, na.rm=T), "\n");
	if(mean(hewa, na.rm=T)>1.3 && mean(itowa, na.rm=T) < 0.7)
		{#cat("TO change, HEWA=2", " snp", j, "\n")	; print(hewa)		
		rowsto0=which(dd[,j]==2); #cat("rowsto0: ",dd[rowsto0, j], "\n")
		rowsto1=which(dd[,j]==0); 
		dd[rowsto0, j]=0
		dd[rowsto1, j]=1
			#cat("after",dd[rowsto0, j], "\n")
			#cat("hewa.after",dd[which(saved$pop.ids=="HEWA"),j],"\n")
			} else if(mean(hewa, na.rm=T)<0.7 && mean(itowa, na.rm=T) > 1.3){rowsto1=which(dd[,j]==2); dd[rowsto1, j]=1; 
				} }
			

#check
for(j in 1:length(c(chr5.cand, chrz.cand)))	
{	hewa=dd[which(saved$pop.ids=="HEWA"),j]
	itowa=dd[which(saved$pop.ids=="inlandTOWA"),j]
	print(j)
	print(hewa); cat("\n")}

library(readxl)
library(reshape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DescTools)
#library(fpc)
library("RColorBrewer")
library(qlcMatrix)
library(splitstackshape)
col.sp <- colorRampPalette(c("turquoise","plum"))(3)
#jpeg(filename="heatmap.jpeg",width=8,height=8,units="in",res=500)
ddfull=data.frame(dd, saved)
ddfull=ddfull[order(ddfull$pop.SP, ddfull$mt.hap),] #careful!!! changed orders for samples by population
#write.csv(ddfull,"chr5.chr.z.cand.hapblock.heatmap.input.csv")

dd.plot=ddfull[c(1:23, 66:104, 56:65,24:55, 182:196,105:181),c(1:length(c(chr5.cand, chrz.cand)), which(colnames(ddfull)=="mt.hap"))]
rownames(dd.plot)=paste(ddfull$pop.SP[c(1:23, 66:104, 56:65,24:55, 182:196,105:181)], ddfull$id[c(1:23, 66:104, 56:65,24:55, 182:196,105:181)], sep=".")
#mar=c(1,1,1,1)
heatmap(as.matrix(dd.plot) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.2)
#dev.off()
mtd=ddfull[which(is.na(ddfull$mt.hap)==0),]
haida=mtd[which(mtd$pop.SP=="hadiagwaii"),]
val=mtd[which(mtd$pop.SP=="valdez"),]
ctowa=mtd[which(mtd$pop.SP=="coastalTOWA"),]
itowa=mtd[which(mtd$pop.SP=="inlandTOWA"),]
hewa=mtd[c(which(mtd$pop.SP=="inlandHEWA"),which(mtd$pop.SP=="coastalHEWA")),]
library(plotrix)
heatmap.popd=function(popd){
	pplot=popd[order(popd$mt.hap),]; plot=pplot[,c(which(colnames(pplot)=="mt.hap"),1:length(chr5.cand))]
	rownames(plot)=pplot$sample.id
	 color2D.matplot(plot,  extremes = c("turquoise", "plum"))	}
heatmap.popd(val)
heatmap.popd(haida)
heatmap.popd(ctowa)
heatmap.popd(itowa)
heatmap.popd(hewa)
all=rbind(hewa,  haida, ctowa, val,itowa)
allp=all[,c(which(colnames(all)=="mt.hap"),1:length(chr5.cand))]
color2D.matplot(allp,  extremes = c("turquoise", "plum"))
  abline(h=nrow(all)-nrow(hewa), col="black", lwd=2)
  abline(h=nrow(all)-(nrow(hewa)+nrow(haida)), col="black",lwd=2)
  abline(h=nrow(all)-(nrow(hewa)+nrow(ctowa)+nrow(haida)), col="black",lwd=2)
  abline(h=nrow(all)-(nrow(hewa)+nrow(ctowa)+nrow(haida)+nrow(val)), col="black",lwd=2)

#only plot 10 randomly sampled HEWA, i-TOWA
set.seed(26);itwoa.srows=sample(1:nrow(itowa), 10)
set.seed(26);hewa.srows=sample(1:nrow(hewa), 10)
all.ps=rbind(hewa[hewa.srows,],  haida, ctowa, val,itowa[itwoa.srows,])
all.ps.p=all.ps[,c(which(colnames(all.ps)=="mt.hap"),1:length(chr5.cand))]
color2D.matplot(all.ps.p,  extremes = c("turquoise", "plum"))
  abline(h=(nrow(all.ps.p)-10), col="black", lwd=2)
  abline(h=nrow(all.ps.p)-(10+nrow(haida)), col="black",lwd=2)
  abline(h=nrow(all.ps.p)-(10+nrow(ctowa)+nrow(haida)), col="black",lwd=2)
  abline(h=nrow(all.ps.p)-(10+nrow(ctowa)+nrow(haida)+nrow(val)), col="black",lwd=2)

  
ddfulls=ddfull[c(1:23, 66:104, 56:65,24:55, 182:196,105:181),]
dd.sortmt=ddfulls[order(ddfulls$mt.hap),] 
dd.plot.mt=dd.sortmt[,c(1:length(c(chr5.cand, chrz.cand)), which(colnames(dd.sortmt)=="mt.hap"))]
rownames(dd.plot.mt)=paste(dd.sortmt$pop.SP, dd.sortmt$id, sep=".")
#mar=c(1,1,1,1)
heatmap(as.matrix(dd.plot.mt) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.2)
#dev.off()
#check if chr5 ancestry is correlated with chr z ancestry in differnat pops
head(ddfull)
chr5g=ddfull[,1:14]; ddfull$chr5.hi=rowMeans(chr5g, na.rm=T)
chrzg=ddfull[,15:17]; ddfull$chrz.hi=rowMeans(chrzg, na.rm=T)

#########           BGC
#need to change EV1 to scale from 0 and 1, so that 0=HEWA, 1=TOWAe
ev1=-ddfull$EV1
ddfull$EV1.HI01=(ev1-min(ev1, na.rm=T))/(max(ev1, na.rm=T)-min(ev1, na.rm=T))
plot(ddfull$EV1.HI01, ddfull$chr5.hi)

p=ddfull$chr5.hi; 
h=ddfull$EV1.HI01
m=nls(p~h+2*(h-h^2)*(a+b*(2*h-1)), start=list(a=0, b=0), trace=T, algorithm="port")
	summary(m)
logLik(m)
#Null models to compare AIC with
m.null.b=nls(p~h+2*(h-h^2)*(a+0*(2*h-1)), start=list(a=0), trace=T, algorithm="port")
summary(m.null.b); logLik(m.null.b); AIC(m.null.b)
m.null.a=nls(p~h+2*(h-h^2)*(0+b*(2*h-1)), start=list(b=0), trace=T, algorithm="port")
summary(m.null.a); logLik(m.null.a); AIC(m.null.a)
plot(h, p)
x<-data.frame(h=seq(min(h, na.rm=T), max(h, na.rm=T), length.out=1000))
 y<-predict(m,x)
lines(x[,1], y, col="cyan") #b=0.7088; 95CI: 0.1312014 1.2864631 
library(nls2)  
confint(m.null.a)
col.to.rgb=function(col, fr)
	{colm=col2rgb(col);output={}; for(i in colm){output=c(output, i/255)}; 
		return(rgb(output[1],output[2],output[3], fr))}
	
#------------------------------------------------------------------------
#MCMC confidence interval
source('predictNLS_function.R', chdir = TRUE)
newdat=data.frame(h=seq(min(ddfull$EV1.HI01, na.rm=T), max(ddfull$EV1.HI01, na.rm=T), length.out=1000))
#newdat=data.frame(h=seq(min(matching$EV1.01, na.rm=T), max(matching$EV1.01, na.rm=T), length.out=100))
predicted.intervals<-predictNLS(m, newdata=newdat)
newdat$lci<-predicted.intervals[,6]
newdat$uci<-predicted.intervals[,7]
newdat$fit<-predicted.intervals[,1]
#png(width=900, height=900, res=200, family="Arial", filename="cheekdarkeningVSgenomicHI.png",pointsize=8)
ylab="Chr5 gene block ancestry"
p=ddfull$chr5.hi; h=ddfull$EV1.HI01
m=nls(p~h+2*(h-h^2)*(a+b*(2*h-1)), start=list(a=0, b=0), trace=T, algorithm="port")
	summary(m)
plot(p~h,  pch=24, col=col.to.rgb("darkorange",0.6), lwd=2.5, ylab=ylab, xlab="Genomic Ancestry", cex.lab=1.3, cex=1.3, xlim=c(0,1), ylim=c(0,1), bg=col.to.rgb("darkorange",0.3))
abline(0, 1, col="black", lwd=2.5,lty=2)
#add polygons for 95% confidence intervals
i.for <- order(newdat$h)
i.back<- order( newdat$h , decreasing = TRUE )
x.polygon <- c( newdat$h[i.for] , newdat$h[i.back] )
#add models in here
y.polygon<- c( newdat$uci[i.for] , newdat$lci[i.back])
polygon( x.polygon , y.polygon, col =col.to.rgb("darkorange", 0.4), border = NA )#(0/255, 191/255, 255/255, 0.4) 
#points(p~h,  col="royalblue3",  pch=21, cex=1.3, lwd=2.5)
with(newdat, lines(x=h, y=fit, col="darkorange", lty=1, lwd=3)) #royalblue3
#dev.off()
#------------------------------------------------------------------------
p=ddfull$mt.hap
h=ddfull$EV1.HI01
m=nls(p~h+2*(h-h^2)*(a+b*(2*h-1)), start=list(a=0, b=0), trace=T, algorithm="port")
	summary(m)
logLik(m)
ylab="mtDNA ancestry"
plot(p~h,  pch=24, col=col.to.rgb("deepskyblue",0.2), lwd=2.5, ylab=ylab, xlab="Genomic Ancestry", cex.lab=1.3, cex=1.3, xlim=c(0,1), ylim=c(0,1), bg=col.to.rgb("deepskyblue",0.2))
abline(0, 1, col="black", lwd=2.5,lty=2)
#add polygons for 95% confidence intervals
i.for <- order(newdat$h)
i.back<- order( newdat$h , decreasing = TRUE )
x.polygon <- c( newdat$h[i.for] , newdat$h[i.back] )
#add models in here
y.polygon<- c( newdat$uci[i.for] , newdat$lci[i.back])
polygon( x.polygon , y.polygon, col =col.to.rgb("deepskyblue", 0.4), border = NA )#(0/255, 191/255, 255/255, 0.4) 
#points(p~h,  col="royalblue3",  pch=21, cex=1.3, lwd=2.5)
with(newdat, lines(x=h, y=fit, col="deepskyblue", lty=1, lwd=3)) #royalblue3
#dev.off()
#------------------------------------------------------------------------

#MCMC confidence interval for BGC
newdat=data.frame(h=seq(min(h, na.rm=T), max(, na.rm=T), length.out=100))
#newdat=data.frame(h=seq(min(matching$EV1.01, na.rm=T), max(matching$EV1.01, na.rm=T), length.out=100))
predicted.intervals<-predictNLS(m, newdata=newdat)
newdat$lci<-predicted.intervals[,6]
newdat$uci<-predicted.intervals[,7]
newdat$fit<-predicted.intervals[,1]
#png(width=900, height=900, res=200, family="Arial", filename="cheekdarkeningVSgenomicHI.png",pointsize=8)
plot(p~h,  pch=24, col="lightblue", lwd=2.5, ylab="co-functioning nuclear HI", xlab="Genomic HI", cex.lab=1.3, cex=1.3, xlim=c(0,1), ylim=c(0,1), bg=rgb(50/100, 100/100, 83/100, 0.2))#ylab="mt HI", xlab="co-functioning nuclear HI"
#add polygons for 95% confidence intervals
i.for <- order(newdat$h)
i.back<- order( newdat$h , decreasing = TRUE )
x.polygon <- c( newdat$h[i.for] , newdat$h[i.back] )
#add models in here
y.polygon<- c( newdat$uci[i.for] , newdat$lci[i.back])
polygon( x.polygon , y.polygon, col = rgb(127/255, 255/255, 212/255, 0.4), border = NA )#(0/255, 191/255, 255/255, 0.4) 
#points(p~h,  col="royalblue3",  pch=21, cex=1.3, lwd=2.5)
abline(0, 1, col="red", lwd=2.5,lty=2)
with(newdat, lines(x=h, y=fit, col="lightseagreen", lty=1, lwd=2.5)) #royalblue3
#dev.off()

#make PREDICTIOn concept plot
plot(ddfull$EV1.HI01, ddfull$EV1.HI01,  pch=24, col="lightblue", lwd=2.5, ylab="Local Ancestry", xlab="Genomic HI", cex.lab=1.3, cex=1.3, xlim=c(0,1), ylim=c(0,1), bg=rgb(50/100, 100/100, 83/100, 0.2))#ylab="mt HI", xlab="co-functioning nuclear HI"
abline(0, 1, col="red", lwd=2.5,lty=2)
bgc=function(a, b, h){p=h+2*(h-h^2)*(a+b*(2*h-1))}
h=seq(0,1, length.out=100)
lines(h,bgc(a=0,b=0.8,h), col="blue", lty=2, lwd=3 ); text(0.2,0.8, expression(paste(alpha, "= 0, ", beta, " = 0.8")))
lines(h,bgc(a=-0.5,b=0,h), col="blue", lty=2, lwd=3 ); text(0.2,0.8, expression(paste(alpha, "= -0.5, ", beta, " = 0")))


