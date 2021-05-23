library(SNPRelate)
library(gdsfmt)
library("GWASTools")
g<-snpgdsOpen("silu.maddie.197.gds",  readonly = TRUE)
snpgdsSummary(g)
set.seed(1000)
# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(g, ld.threshold=0.2)
names(snpset)
head(snpset$chr1)

#PCA, 
pca <- snpgdsPCA(g)
pca$varprop*100
library(stringr)
id<-str_split_fixed(pca$sample.id, "_", 3)[,3]
tab <- data.frame(sample.id = id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)
#write.csv(tab, "~/Desktop/GBS_TOWAHEWA/Sep_17_analysis/tohe_plate1_pca.csv")
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
text(tab$EV2, tab$EV1, tab$sample.id)
 tab$sample.id=as.character( tab$sample.id)
s.name=str_split_fixed(tab$sample.id, "B", 2)[,1:2]
tab$id=c(s.name[1:81,2], s.name[82:196,1])
head(d); tabd={}; d$short.name=as.character(d$short.name); 
erow=rep(NA, ncol(d));names(erow)=colnames(d)
for(i in tab$sample.id)
{row=which(d$short.name==i)
	if(length(row)>0)
		{tabd=rbind(tabd, d[row,])}
	else{tabd=rbind(tabd, erow)}}

tabd$pop.SP=tabd$SP
tabd$pop.SP[which(tabd$siteID=="Queen_Charlotte Is._BC")]="hadiagwaii"
tabd$pop.SP[which(tabd$siteID=="Haida_Gwaii")]="hadiagwaii"
tabd$pop.SP[which(tabd$siteID=="Valdez_AK")]="valdez"
tabd=data.frame(tab$id, tabd)
#missing ones find match in another datasheet
names=as.character(tab[which(is.na(tabd$pop.SP)==1),1]); missingd={}
f=read.csv("TOHE_silu1234.maddie12.site.info.csv")
na.row=rep(NA, ncol(f)); names(na.row)=colnames(f)
f$short.name=as.character(f$short.name); 
for(i in names)
{row=which(i==f$short.name)
if(length(row)>0)
	{missingd=rbind(missingd, f[row,])}
else{print(i); missingd=rbind(missingd,na.row); mt.miss=c(mt.miss, i)}
}

#fill it in tabd for N, W, pop.SP, site.name
tabd$pop.SP=as.character(tabd$pop.SP); tabd$siteID=as.character(tabd$siteID); missingd$pop.SP=as.character(missingd$pop.SP); missingd$site.name=as.character(missingd$site.name)
for(i in which(is.na(tabd$sample.id)==1))
	{row=which(as.character(missingd$short.name)==tab$sample.id[i])
		if(length(row)>0)
			{tabd$N[i]=missingd$N[row];tabd$W[i]=missingd$W[row]; 			tabd$pop.SP[i]=missingd$pop.SP[row];tabd$siteID[i]=missingd$site.name[row];}
		else{print(tab$sample.id[i])}}
tabd[which(is.na(tabd$sample.id)==1),]
table(tabd$pop.SP); table(tabd$siteID)
tabd$pop.ids=as.character(tabd$pop.SP)
tabd$pop.ids[which(tabd$pop.SP=="inlandHEWA")]="HEWA"
tabd$pop.ids[which(tabd$pop.SP=="coastalHEWA")]="HEWA"
tabd$pop.ids[which(tabd$pop.SP=="coastalTOWA")]=tabd$siteID[which(tabd$pop.SP=="coastalTOWA")]
table(tabd$pop.ids)
tabd$pop.ids=factor(tabd$pop.ids, levels=c("HEWA", "inlandTOWA","valdez", "hadiagwaii", "Bella_Coola", "Haines_AK", "N._Vancouver_Island_BC","Prince_Rupert"))
fc <- c("turquoise", "magenta","royalblue","seagreen","coral1", "gold","limegreen","khaki4")
Col <- fc[as.numeric(as.factor(tabd$pop.ids))]

#add missing mt hap
mt=read.csv("mt.pop.226indv.towa.hewa.burke.Jan2019.csv")
names=as.character(tabd[which(is.na(tabd$mt.hap)==1),1]); 
for(i in names)
{row=which(tabd[,1]==i)
rowadd=which(mt$Burke.rec==i);
if(length(rowadd)==1) {print(i); tabd$mt.hap[row]=mt$mt.hi[rowadd]}}



plot(-tab$EV1,tab$EV2,  xlab="Genomic EV1 (1.70%)", cex=1.5,ylab="Genomic EV2 (1.10%)", bg=Col, pch=21)
# points(tab$EV2[which(tab$sample.id=="QF10W06")], tab$EV1[which(tab$sample.id=="QF10W06")],col="red", pch=16)
# points(tab$EV2[which(tab$sample.id=="QE10W06")], tab$EV1[which(tab$sample.id=="QE10W06")],col="red", pch=16)
legend("topleft", title="Coastal TOWA", legend=c("HEWA", "inlandTOWA","valdez", "hadiagwaii", "Bella_Coola", "Haines_AK", "N._Vancouver_Island_BC","Prince_Rupert"), pch=21, pt.bg=fc, cex=0.6, pt.cex=1.5)

inhewa.tab=tabd[which(tabd$pop.SP=="inlandHEWA"),]
cohewa.tab=tabd[which(tabd$pop.SP=="coastalHEWA"),]

inhewa.tab[which((-inhewa.tab$EV1)>-0.01),]
cohewa.tab[which((-cohewa.tab$EV1)>-0.01),]
cohewa.tab[which((cohewa.tab$EV2)>0.05),]
inhewa.tab[which((inhewa.tab$EV2)>0.05),]

saved=data.frame(pca$sample.id, tab, tabd)
#write.csv(saved, "silu.maddie.197.genomicEV.mt.csv", row.names=F)

############
################ SAMPLE SIZE
saved=read.csv("silu.maddie.197.genomicEV.mt.csv")
table(saved$pop.SP ); table(saved$pop.ids)
hewa.ids=saved$pca.sample.id[which(saved$pop.id=="HEWA")]
hewaid<-str_split_fixed(hewa.ids, "_", 3)[,3]
sort(hewaid)
itowa.ids=saved$pca.sample.id[which(saved$pop.id=="inlandTOWA")]; 
itowaid=str_split_fixed(itowa.ids, "_", 3)[,3]; sort(itowaid)
itowad=saved[which(saved$pop.id=="inlandTOWA"),]
itowad$recnumber=itowaid
itowad=itowad[order(itowaid),]
tail(itowad,8)

table(str_split_fixed(saved$pca.sample.id, "_", 3)[,1])

