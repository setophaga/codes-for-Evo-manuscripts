library(googleway)
library(vegan)
library(ggplot2);library(ggbiplot)
library("ggmap")
cc=read.csv("mitonuclear.stow.socc.climate.csv")

#########associating climate with Genotypes
#make sure the site names are right
#first, run climate PC and compare across socc, inland and coastal stow habitat
climd=cc[,which(colnames(cc)=="AHM"):which(colnames(cc)=="Tmin01")]
clim.pca=prcomp(climd, center=T, scale=T);biplot(clim.pca)
summary(clim.pca); (clim.pca$rotation)
ggbiplot(clim.pca, obs.scale = 1, var.scale = 1,groups = cc$sp.pop, ellipse = TRUE, circle = F)+theme(legend.direction = 'horizontal',legend.position = 'top')+ guides(fill=FALSE)+theme_classic()+ scale_color_manual(values=c("darkolivegreen3", "magenta", "turquoise"))
cc$clim.pc1=clim.pca$x[,1];cc$clim.pc2=clim.pca$x[,2];
library(vioplot)
pop=factor(cc$sp.pop, levels=c("SOCC", "coastal.STOW", "inland.STOW"))
pop.col=c( "turquoise","darkolivegreen3", "magenta")
vioplot(cc$clim.pc1~pop,col=rgb.convert(pop.col, 0.3), xaxt="n",border = pop.col,horizontal = F, las = 1, cex=1.2, ylab="Climate PC1", xlab="", lwd=2)
stripchart(cc$clim.pc1~pop,vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1.5, col=rgb.convert(pop.col, 1), jitter=0.3)  
vioplot(clim.pc1~pop,col=rgb.convert(pop.col, 0.1), xaxt="n",border = pop.col,horizontal = F, las = 1, cex=1.2, ylab="Climate PC1", xlab="", lwd=2, add=T)
kruskal.test(mt.nu.comb$clim.pc1~mt.nu.comb$sp.pop)
pairwise.wilcox.test(mt.nu.comb$clim.pc1,mt.nu.comb$sp.pop,p.adjust.method = "BH")

vioplot(cc$clim.pc2~pop,col=rgb.convert(pop.col, 0.3), xaxt="n",border = pop.col,horizontal = F, las = 1, cex=1.2, ylab="Climate PC2", xlab="", lwd=2)
stripchart(cc$clim.pc2~pop,vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1.5, col=rgb.convert(pop.col, 1), jitter=0.3)  
vioplot(cc$clim.pc2~pop,col=rgb.convert(pop.col, 0.1), xaxt="n",border = pop.col,horizontal = F, las = 1, cex=1.2, ylab="Climate PC2", xlab="", lwd=2, add=T)
kruskal.test(cc$clim.pc2~cc$sp.pop)
pairwise.wilcox.test(cc$clim.pc2,cc$sp.pop,p.adjust.method = "BH")


#first need to find mtDNA and nDNA ancestry for each sie
mt$mt.hi=as.character(mt$mt.hap)
mt$mt.hi[which(mt$mt.hi=="h")]=0; mt$mt.hi[which(mt$mt.hi=="t")]=1; 
mt$mt.hi=as.numeric(as.character(mt$mt.hi))
mt$site=as.character(mt$site)
#John Day_OR for mt into John_Day_OR ;  Queen_Charlotte_Is._BC into  Haida_Gwaii 
mt$site[which(mt$site=="John Day_OR")]="John_Day_OR"
mt$site[which(mt$site=="Queen_Charlotte_Is._BC")]="Haida_Gwaii"
mt.dd=aggregate(mt, list(mt$site),  mean, na.rm=T); mt.dd$W=-mt.dd$W
mtd=data.frame(mt.dd$Group.1, mt.dd$mt.hi, mt.dd$N, mt.dd$W)
colnames(mtd)=c("siteID", "mt.hi", "N", "W")

#find which mtd site is not in cc 
for(i in 1:length(mtd$siteID)){if(length(which(mtd$siteID[i]==cc$site ))==0){print(mtd[i,])}}
############################
#Figure 1.   A mtDNA map 
#int mt map
map=get_map(location =c(-128, 52), maptype="hybrid" , zoom = 4)
ggmap(map) + geom_point(aes(x = W, y = N, col=mt.hi), data =mtd, alpha = 0.7, size = 4)+labs(x="Longitude", y="Latitude",colour="mt ancestry")+scale_color_gradient2(low="turquoise",midpoint=0.5, mid="gold", high="magenta") #, 
#Figure 6 A ch5 peak 1 map
ggmap(map) + geom_point(aes(x = W, y = N, col=chr5.hi), data =nu, alpha = .7, size = 4)+labs(x="Longitude", y="Latitude",colour="chr 5 SNPs")+scale_color_gradient2(low="turquoise",midpoint=0.5, mid="gold", high="magenta") 
#EV1
ggmap(map) + geom_point(aes(x = W, y = N, col=EV1.HI01), data =nu, alpha = .7, size = 4)+labs(x="Longitude", y="Latitude",colour="Genomic EV1")+scale_color_gradient2(low="turquoise",midpoint=0.5, mid="gold", high="magenta") 

#match mt dataset and nu dataset
mtd$siteID=as.character(mtd$siteID)
colnames(nu)[which(colnames(nu)=="Group.1")]="site"
nu$site=as.character(nu$site)
nu$site[which(nu$site=="Bella_Coola")]="Bella_Coola_BC"
nu$site[which(nu$site=="Dease_Lake")]="Dease_Lake_BC"
nu$site[which(nu$site=="Prince_Rupert")]="Prince_Rupert_BC"
nu$site[which(nu$site=="Valemount_McBride")]="Valemount_McBride_BC"
sites=intersect(unique(mtd$siteID), unique(nu$site))
for(s in unique(mtd$siteID)){if(length(which(s==sites))==0){}}
mtnud={}
for(s in sites)
{mt.row=which(mtd$siteID==s)
nu.row=which(nu$site==s)
rowtoadd=c(mtd[mt.row,], nu[nu.row,])
mtnud=rbind(mtnud, rowtoadd)}
mtnud=data.frame(mtnud)
mtnud$siteID=sites
plot(mtnud$W, mtnud$N, col="green", pch=16)
text(mtnud$W, mtnud$N, mtnud$siteID)

#NOW match climate data frome cc with the genotype files
climpc1={}; climpc2={}; n={};w={}; sp={}; 
for(s in mtnud$siteID)
{row=which(cc$site==s)
if(length(row)==1)
{climpc1=c(climpc1, cc$clim.pc1[row]);climpc2=c(climpc2, cc$clim.pc2[row])
w=c(w, cc$W[row]);n=c(n, cc$N[row]); sp=c(sp, as.character(cc$sp.pop[row]))}
else{print(s)
climpc1=c(climpc1, NA); climpc2=c(climpc2, NA);w=c(w, NA);n=c(n, NA); sp=c(sp, NA); }
}
mt.nu.comb=data.frame(mtnud, w,n, sp, climpc1, climpc2)

for(j in 1:ncol(mt.nu.comb)){mt.nu.comb[,j]=as.character(mt.nu.comb[,j])}
mt.nu.comb$chr5.mt=(as.numeric(mt.nu.comb$mt.hi)+as.numeric(mt.nu.comb$chr5.hi))/2


#plot climate data on MAP
library(ggmap)
map=get_map(location =c(-128, 52), maptype="hybrid" , zoom = 4) #terrain-background
plotd=data.frame(mt.nu.comb$W, mt.nu.comb$N, mt.nu.comb$climpc1, mt.nu.comb$climpc2, mt.nu.comb$chr5.mt); colnames(plotd)=c("W", "N", "clim.pc1", "clim.pc2", "chr5.mt")
plotd$N=as.numeric(as.character(plotd$N));plotd$W=as.numeric(as.character(plotd$W))
plotd$clim.pc2=as.numeric(as.character(plotd$clim.pc2));plotd$clim.pc1=as.numeric(as.character(plotd$clim.pc1)); plotd$chr5.mt=as.numeric(as.character(plotd$chr5.mt))

ggmap(map) + geom_point(aes(x = W, y = N, col=clim.pc1), data =plotd, alpha = .9,size=4)+scale_color_gradient2(low="yellow",midpoint=0, mid="green", high="dodgerblue")+labs(x="Longitude", y="Latitude", colour="Climate PC1")
ggmap(map) + geom_point(aes(x = W, y = N, col=clim.pc2), data =plotd, alpha = .9,size=4)+scale_color_gradient2(low="dodgerblue",midpoint=-0.46, mid="green", high="yellow")+labs(x="Longitude", y="Latitude", colour="Climate PC2")
ggmap(map) + geom_point(aes(x = W, y = N, col=chr5.mt), data =plotd, alpha = .9,size=4)+scale_color_gradient2(low="yellow",midpoint=0.5, mid="green", high="dodgerblue")+labs(x="Longitude", y="Latitude", colour="Mitonuclear HI")


#### MANTEL TEST
station.dists <- dist(cbind(mt.nu.comb$w,mt.nu.comb$n)) #y is w; x is n
climpc1.dists=dist(mt.nu.comb$clim.pc1) 
climpc2.dists=dist(mt.nu.comb$clim.pc2)

#######test correlation of mt-nu combo ~climate matrix
mantel_partial.test=function(clim, geno)
{d=data.frame(clim, geno, mt.nu.comb$n, mt.nu.comb$w); d=na.omit(d)
clim.dist=dist(d$clim); geno.dist=dist(d$geno); station.dist=dist(cbind(d$mt.nu.comb.w, d$mt.nu.comb.n))
mantel.partial( geno.dist,clim.dist, station.dist,permutations=10000)}
#test chrz marker and mtDNA ancestry 
mantel_partial.test(mt.nu.comb$mt.hi, mt.nu.comb$chr5.hi)
mantel_partial.test(mt.nu.comb$climpc1, mt.nu.comb$chr5.mt)
mantel_partial.test(mt.nu.comb$climpc2, mt.nu.comb$chr5.mt)

################################
#without Valdez_AK; Haida_Gwaii
combo.d=mt.nu.comb[-c(which(mt.nu.comb$siteID=="Haida_Gwaii"), which(mt.nu.comb$siteID=="Valdez_AK")),]
station.dists <- dist(cbind(combo.d$w,combo.d$n))
clim.dist=dist(combo.d$clim.pc1); geno.dist=dist(combo.d$chr5.mt); 
chr5.dist=dist(combo.d$chr5.hi); mt.dist=dist(combo.d$mt.hi)
mantel.partial( geno.dist,clim.dist, station.dists,permutations=10000)
mantel.partial( chr5.dist,mt.dist, station.dists,permutations=10000)
################################
################################
################################
mantel_partial.test(mt.nu.comb$clim.pc1, mt.nu.comb$chr5.mt)


###.             Figure 7 B
#####CLIMATE PC1 & CHR5.MT
col2rgb("magenta")/255; col2rgb("turquoise")/255;col2rgb("khaki4")/255;
colors=c(rgb(0.2509804,0.8784314,0.8156863,1), rgb(0.5450980,0.5254902,0.3058824,1),rgb(0.2509804,0.8784314,0.8156863,1),rgb(1,0,1,1))
plot(mt.nu.comb$clim.pc1, mt.nu.comb$chr5.mt,  cex=1.5, xlab="Climate PC1", ylab="mt-nDNA ancestry",col=colors[as.numeric(as.factor(mt.nu.comb$sp))], pch=16)
points(mt.nu.comb$clim.pc1[which(mt.nu.comb$siteID=="Bella_Coola_BC")], mt.nu.comb$chr5.mt[which(mt.nu.comb$siteID=="Bella_Coola_BC")],bg="coral", pch=21, cex=1.5)
points(mt.nu.comb$clim.pc1[which(mt.nu.comb$siteID=="Haida_Gwaii")], mt.nu.comb$chr5.mt[which(mt.nu.comb$siteID=="Haida_Gwaii")],bg="seagreen", pch=21, cex=1.5)
points(mt.nu.comb$clim.pc1[which(mt.nu.comb$siteID=="Haines_AK")], mt.nu.comb$chr5.mt[which(mt.nu.comb$siteID=="Haines_AK")],bg="gold", pch=21, cex=1.5)
points(mt.nu.comb$clim.pc1[which(mt.nu.comb$siteID=="N._Vancouver_Island_BC")], mt.nu.comb$chr5.mt[which(mt.nu.comb$siteID=="N._Vancouver_Island_BC")],bg="limegreen", pch=21, cex=1.5)
points(mt.nu.comb$clim.pc1[which(mt.nu.comb$siteID=="Prince_Rupert_BC")],mt.nu.comb$chr5.mt[which(mt.nu.comb$siteID=="Prince_Rupert_BC")],bg="Khaki4", pch=21, cex=1.5)
points(mt.nu.comb$clim.pc1[which(mt.nu.comb$siteID=="Valdez_AK")], mt.nu.comb$chr5.mt[which(mt.nu.comb$siteID=="Valdez_AK")],bg="royalblue",  cex=1.5, pch=21)
legend(0.5, 0.87, title="Coastal STOW", pch=21, pt.bg=c("coral1","seagreen", "gold","limegreen","khaki4","royalblue"), c("Bella_Coola", "Haida_Gwaii", "Haines","N._Vancouver_Island","Prince_Rupert", "Valdez"), cex=0.6, pt.cex=1)
#fitting 


