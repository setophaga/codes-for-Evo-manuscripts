
library(lattice)
library(RColorBrewer)
setwd("~/Desktop/1_PhD.projects/1.0.CoastalTOWA/1.Oct.coastal.towa.2019.MS.version/Ref.yrwa.Dec.2020/pi.ndiversity.ld/")
hai=read.csv("~/Desktop/1_PhD.projects/1.0.CoastalTOWA/1.Oct.coastal.towa.2019.MS.version/Ref.yrwa.Dec.2020/pi.ndiversity.ld//haidagwaii.ld.geno.ld", sep="\t")

val=read.csv("~/Desktop/1_PhD.projects/1.0.CoastalTOWA/1.Oct.coastal.towa.2019.MS.version/Ref.yrwa.Dec.2020/pi.ndiversity.ld//valdez.ch5.r2.geno.ld", sep="\t")

oc=read.csv("~/Desktop/1_PhD.projects/1.0.CoastalTOWA/1.Oct.coastal.towa.2019.MS.version/Ref.yrwa.Dec.2020/pi.ndiversity.ld//oc.stow.ld.geno.ld", sep="\t")
all=read.csv("~/Desktop/1_PhD.projects/1.0.CoastalTOWA/1.Oct.coastal.towa.2019.MS.version/Ref.yrwa.Dec.2020/pi.ndiversity.ld/c.stow.ld.geno.ld", sep="\t")


coul <- colorRampPalette(c('dimgrey','deepskyblue','darkolivegreen1'))(25)
levelplot(R.2 ~ POS1*POS2, data=hai  ,xlab="X",main="", col.regions =coul, panel = function(...){panel.levelplot(...); panel.abline(h = c(3726725, 5172249)); panel.abline(v = c(3726725, 5172249))})  
 
levelplot(R.2 ~ POS1*POS2, data=val  ,xlab="X",main="", col.regions =coul, panel = function(...){panel.levelplot(...); panel.abline(h = c(3726725, 5172249)); panel.abline(v = c(3726725, 5172249))})  
  
levelplot(R.2 ~ POS1*POS2, data=oc  ,xlab="X",main="", col.regions =coul, panel = function(...){panel.levelplot(...); panel.abline(h = c(3726725, 5172249)); panel.abline(v = c(3726725, 5172249))})  
 
levelplot(R.2 ~ POS1*POS2, data=all  ,xlab="X",main="", col.regions =coul, panel = function(...){panel.levelplot(...); panel.abline(h = c(3726725, 5172249)); panel.abline(v = c(3726725, 5172249))})  
 
 col.to.rgb=function(col, fr)
	{colm=col2rgb(col);output={}; for(i in colm){output=c(output, i/255)}; 
		return(rgb(output[1],output[2],output[3], fr))}
	
boot.cand.rest=function(d, iterations)
{
cand=intersect(intersect(which(d$POS1 < 5172249), which(d$POS1 >3726725)), intersect(which(d$POS2 < 5172249), which(d$POS2 >3726725))) ; n.cand=length(cand); 
rest=seq(1:nrow(d))[-cand];n.rest=length(rest)
cand.boot={}; rest.boot={}
for(i in 1:iterations)
{cand.boot=c(cand.boot,mean(d$R.2[sample(cand, n.cand, replace=T)], na.rm=T, rm.nan=T))
rest.boot=c(rest.boot, mean(d$R.2[sample(rest, n.rest, replace=T)], na.rm=T, rm.nan=T))}
return(data.frame(cand.boot, rest.boot))}salmon

 #other coastal STOW
oc.boot=boot.cand.rest(oc, 10000)
#write.csv(oc.boot, "LD.boot.ocSTOW.csv")
oc.boot=read.csv("LD.boot.ocSTOW.csv"); oc.boot=data.frame(oc.boot$cand.boot, oc.boot$rest.boot)
hist(oc.boot[,1], xlim=c(0.04, 0.115), col=col.to.rgb("salmon",0.5), breaks=200,border= col.to.rgb("salmon",1),xlab=expression(LD (r^2)), main="")
hist(oc.boot[,2], add=T, col=col.to.rgb("gold",0.5), breaks=300,border= "gold")
round(quantile(oc.boot[,1],c(0.025, 0.975)),4); round(quantile(oc.boot[,2],c(0.025, 0.975)),4)

#haida gwaii
haida.boot=boot.cand.rest(hai, 10000)
#write.csv(haida.boot, "LD.boot.haidagwaiiSTOW.csv")
haida.boot=read.csv("LD.boot.haidagwaiiSTOW.csv"); haida.boot=data.frame(haida.boot$cand.boot, haida.boot$rest.boot)
hist(haida.boot[,2],  xlim=c(0.125, 0.19), col=col.to.rgb("seagreen",0.6), breaks=200,border= "seagreen",xlab=expression(LD (r^2)), main="")
hist(haida.boot[,1], add=T,col=col.to.rgb("salmon",0.5), breaks=400,border= col.to.rgb("salmon",1))
round(quantile(haida.boot[,1],c(0.025, 0.975)),4); round(quantile(haida.boot[,2],c(0.025, 0.975)),4)

#valdez
val.boot=boot.cand.rest(val, 10000)
#write.csv(val.boot, "LD.boot.valdezSTOW.csv")
val.boot=read.csv("LD.boot.valdezSTOW.csv"); val.boot=data.frame(val.boot$cand.boot, val.boot$rest.boot)
hist(val.boot[,1], xlim=c(0.0556, 0.22), col=col.to.rgb("salmon",0.5), breaks=200,border= col.to.rgb("salmon",1),xlab=expression(LD (r^2)), main="")
hist(val.boot[,2], add=T, col=col.to.rgb("royalblue",0.5), breaks=200,border= "royalblue")
round(quantile(val.boot[,1],c(0.025, 0.975)),4); round(quantile(val.boot[,2],c(0.025, 0.975)),4)

#all
all.boot=boot.cand.rest(all, 10000)
hist(all.boot[,1], xlim=c(0.02, 0.11), col=col.to.rgb("salmon",0.5), breaks=200,border= col.to.rgb("salmon",1),xlab=expression(LD (r^2)), main="")
hist(all.boot[,2], add=T, col=col.to.rgb("cadetblue4",0.5), breaks=200,border= "cadetblue4")
round(quantile(all.boot[,1],c(0.025, 0.975)),4); round(quantile(all.boot[,2],c(0.025, 0.975)),4)

