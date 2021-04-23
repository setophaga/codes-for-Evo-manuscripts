library(ggplot2)
setwd("/Users/siluwang/Desktop/1_PhD.projects/1.0.CoastalTOWA/1.Oct.coastal.towa.2019.MS.version/Ref.yrwa.Dec.2020/pi.diversity")
cstow.p=read.csv("pixy_pi.socc.stow.cstow.sub.filter.txt",header=T, sep="\t")
ocstow=cstow.p[which(cstow.p$pop=="ocstow"),]
haida.chr5=cstow.p[which(cstow.p$pop=="haidagwaii"),]
val.chr5=cstow.p[which(cstow.p$pop=="valdez"),]

pixy=read.csv("pixy_pi.socc.stow.cstow.filter.txt",header=T, sep="\t")
istow.chr5=pixy[which(pixy$pop=="istow"),]
cstow=pixy[which(pixy$pop=="cstow"),]
socc=pixy[which(pixy$pop=="socc"),]

col.to.rgb=function(col, fr)
	{colm=col2rgb(col);output={}; for(i in colm){output=c(output, i/255)}; 
		return(rgb(output[1],output[2],output[3], fr))}
		
plot(socc$window_pos_1 , socc$avg_pi, pch=21, bg=col.to.rgb("turquoise", 0.6),xlim=c(2.5e6, 7e6),  ylab=expression(pi),xlab="Positions" ,  ylim=c(0,0.02))#[1:5000]
lines(socc$window_pos_1 , socc$avg_pi, col=col.to.rgb("turquoise", 1), lwd=2, lty=2)
#inland stow
plot(istow.chr5$window_pos_1 , istow.chr5$avg_pi, pch=21, bg=col.to.rgb("magenta", 0.1),xlim=c(2.5e6, 7e6),  ylab=expression(pi),xlab="Positions" , ylim=c(0,0.02))
	lines(istow.chr5$window_pos_1 , istow.chr5$avg_pi, col=col.to.rgb("magenta", 1), lwd=2, lty=2)
abline( v=c(3726725, 5172249),col="salmon",lwd=6,  lty=3)
abline(v=cc(2.5e6, 7e6))
#other coastal stow
plot(ocstow$window_pos_1 , ocstow$avg_pi ,bg=col.to.rgb("gold", 0.5), pch=21,xlim=c(2.5e6, 7e6),  ylab=expression(pi),xlab="Positions" , ylim=c(0,0.02))
lines(ocstow$window_pos_1 , ocstow$avg_pi, col=col.to.rgb("gold", 1), lwd=2, lty=2)
abline( v=c(3726725, 5172249),col="salmon",lwd=6,  lty=3)

#haida gwaii
plot(haida.chr5$window_pos_1 , haida.chr5$avg_pi ,bg=col.to.rgb("seagreen", 0.5), pch=21,xlim=c(2.5e6, 7e6),  ylab=expression(pi),xlab="Positions" , ylim=c(0,0.02))
lines(haida.chr5$window_pos_1 , haida.chr5$avg_pi, col=col.to.rgb("seagreen", 1), lwd=2, lty=2)
abline( v=c(3726725, 5172249),col="salmon",lwd=6,  lty=3)

#other valdez
plot(val.chr5$window_pos_1 , val.chr5$avg_pi ,bg=col.to.rgb("royalblue", 0.5), pch=21,xlim=c(2.5e6, 7e6),  ylab=expression(pi),xlab="Positions" , ylim=c(0,0.02))
lines(val.chr5$window_pos_1 , val.chr5$avg_pi, col=col.to.rgb("royalblue", 1), lwd=2, lty=2)
abline( v=c(3726725, 5172249),col="salmon",lwd=6,  lty=3)

#all coastal stow
plot(cstow$window_pos_1 , cstow$avg_pi ,bg=col.to.rgb("chartreuse2", 0.5), pch=21,xlim=c(2.5e6, 7e6),  ylab=expression(pi),xlab="Positions" , ylim=c(0,0.02))
lines(cstow$window_pos_1 , cstow$avg_pi, col=col.to.rgb("chartreuse", 1), lwd=1, lty=2)
abline(v=c(3726725,5172249), col="salmon", lty=3, lwd=6)
abline(v=2.5e6, lwd=0.4)

#all coastal stow, all chr5
plot(cstow$window_pos_1 , cstow$avg_pi ,bg=col.to.rgb("chartreuse2", 0.5), pch=21,ylab=expression(pi),xlab="Positions" )
lines(cstow$window_pos_1 , cstow$avg_pi, col=col.to.rgb("chartreuse", 1), lwd=1, lty=2)
abline(v=c(3726725,5172249), col="salmon", lty=3)
abline(v=c(2.5e6, 7e6), lwd=0.4)
#zoom in around 2.5e6, 7e6
plot(cstow$window_pos_1 , cstow$avg_pi ,bg=col.to.rgb("chartreuse2", 0.5), pch=21,ylab=expression(pi),xlab="Positions" , xlim=c(2.5e6, 7e6), ylim=c(0, 0.03))
lines(cstow$window_pos_1 , cstow$avg_pi, col=col.to.rgb("chartreuse", 1), lwd=1, lty=2)
abline(v=c(3726725,5172249), col="salmon", lty=3, lwd=6)
abline(v=c(2.5e6, 7e6), lwd=0.4)

boot.cand.rest=function(d, iterations)
{cand=intersect(which(d$window_pos_1  < 5172249), which(d$window_pos_1  >3726725)) ; n.cand=length(cand); 
rest=c(which(d$window_pos_1  > 5172249), which(d$window_pos_1  < 3726725));n.rest=length(rest)
cand.boot={}; rest.boot={}
for(i in 1:iterations)
{cand.boot=c(cand.boot,mean(d$avg_pi[sample(cand, n.cand, replace=T)], na.rm=T))
rest.boot=c(rest.boot, mean(d$avg_pi[sample(rest, n.rest, replace=T)], na.rm=T))}
return(data.frame(cand.boot, rest.boot))}
#other coastal STOW
oc.boot=boot.cand.rest(ocstow, 10000)
hist(oc.boot[,1], xlim=c(0, 0.005), col=col.to.rgb("salmon",0.5), breaks=200,border= col.to.rgb("salmon",1))
hist(oc.boot[,2], add=T, col=col.to.rgb("gold",0.5), breaks=200,border= "gold")
round(quantile(oc.boot[,1],c(0.025, 0.975)),4); round(quantile(oc.boot[,2],c(0.025, 0.975)),4)

#haida gwaii
haida.boot=boot.cand.rest(haida.chr5, 10000)
hist(haida.boot[,2], xlim=c(0, 0.005),  col=col.to.rgb("seagreen",0.6), breaks=200,border= "seagreen")
hist(haida.boot[,1], add=T,col=col.to.rgb("salmon",0.5), breaks=200,border= col.to.rgb("salmon",1))
round(quantile(haida.boot[,1],c(0.025, 0.975)),4); round(quantile(haida.boot[,2],c(0.025, 0.975)),4)

#valdez
val.boot=boot.cand.rest(val.chr5, 10000)
hist(val.boot[,1], xlim=c(0, 0.005), col=col.to.rgb("salmon",0.5), breaks=200,border= col.to.rgb("salmon",1))
hist(val.boot[,2], add=T, col=col.to.rgb("royalblue",0.5), breaks=200,border= "royalblue")
round(quantile(val.boot[,1],c(0.025, 0.975)),4); round(quantile(val.boot[,2],c(0.025, 0.975)),4)

#all coastal STOW
c.boot=boot.cand.rest(cstow, 10000)
hist(c.boot[,1], xlim=c(0.001, 0.005), col=col.to.rgb("darkorange",0.5), breaks=200,border= col.to.rgb("darkorange",1), main="", xlab=expression(pi))
hist(c.boot[,2], add=T, col=col.to.rgb("chartreuse2",0.5), breaks=200,border= "chartreuse2")
round(quantile(c.boot[,1],c(0.025, 0.975)),4); round(quantile(c.boot[,2],c(0.025, 0.975)),4)
