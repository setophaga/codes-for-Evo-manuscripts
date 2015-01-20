bootstrap_corr<-function(a, b, n)
{r.boot<-{}
for (i in 1:n)
	{a.s<-sample(a, length(a), replace=T)
	b.s<-sample(b, length(b), replace=T)
	if(sd(a.s)==0||sd(b.s)==0)
		{r1<-0}
	else {r1<-cov(a.s,b.s)/(sd(a.s)*sd(b.s))}
	r.boot<-c(r.boot,r1)}
r<-cov(a,b)/(sd(a)*sd(b))
m<-0
for (i in 1: n)
	{if(r.boot[[i]]>r||r.boot[[i]]==r)
		{m<-m+1}}
p<-m/n
se.rboot<-sd(r.boot)/sqrt(n)
return(c(r, se.rboot, p))
}
