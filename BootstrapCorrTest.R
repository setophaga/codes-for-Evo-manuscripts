#this function testst correlation of two variables: a and b, with n bootstrap samples
bootstrap_corr<-function(a, b, n)
{r.boot<-{}
for (i in 1:n)
	{a.s<-sample(a, length(a), replace=T) #sample a with replacement, so that this sample vector is the same length as a
	b.s<-sample(b, length(b), replace=T)  #sample b with replacement, so that this sample vector is the same length as b
	if(sd(a.s)==0||sd(b.s)==0)            #if the standard deviation of a or b is zero, 
		{r1<-0}                           #there is no correlation, correlation coefficient is zero
	else {r1<-cov(a.s,b.s)/(sd(a.s)*sd(b.s))} #correlation coefficient is a function of covrianace and standard deviation of a, b
	r.boot<-c(r.boot,r1)}        #record correlation coefficient calculated from each bootstrap sample, to form r.boot of n values 
r<-cov(a,b)/(sd(a)*sd(b))        #calcualte correlation coefficient of the data as "r"
m<-0  
for (i in 1: n)                 
	{if(r.boot[[i]]>r||r.boot[[i]]==r)
		{m<-m+1}}                 #this loop go thorugh r.boot and count the cases when a value is equal to or greater than r 
p<-m/n                            #p-value is the proportion of such cases out of total number of cases
se.rboot<-sd(r.boot)/sqrt(n)      #calculate bootstrap sample error of correlation coefficient
return(c(r, se.rboot, p))         #return 1) correlation coeffient of the data, 2)bootstrap sample error, 3)p-value
}
