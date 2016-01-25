#function for probability matrix
prob.matrix<-function(m)
	{for(i in 1:length(m[,1]))
		{rowsum<-sum(m[i,])
			for(j in 1:length(m[1,]))
				{m[i,j]<-m[i,j]/rowsum}
		}
		n<-length(m[,1])
		return(m[-n, -n])
	}

#function to get a transition vector with the probability matrix
steps.to.hhi<-function(pmatrix)
{pm<-as.matrix(pmatrix)
diag(pm)<-0
t<-solve(diag(length(pm[,1]))-pm)
one<-matrix(rep(1, length(pm[,1])), nrow=length(pm[,1]), ncol=1)
t<-t%*%one
return(t)
}

#a function combining the above two functions
transition.calc<-function(m)
	{x<-prob.matrix(m)
	y<-steps.to.hhi(x)
	return(y)
	}
