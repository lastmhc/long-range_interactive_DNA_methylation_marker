# R code for log-rank test of all probes (survival analysis)

# usage: Rscript R_logrank.R

# Before run, the methylation (or expression) values must be changed to index (0 or 1) 
# 0 means "over the median", 1 means "below the median". 

# censor: 0 - dead, 1 - alive

#input.txt
#ID	ID_name1	ID_name2	ID_name3	ID_name4	ID_name5	...
#censor	0	0	1	1	1	...
#survival_day	461	118	81	146	4638	...
#cg00000029	0	0	1	1	0	...
#cg00000030	1	0	0	1	0	...
#				.
#				.

####################################################
## Step 0 | Install survival package
####################################################
library(survival)

####################################################
## Step 1 | read input file and format the data 
####################################################
a=read.delim("input.txt", sep="\t", header=F)
a=as.matrix(a)
row.names(a)=a[,1]
b=a[-c(1:3),]
t=matrix(0, ncol=2, nrow=dim(b)[1])
t[,1]=b[,1]
a=a[,-1]

####################################################
## Step 2 | calculate survival assocation of each probe (log-rank test)
####################################################
for (i in 4:dim(a)[1])
{
	test=survdiff(Surv(as.numeric(a[3,]),as.numeric(a[2,]))~as.numeric(a[i,]))
	pval=round(1-pchisq(test$chisq, 1),3)
	k=i-3
	t[k,2]=pval
}

####################################################
## Step 3 | write the results 
####################################################
t=as.matrix(t)
write.table(t, file="surv.txt", quote=F, row.names=F, col.names=F, sep="\t")
