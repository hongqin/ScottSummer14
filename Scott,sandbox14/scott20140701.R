setwd("~/github/ScottSummer14/Scott,sandbox14")
rm(list=ls())

#RLS table
RLS.tb = read.csv("data/lifespan.csv")
RLS.tb$ORF = as.character(RLS.tb$ORF)
str(RLS.tb)
#'data.frame':  584 obs. of  17 variables:

#these are robustness from gene expression
#robustness is inversely represented by CV=signma/mean
list.files(, path="data", pattern='GSE') 
#GSE1221 CV !~ RLS
#GSE18334_log2CV.csv CV~ RLS p=0.45
#GSE33276_log2CV.csv, p=0.9
# GSE3821_log2CV.csv, p=0.97
#GSE9840_log2CV.csv, p=0.052

CV.tb = read.csv("data/GSE9985_log2CV.csv")  #CV table  
CV.tb$ORF = as.character( CV.tb$ORF )
str(CV.tb)

### merged RLS, CV based on shared 'ORF' column
RLS.tb2 = merge(RLS.tb, CV.tb)

m = lm(RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myCV)  
summary(m)

m = lm(RLS.tb2$RLS_Del.pooled  ~ RLS.tb2$myCV)
summary(m)

plot( RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myCV  )
abline(m, col='red')
#large CV or Stddev -> more noisy -> less robust -> slower aging rate -> long RLS (?!)

#######################
#Now, let's look at morphological robustness

#Read cellular morpholy mutant table with stddev, mean, CV into R
# scmd_CVbyRow20140701.csv was an output of "ScottSummer14/scmd_CV_20140701.R"
scmd = read.csv( "scmd_CVbyRow20140701.csv")
names(scmd)[c(1:3, 500:505)]; #check the column names

#scmd$CV = scmd$stddev / (scmd$mean - min(scmd$mean) + 0.01) #make sure CV is positive
summary(scmd[,c("CV", "stddev", "mean")])


############## Analyze morphology robustness, lifespan, and gene expression robustness
# merge scmd with RLS.tb2 
#lifespan$scmdstddev = scmd$stddev[match(lifespan$ORF, scmd$name)];
RLS.tb2$scmdCV = scmd$CV[match(RLS.tb2$ORF, scmd$name)];
RLS.tb2$scmdMean = scmd$mean[match(RLS.tb2$ORF, scmd$name)];
RLS.tb2$scmdStddev = scmd$stddev[match(RLS.tb2$ORF, scmd$name)];

#check the merged results
#head(RLS.tb2)

#now, regression analysis
m = lm( RLS.tb2$RLS_Del_alpha ~ RLS.tb2$scmdCV )
summary(m)

m = lm( RLS.tb2$RLS_Del_alpha ~ sqrt(RLS.tb2$scmdCV ))
summary(m)


m = lm( RLS.tb2$RLS_Del_alpha ~ RLS.tb2$scmdMean )
summary(m)

m = lm( RLS.tb2$RLS_Del_alpha ~ RLS.tb2$scmdStddev )
summary(m)

summary( lm( RLS.tb2$myCV ~ sqrt(RLS.tb2$scmdCV )) )#p 0.78
summary( lm( RLS.tb2$myCV ~ sqrt(RLS.tb2$scmdMean )) )#p 0.03 ?!
summary( lm( RLS.tb2$myCV ~ sqrt(RLS.tb2$scmdStddev )) )#p 0.02 ?!


