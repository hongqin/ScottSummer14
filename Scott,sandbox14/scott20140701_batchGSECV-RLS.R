# batch_GSELog2CV - RLS regression, 2014 July 1

setwd("~/github/ScottSummer14/Scott,sandbox14")
rm(list=ls())

#RLS table
RLS.tb = read.csv("data/lifespan.csv")
RLS.tb$ORF = as.character(RLS.tb$ORF)
str(RLS.tb)
#'data.frame':  584 obs. of  17 variables:

files = list.files(, path="data/GSElog2CV", pattern='GSE') 

output = data.frame(files)

# i = 1
for ( i in 1:length(files) ) {
  filename = paste( "data/GSElog2CV/", files[i], sep='')
  CV.tb = read.csv(filename)  
  CV.tb$ORF = as.character( CV.tb$ORF )

  RLS.tb2 = merge(RLS.tb, CV.tb)
  
  m = lm(RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myCV)  
  s = summary(m)
  output$p.RLS.CV[i] = 1-pf( s$fstat[1], s$fstat[2], s$fstat[3])
  
  m = lm( log(RLS.tb2$RLS_Del_alpha) ~ RLS.tb2$myCV)  
  s = summary(m)
  output$p.logRLS.CV[i] = 1-pf( s$fstat[1], s$fstat[2], s$fstat[3])
 
  m = lm( RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myMean)  
  s = summary(m)
  output$p.RLS.GSEMean[i] = 1-pf( s$fstat[1], s$fstat[2], s$fstat[3])

  m = lm( RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myStddev)  
  s = summary(m)
  output$p.RLS.GSEStddev[i] = 1-pf( s$fstat[1], s$fstat[2], s$fstat[3])
  
}

summary(output)
output[ output$p.RLS.GSEStddev<0.05, ]
output[ output$p.RLS.GSEMean<0.05, ]

