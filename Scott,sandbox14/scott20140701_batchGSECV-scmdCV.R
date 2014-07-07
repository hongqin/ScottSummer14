# batch GSE stddev ~ scmdStddev, CV, 2014 July 1

setwd("~/github/ScottSummer14/Scott,sandbox14")
rm(list=ls())

#Read cellular morpholy mutant table with stddev, mean, CV into R
# scmd_CVbyRow20140701.csv was an output of "ScottSummer14/scmd_CV_20140701.R"
scmd = read.csv( "scmd_CVbyRow20140701.csv")
names(scmd)[c(1:2, 500:506)]
scmd$ORF = as.character(scmd$name)

#GSE files 
files = list.files(, path="data/GSElog2CV", pattern='GSE') 

output = data.frame(files)

# i = 15
for ( i in 1:length(files) ) {
  filename = paste( "data/GSElog2CV/", files[i], sep='')
  CV.tb = read.csv(filename)  
  CV.tb$ORF = as.character( CV.tb$ORF )

  tb = merge(CV.tb, scmd[, c(1:2,502:506)] )#based on ORF 

  m = lm(tb$myCV ~ tb$CV )  
  s = summary(m)
  output$p.gseCV.scmdCV[i] = 1-pf( s$fstat[1], s$fstat[2], s$fstat[3])
  
  m = lm(tb$myStddev ~ tb$stddev)  
  s = summary(m)
  output$p.gseStddev.scmdStddev[i] = 1-pf( s$fstat[1], s$fstat[2], s$fstat[3])

  picname = paste( 'plots/', files[i], '.stddev.pdf', sep='' )
  pdf(picname, width=5, height=5)
  plot(tb$myStddev ~ tb$stddev, xlab='morphology plasticity', ylab='expression plasticity', main=files[i])
  myY = max(tb$myStddev)*0.8; myX=max(tb$stddev)*0.7;
  mytext = paste("p=",1-pf( s$fstat[1], s$fstat[2], s$fstat[3]) );
  text(myX,myY, mytext)
  abline(m, col='red')
  dev.off()  
}


summary(output)
hist(output$p.gseCV.scmdCV, br=20)
hist(output$p.gseStddev.scmdStddev, br=20)
tmp = output[ output$p.gseStddev.scmdStddev<0.05, ]
length(tmp[,1])

# i = 8

