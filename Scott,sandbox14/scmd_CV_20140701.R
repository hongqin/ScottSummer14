
#######################
#Now, let's look at morphological robustness

#Read cellular morpholy mutant table into R
scmd = read.table( "data/scmd.tab", sep="\t", header=T)

#Normalize the data by column. This took <1min in my apple laptop "byte"
row.names(scmd) = as.character( scmd$name )
for( j in 2:502 ){
  scmd[,j] = ( scmd[,j] - mean(scmd[,j],na.rm=T) )/ sqrt( var( scmd[,j], na.rm=T))
};
head(scmd)


#Calculate sigma, standard deviation, CV by row. This took about 10 min in apple laptop. 
for ( i in 1:4718){
  scmd$stddev[i] = sqrt(var( t(scmd[i,2:502]), na.rm=T ) )
  scmd$mean[i] = mean( t(scmd[i,2:502]), na.rm=T ) 
}; 
scmd$CV = scmd$stddev / (scmd$mean - min(scmd$mean) + 0.01) #make sure CV is positive
head(scmd[,c("CV", "stddev", "mean")])
summary(scmd[,c("CV", "stddev", "mean")])

names(scmd)[c(1:3, 500:505)]; #check the column names

write.csv(scmd, "scmd_CVbyRow20140701.csv", quote=F, row.names=F)
