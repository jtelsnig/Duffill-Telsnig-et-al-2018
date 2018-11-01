

#Libraries

library(BEST)
library(modeest)
library(devtools)
library(siar)
library(rjags)
library(plyr)
library(coda)

# Example of model using carbon and sulphur isotopes for 2006 data 

modelstring =' #jaggs language
model {
  for(i in 1:N) { 
    for(j in 1:J) {
      y[i,j] ~ dnorm(inprod(p,s_mean[,j]),1/var_y[j])
    }
  }
  p ~ ddirch(alpha)
  for(j in 1:J) {
    var_y[j] <- inprod(pow(p,2),1/s_prec[,j])
    + pow(sigma[j],2)
  }
  for(j in 1:J) { sigma[j] ~ dunif(0,100) }
}
'
northall_consumer<-subset(north, select=c("commonname","TLX34S","TLncorr_d13C")) 
northall_consumer<-northall_consumer[!northall_consumer$commonname %in% "mackerel",] #take out the source indicator species
northall_consumer<-northall_consumer[!northall_consumer$commonname %in% "plaice",]
sources = as.matrix(north_sources[,2:5]) #carbon and sulphur mean and sd for source indicator species 
o<-rep(NA,13)
p0.9l<-rep(NA,13)
p0.9u<-rep(NA,13)
p0.5l<-rep(NA,13)
p0.5u<-rep(NA,13)
pmode<-rep(NA,13)
spnamesn<- as.character(unique(northall_consumer$commonname))
for (i in 1:length(spnamesn)){
specdn<- northall_consumer[northall_consumer$commonname==paste(spnamesn[i]),]
data=list(y=specdn[,2:3],s_mean=sources[,c(1,3)],s_prec=1/sources[,c(2,4)]^2, #data object includes consumers, mean and sd of sources, precision, N is number of consumers, k is the number of sources, and alpha is the number of alpha values 
N=nrow(specdn),J=2,
alpha=rep(1,nrow(sources)))
model_2=jags.model(textConnection(modelstring), data=data) #calling the jags.model and put in the modelstring and data
output=coda.samples(model=model_2,variable.names=c("p",'sigma'),n.iter=100000,nthin=10) #get the output of the model which is the posterior distribution and only saving p at the moment with number of iterations, can also have multiple chains 
o<-output[[1]][,1]
hdi0.9<-hdi(o, credMass=0.9) #using hdi to get the highest density interval
p0.9l[i]<-hdi0.9[1]
p0.9u[i]<-hdi0.9[2]
hdi0.5<-hdi(o, credMass=0.5)
p0.5l[i]<-hdi0.5[1]
p0.5u[i]<-hdi0.5[2]
mode<-mlv(o, method = "HSM")
pmode[i]<-mode$M
}

nsc_credi<-cbind.data.frame(spnamesn,p0.9l,p0.5l,pmode,p0.5u,p0.9u) #data frame with species name and 0.05,0.25,0.5,0.75,0.95 credible intervals 


#Example of model using carbon isotope for 2002 data: 

modelstring =' #jaggs language
model {
  for(i in 1:N) { y[i] ~ dnorm(inprod(p,s),1/pow(sigma,2)) } #likelihood - for consumers of 1:i, it is normally distributed for p and s up to K with some precision (which here is shown as 1/sd^2)
  p ~ ddirch(alpha)
  for(k in 1:K) { s[k] ~ dnorm(s_mean[k],s_prec[k]) } #prior - have four priors with a dirichlet distribution. Each source k has a normal distribution and has a mean k and precision k value
  sigma ~ dunif(0,100) # some uncertainity which is vague at the moment so has a uniform distribution from 0 to 100. 
}
'
y02_consumer<-subset(y02, select=c("commonname","TLncorr_d13C")) #species name and carbon isotope
y02_consumer<-y02_consumer[!y02_consumer$commonname %in% "mackerel",]
y02_consumer<-y02_consumer[!y02_consumer$commonname %in% "plaice",]

o<-rep(NA,13)
p0.9l<-rep(NA,13)
p0.9u<-rep(NA,13)
p0.5l<-rep(NA,13)
p0.5u<-rep(NA,13)
pmode<-rep(NA,13)

spnamesn2<- as.character(unique(y02_consumer$commonname)) 

for (i in 1:length(spnamesn2)){
specdn<- y02_consumer[y02_consumer$commonname==paste(spnamesn2[i]),]
data=list(y=specdn$TLncorr_d13C,s_mean=sources_y02[,1],s_prec=1/sources_y02[,2]^2, #data object includes consumers, mean and sd of sources, precision, N is number of consumers, k is the number of sources, and alpha is the number of alpha values 
N=length(specdn$TLncorr_d13C),K=nrow(sources_y02),
alpha=rep(1,nrow(sources_y02)))
model_2=jags.model(textConnection(modelstring), data=data) #calling the jags.model and put in the modelstring and data
output=coda.samples(model=model_2,variable.names=c("p"),n.iter=100000) #get the output of the model which is the posterior distribution and only saving p at the moment with number of iterations, can also have multiple chains 
o<-output[[1]][,1]
hdi0.9<-hdi(o, credMass=0.9) #using hdi to get the highest density interval
p0.9l[i]<-hdi0.9[1]
p0.9u[i]<-hdi0.9[2]
hdi0.5<-hdi(o, credMass=0.5)
p0.5l[i]<-hdi0.5[1]
p0.5u[i]<-hdi0.5[2]
mode<-mlv(o, method = "HSM")
pmode[i]<-mode$M
}

n06_credi<-cbind.data.frame(spnamesn2,p0.9l,p0.5l,pmode,p0.5u,p0.9u)



