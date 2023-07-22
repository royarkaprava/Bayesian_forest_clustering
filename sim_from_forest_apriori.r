# simulate y from the forest process, using the ground-up distribution and scale parameters simulated from the shrinkage prior
y<- matrix(0,nrow=0,ncol=2)

n<- 100

i=1
alpha=0.5

gamma = sqrt(1/rgamma(1,2,1))

a_sigma = 100
b_sigma = 10

xi_sigma = 1
eta_sigma = 1/rgamma(1,a_sigma,xi_sigma)
beta_sigma = rexp(1, rate=1/eta_sigma)

y<- rbind(y,rt(2,df=1)*gamma)
sigma<- numeric()
sigma<- c(sigma, 1/rgamma(1,b_sigma,beta_sigma))

edges<- numeric()
edges<- c(edges,1)

for(i in 2:n){
  sigma<- c(sigma, 1/rgamma(1,b_sigma,beta_sigma))
  
  j <- sample(i,1, prob= c(rep(1,i-1),alpha))
  edges<- c(edges,j)
  
  if(j == i){
    y<- rbind(y,rt(2,df=1)*gamma)
  }else{
    sigma_sd<- sqrt(sigma[i]*sigma[j])
    y<- rbind(y,rnorm(2,0,sd = sigma_sd)+y[j,])
  }
}

png("randomly_gen_data_shrunk_sigma.png",4,4,units = "in",res =300)
plot(y,xlab="",ylab="")
for(i in 1:n){
  if(i != edges[i]){
    lines(y[c(i,edges[i]),1],y[c(i,edges[i]),2])
  }
}
dev.off()


# simulate y from the forest process, using the ground-up distribution and scale parameters fixed to small value
y<- matrix(0,nrow=0,ncol=2)

n<- 100

i=1
alpha=0.5

gamma = 0.1
sigma_sd =  0.01


y<- rbind(y,rt(2,df=1)*gamma)
sigma<- numeric()
sigma<- c(sigma, 1/rgamma(1,b_sigma,beta_sigma))
edges<- numeric()
edges<- c(edges,1)

for(i in 2:n){
  sigma<- c(sigma, 1/rgamma(1,b_sigma,beta_sigma))
  
  j <- sample(i,1, prob= c(rep(1,i-1),alpha))
  edges<- c(edges,j)
  
  if(j == i){
    y<- rbind(y,rt(2,df=1)*gamma)
  }else{
    y<- rbind(y,rnorm(2,0,sd = sigma_sd)+y[j,])
  }
}



png("randomly_gen_data_fixed_sigma.png",4,4,units = "in",res =300)
plot(y,xlab="",ylab="")
for(i in 1:n){
  if(i != edges[i]){
   lines(y[c(i,edges[i]),1],y[c(i,edges[i]),2])
  }
}
dev.off()
