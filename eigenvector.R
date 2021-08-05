y<- c(rnorm(30),rnorm(30)+ 0.1)

p<- length(y)

hist(y)

D<- abs(outer(y,y,"-"))


W<- exp(-D)


diag(W)=0
Dw = rowSums(W)

Lw = diag(Dw) - W  

# normalized Laplacian
Nw = diag(sqrt(1/Dw))%*% Lw %*%diag(sqrt(1/Dw))

# Marginal inclusion
Omega = solve(Lw + 1/p + 1E-12* diag(1,p))

DOmega1 = matrix(diag(Omega),p,p)

M = (DOmega1+ t(DOmega1)-2*Omega)*W

#eigendecomp of M, vectors sorted in descending order of values
ev_decomp_M = eigen(M)

#eigendecomp of Nw, vectors sorted in ascending order of values
ev_decomp_Nw = eigen(Nw)
ev_decomp_Nw$vectors = ev_decomp_Nw$vectors[,order(ev_decomp_Nw$values)]


align<- function(x,y){
  a = sum(abs(x-y))
  b = sum(abs(x+y))
  ifelse(a<b,1,-1)
}


#par(mfrow=c(3,2))
for( i in 1:p){
  #plot(ev_decomp_M$vectors[,i], main=i)
  # align the signs
  s = align(ev_decomp_M$vectors[,i], ev_decomp_Nw$vectors[,i])
  #lines(s*ev_decomp_Nw$vectors[,i])
  print(sum((ev_decomp_M$vectors[,i]-s*ev_decomp_Nw$vectors[,i])^2))
}