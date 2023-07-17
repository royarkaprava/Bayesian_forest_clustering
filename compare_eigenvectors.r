sep_b =2

n <- 60
p <- 2

y<- matrix(rnorm(n*p),n,p)

y[21:40,]<- y[21:40,]+ sep_b
y[41:60,]<- y[41:60,]+ sep_b*2

plot(y[,1],y[,2],xlab="",ylab="")


D<- as.matrix(dist(y, diag=TRUE,upper = TRUE) )

# get the weighted adjacency matrix
A<- exp(-D**2/2)/(2*pi)
diag(A)<- 0 
degree<- rowSums(A)

L<- diag(degree) - A

# get the marginal connecting probability
Omega<- solve(L+1/n/n)
Omega_ii<- diag(Omega)
M<- (outer(Omega_ii,Omega_ii,"+")-2*Omega)*A


# get the normalized Laplacian
N<- diag(degree**(-0.5))%*%L%*%diag(degree**(-0.5))


ev_M<- eigen(M)$vector[,1:4]
ev_N<- eigen(-N)$vector[,1:4]

# compare the eigenvectors side-by-side
#dev.new()
par(mfrow=c(3,2))
mainvec <- c("Eigen vector 1","Eigen vector 2","Eigen vector 3")
for(i in 1:3){
  # check sign:
  if (ev_M[,i]%*%ev_N[,i]<0)  ev_N[,i]<- -ev_N[,i]
  if(i==1){
    plot(ev_M[,i], type = 'l', main = "M", ylab = mainvec[i],xlab="")
    plot(ev_N[,i], type = 'l', main = "N", ylab="",xlab="")
  }
  if(i>1){
    plot(ev_M[,i], type = 'l', ylab = mainvec[i],xlab="")
    plot(ev_N[,i], type = 'l', ylab = "",xlab="")
  }
}
