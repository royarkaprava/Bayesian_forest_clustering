if(Sys.info()["sysname"] %in% c("Linux","Darwin") ){
  setwd("~/git/Spectral-clustering/")
}

if(Sys.info()["sysname"] %in% c("Windows") ){
  setwd("\\\\file.phhp.ufl.edu\\home\\ark007\\My Documents\\GitHub\\Spectral-clustering")
}




Rcpp::sourceCpp('updateTreedefusednormal.cpp')

data <- read.csv("./data.txt", header = F)

loc <- as.matrix(data[, 1:2])
label <- data$V3
label[which(label=="Class 1")] <- 1
label[which(label=="Class 2")] <- 2
label <- as.numeric(label)
label

ind1 <- which(label==1)
ind2 <- which(label==2)
idx <- c(sample(ind1,100), sample(ind2, 100))

locR <- loc[idx, ]
label_set = label[idx]
X    <- locR

#X <- X - rowMeans(X)
######################################################################################

#I am initializing using ComputeMST(), there may be other ways too.

a0 =  1
b0 =  1E-3

# source('forestProcessDirichlet.R')
# fit <- clusteringFP(X,alpha0 =  0.001, a0= a0, b0 = b0, K= 20, Total_itr = 1000, burn=500, gamma=1000)

source('forestProcessQuasiBernoullidefusednormal.R')
fit <- clusteringFP(X,p_b =  0.1, a0= a0, b0 = b0, K= 10,lam=10, Total_itr = 1000, burn=500, gamma=1000,random_scan_n = 0)


ts.plot(fit$sig_ls)
acf(fit$sig_ls,lag.max = 40)

min(dist(X)[dist(X)>0])

rowMeans(fit$clssize_ls) #estimated class sizes
#fit$estiadja #estimated adjacency
rowMeans(fit$stickbrkwts) #Estimated components for stick breaking


barplot(table(apply(fit$clslb_ls,2, function(x)length(unique(x)))))

image(fit$estiadja,col = topo.colors(100, rev=F))

plot(fit$estiadja[1,])
# # fit$estiadja
# 
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

modelabel <- apply(fit$clslb_ls, 1, getmode)

plot(X[,1],X[,2], col = modelabel)

plot(X[,1],X[,2], col = fit$clslb_ls[,100])
plot(X[,1],X[,2], col = fit$clslb_ls[,500])
plot(X[,1],X[,2], col = fit$clslb_ls[,300])
