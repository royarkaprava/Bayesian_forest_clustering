

require("fields")

require("matlab")

# plot(eta_star)

load("fMRI_application/fmri_labels.RDa")


CN_idx <-  (1:S)[labels=="CN"]
LCMI_idx <-  (1:S)[labels=="LMCI" |labels=="LMCI3"]



num_class_CN<- sapply(CN_idx,function(s) length(unique((extractC(forest_objs[[s]]$A_T)))))
num_class_LCMI<- sapply(LCMI_idx,function(s) length(unique((extractC(forest_objs[[s]]$A_T)))))



Z_mat<- sapply(forest_objs, function(x){c(x$Z)})
Z_dist<- as.matrix(dist(t(Z_mat),upper = TRUE, diag=TRUE))


# png("z_dist.png",600,600,res = 100)
par(cex.axis=2.0)

# Z_dist[Z_dist<2.5]<- 2.5

Z_dist<- t(apply(Z_dist,1,rev))
# Z_dist<- apply(Z_dist,2,rev)

image.plot(Z_dist,col = jet.colors(30), xaxt= "n", yaxt= "n" )
axis(1, at = seq(0, 160/S, length.out =  5), labels= c(0:4)*40, las=1,
     cex.lab=30)
axis(2, at = 1-seq(0, 160/S, length.out =  5), labels= c(0:4)*40, las=1,
     cex.lab=3
)

# dev.off()



node_wise_total_var <- sapply(c(1:n),function(i){
  
  l1<- c(1:S)
  Z_array<- t(sapply(l1, function(x){
    forest_objs[[x]]$Z[i,]
  }))
  SS1<- sum(sum((t(Z_array) - colMeans(Z_array))**2))
  
  SS1/S/2
})



node_wise_resid_var <- sapply(c(1:n),function(i){
  
  l1<- CN_idx
  Z_array<- t(sapply(l1, function(x){
    forest_objs[[x]]$Z[i,]
  }))
  SS1 <-  sum(sum((t(Z_array) - colMeans(Z_array))**2))
  
  l2<- LCMI_idx
  Z_array<- t(sapply(l2, function(x){
    forest_objs[[x]]$Z[i,]
  }))
  SS2 <-  sum(sum((t(Z_array) - colMeans(Z_array))**2))
  
  (SS1+SS2)/S/2
})



node_wise_resid_var_CN <- sapply(c(1:n),function(i){
  
  l1<- CN_idx
  Z_array<- t(sapply(l1, function(x){
    forest_objs[[x]]$Z[i,]
  }))
  SS1 <-  sum(sum((t(Z_array) - colMeans(Z_array))**2))
  
  SS1/nrow(Z_array)/2
})



node_wise_resid_var_LCMI <- sapply(c(1:n),function(i){
  
  l1<- LCMI_idx
  Z_array<- t(sapply(l1, function(x){
    forest_objs[[x]]$Z[i,]
  }))
  SS1 <-  sum(sum((t(Z_array) - colMeans(Z_array))**2))
  
  SS1/nrow(Z_array)/2
})







node_wise_resid_var_LCMI <- sapply(c(1:n),function(i){
  
  l1<- LCMI_idx
  Z_array<- t(sapply(l1, function(x){
    forest_objs[[x]]$Z[i,]
  }))
  SS1 <-  sum(sum((t(Z_array) - colMeans(Z_array))**2))
  
  SS1/nrow(Z_array)/2
})



coords<- read.table("fMRI_application/centroid.txt")



hippocampus_sel <- c(1:n)[(coords$V3<60 &coords$V3>50 & coords$V4<35) | (coords$V3<23 & coords$V4>49)]
pcl_sel <- c(1:n)[(coords$V3>0 & coords$V4>63) | (coords$V3>75 & coords$V4>55)]



node_wise_hippocampus_cn <- lapply(hippocampus_sel,function(i){
  
  l1<- CN_idx
  Z_array<- t(sapply(l1, function(x){
    forest_objs[[x]]$Z[i,]
  }))
  
  Z_array
})

node_wise_hippocampus_lcmi <- lapply(hippocampus_sel,function(i){
  
  l1<- LCMI_idx
  Z_array<- t(sapply(l1, function(x){
    forest_objs[[x]]$Z[i,]
  }))
  
  Z_array
})


node_wise_pcl_cn <- lapply(pcl_sel,function(i){
  
  l1<- CN_idx
  Z_array<- t(sapply(l1, function(x){
    forest_objs[[x]]$Z[i,]
  }))
  
  Z_array
})

node_wise_pcl_lcmi <- lapply(pcl_sel,function(i){
  
  l1<- LCMI_idx
  Z_array<- t(sapply(l1, function(x){
    forest_objs[[x]]$Z[i,]
  }))
  
  Z_array
})



require("reshape")
require("ggplot2")



# hippocampus regions
roi_names = c("HIP.L","PHG.L","SOG.L", "SOG.R", "HIP.R", "PHG.R")

df1<- sapply(c(1:6), function(k){ 
  mat<- node_wise_hippocampus_cn[[k]]
  rowSums(mat)
})
colnames(df1)<- roi_names

df2<- sapply(c(1:6), function(k){ 
  mat<- node_wise_hippocampus_lcmi[[k]]
  rowSums(mat)
})
colnames(df2)<- roi_names


df1<- melt(df1)
df1$X1<- "Healthy"

df2<- melt(df2)
df2$X1<- "Diseased"

df3<- rbind(df1,df2)

colnames(df3)<- c("Group","Region","Value")

df3$Region<- factor(df3$Region, levels = c("HIP.L", "HIP.R","PHG.L", "PHG.R", "SOG.L", "SOG.R"))
df3$Group<- factor(df3$Group, levels = c("Healthy","Diseased"))

# png("brain_boxplot_hippo.png",1000,600,res = 100)

ggplot(df3, aes(x = Region, y = Value, fill = Group)) +
  # geom_violin() +
  geom_boxplot(alpha = 0.5, outlier.shape = "") +
  scale_fill_manual(values = c("dodgerblue2", "orange")) +
  xlab("Region of Interest") +
  ylab("Value")+theme_bw()+ theme(
    text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust = 1)
  )

# dev.off()



validation_X1<- rbind(do.call("cbind",node_wise_hippocampus_cn), do.call("cbind",node_wise_hippocampus_lcmi))
validation_X2<- rbind(do.call("cbind",node_wise_pcl_cn), do.call("cbind",node_wise_pcl_lcmi))



validation_X<- cbind(validation_X1,validation_X2)
validation_y<- c( rep(1, length(CN_idx)),rep(0, length(LCMI_idx)))



library(pROC)



glm_fit <- glm(validation_y ~ validation_X, family = binomial)
prob <- predict(glm_fit, type = "response")
roc <- roc(validation_y, prob)
auc <- auc(roc)

auc



glm_fit <- glm(validation_y ~ validation_X1, family = binomial)
prob <- predict(glm_fit, type = "response")
roc <- roc(validation_y, prob)
auc <- auc(roc)

auc



glm_fit <- glm(validation_y ~ validation_X2, family = binomial)
prob <- predict(glm_fit, type = "response")
roc <- roc(validation_y, prob)
auc <- auc(roc)

auc



library(glmnet)


validation_full_z<- t(Z_mat)

glm_fit <- glmnet(validation_full_z,validation_y, family = binomial, alpha = 1,intercept = TRUE)




prob <- predict(glm_fit, newx = validation_full_z, type = "response")




roc <- roc(validation_y, prob[,23])
auc <- auc(roc)

auc







# pcl regions
roi_names = c("SPG.L","SFGdor.L","PCL.L","SMA.L","PCL.R" ,"SMA.R", "SPG.R", "SFGdor.R")

df1<- sapply(c(1:8),function(k){ 
  mat<- node_wise_pcl_cn[[k]]
  rowSums(mat)
})
colnames(df1)<- roi_names

df2<- sapply(c(1:8),  function(k){ 
  mat<- node_wise_pcl_lcmi[[k]]
  rowSums(mat)
})
colnames(df2)<- roi_names


df1<- melt(df1)
df1$X1<- "Healthy"

df2<- melt(df2)
df2$X1<- "Diseased"

df3<- rbind(df1,df2)

colnames(df3)<- c("Group","Region","Value")

df3$Region<- factor(df3$Region, levels = c("PCL.L","PCL.R" , "SPG.L","SPG.R", "SFGdor.L","SFGdor.R","SMA.L","SMA.R"))
df3$Group<- factor(df3$Group, levels = c("Healthy","Diseased"))

# png("brain_boxplot_pcl.png",1000,600,res = 100)

ggplot(df3, aes(x = Region, y = Value, fill = Group)) +
  # geom_violin() +
  geom_boxplot(alpha = 0.5, outlier.shape = "") +
  scale_fill_manual(values = c("dodgerblue2", "orange")) +
  xlab("Region of Interest") +
  ylab("Value")+theme_bw()+ theme(
    text = element_text(size = 30),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# dev.off()



full_roi_names<- read.table("fMRI_application/regionname.txt")



library("ggrepel")  



require("brainGraph")



coords_aal<- cbind(aal116$name,aal116$x.mni,aal116$y.mni,aal116$z.mni)



coords_aal[,2:4]<- as.numeric(coords_aal[,2:4])



df<- data.frame(coords_aal)



node_wise_icc<-  1 - node_wise_resid_var/node_wise_total_var

# png("node_z_brain_icc.png",700,600,res = 100)

ggplot()+ geom_point(aes(x=coords$V3,y=coords$V4,col=node_wise_icc),size=5)+
  theme_bw()+ theme(legend.title=element_blank(),
                    text = element_text(size = 30),
  )+ xlab("") + ylab("")  +
  scale_color_gradientn(colors = jet.colors(60),limits=c(0., .4))+
  theme(plot.margin = margin(1,1,1,1, "cm"))

# dev.off()



# png("node_z_var_brain_CN.png",700,600,res = 100)

ggplot()+ geom_point(aes(x=coords$V3,y=coords$V4,col= (node_wise_resid_var_CN)),size=5)+
  theme_bw()+ theme(legend.title=element_blank(),
                    text = element_text(size = 30),
  )+ xlab("") + ylab("")  +
  scale_color_gradientn(colors = jet.colors(60),limits=c(0,0.7))+
  theme(plot.margin = margin(1,1,1,1, "cm"))

# dev.off()



# png("node_z_var_brain_LCMI.png",700,600,res = 100)

ggplot()+ geom_point(aes(x=coords$V3,y=coords$V4,col= (node_wise_resid_var_LCMI)),size=5)+
  theme_bw()+ theme(legend.title=element_blank(),
                    text = element_text(size = 30),
  )+ xlab("") + ylab("")  +
  scale_color_gradientn(colors = jet.colors(60) ,limits=c(0., 0.7))+
  theme(plot.margin = margin(1,1,1,1, "cm"))

# dev.off()


Z_dist<- as.matrix(dist(t(Z_mat),upper = TRUE, diag=TRUE))



diag(Z_dist)<- Inf



Z_dist_cn<- Z_dist[CN_idx,CN_idx]



sel_subs_cn<- CN_idx[apply(Z_dist_cn== min(Z_dist_cn),1,any)]

