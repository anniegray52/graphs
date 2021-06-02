## packages 
library(data.table) # fread
library(tidyverse) 
library(RSpectra)   # svds
library(igraph)     # graphs
library(ggplot2)    # plotting
library(gridExtra)  # grid.arrange 
library(rgl)        # scatter3D
library(plot3D)     # scatter3D
library(Rtsne)      # Rtsne
library(umap)       # umap
set.seed(123)


## Functions: 
sym <- function(s){s[lower.tri(s)] = t(s)[lower.tri(s)]; s}
doublecentre <- function(B){
  colmeans = matrix(rep(colMeans(B),each=nrow(B)),nrow=nrow(B),ncol=ncol(B), byrow=FALSE)
  rowmeans = matrix(rep(rowMeans(B),each=nrow(B)),nrow=nrow(B),ncol=ncol(B), byrow=TRUE)
  grandmean = matrix(mean(B), nrow=nrow(B), ncol=nrow(B))
  B - colmeans - rowmeans + grandmean
}


## make Z, P and A matrix 
## change m to change n - n = m^2
m = 80
reach = .25
x1 = seq(from=-pi+reach, to=pi-reach, length.out=m+1)[-(m+1)]
Z = expand.grid(x1, x1)
x = Z[,1]
y = Z[,2]
x = expand.grid(x ,x)
y = expand.grid(y ,y)
n = nrow(Z)
Z = as.matrix(Z)
x = matrix(x[,1]-x[,2], nrow=n)
y = matrix(y[,1]-y[,2], nrow=n)
cosx = (cos(x)+1)/2
cosy = (cos(y)+1)/2

##if sparse factor is being used: 
# c = 100/(log(100)^4) 
# rho = c*log(n)^4/n
rho = 1
P = rho*(cosx + cosy)
A = sym(matrix(runif(n^2)<P, nrow=n, ncol=n))
print(paste0("Sparsity: ", 1 - (sum(A, na.rm = TRUE)/n^2)))


## Estimated X 
d = 5 # rank of kernel 
s = svds(A,d)
X = s$u %*% diag(sqrt(s$d))


## to plot first three dimensions of X
pal = c("#9C964A", "#E2D200", "#46ACC8", "#E58601", "#B40F20")
level = matrix(c(rep(1,n/5),rep(2,n/5),rep(3,n/5),rep(4,n/5),rep(5,n/5)), ncol=1)
df = as.data.frame(cbind(X,level))
scatter3D( 
  x=df$V1, y=df$V2, z=df$V3, 
  colvar = df$V6, col = pal, colkey=FALSE,
  pch = 19, cex = 0.5, theta = 180, phi = 45, bty="b2", ticktype="detailed",
  xlab="", ylab="", zlab="")


## create graph
DistX = as.matrix(dist(X))
epsilon = quantile(DistX, 0.006)
newdist=DistX
#set weights equal to distance
newdist[DistX > epsilon] = 0 
G = graph_from_adjacency_matrix(newdist, weighted=TRUE,mode="min")
##check it is connected
print(paste0("Number of graph components: ", components(G)$no))

## isomap
ManifoldDist=distances(G)
InnerProductMatrix = doublecentre(ManifoldDist^2)
s = svds(InnerProductMatrix,k=2)
Xhat = s$u[,1:2]%*% diag(sqrt(s$d[1:2]))

## procrustes
s=svd(t(Xhat)%*%Z)
c=sum(s$d)/sum(diag(Xhat%*%t(Xhat)))
T = s$u %*% t(s$v)
meanXhat=as.matrix(c(mean(Xhat[,1]),mean(Xhat[,2])))
meanZ = as.matrix(c(mean(Z[,1]),mean(Z[,2])))
b = meanZ - T %*% meanXhat
Y = Xhat%*%T
Y = c*(Xhat%*%T)
Y = cbind(Y[,1] + b[1,], Y[,2] + b[2,])


## Z for n = 100
## to find subset
m10 = 10
x_10 = seq(from=-pi+reach, to=pi-reach, length.out=m10+1)[-(m10+1)]
Z10 = expand.grid(x_10, x_10)
x10 = Z[,1]
y10 = Z[,2]
x10 = expand.grid(x10,x10)
y10 = expand.grid(y10,y10)
Z10 = as.matrix(Z10)


## subset 
## find indices needed for subset
index = rep(0, dim(Z10)[1])
for (i in 1:dim(Z10)[1]){
  index[i] = which(apply(Z, 1, function(x) all(x == Z10[i,])))
}
lvl = matrix(c(rep(1,100/5),rep(2,100/5),rep(3,100/5),rep(4,100/5),rep(5,100/5)), ncol=1)
subset =  as.data.frame(cbind(Y[index,],lvl))

## plot subset
fig = ggplot(subset, aes(x = V1, y = V2, colour = factor(V3))) +
  geom_point() + theme_bw() + ggtitle("n=6400") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_color_manual(values = pal)
fig



### approximate geodesic distance and ambient distance
par(mfrow=c(1,2))
dists = c(as.matrix(dist(Z)))
graphdists=c(ManifoldDist)
ambientdists=c(as.matrix(dist(X)))
## take a sample to plot
s = sample(length(dists), 5000)
dists = dists[s]
graphdists = graphdists[s]
ambientdists = ambientdists[s]

man_amb = data.frame(cbind(dists, graphdists, ambientdists))
man = ggplot(man_amb, aes(x = dists, y = graphdists)) +
  geom_point(size = .5, color = "grey4") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))  +
  xlab("Latent distance") + ylab("Approximate geodesic distance") + theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14)) + theme(legend.position = "none")
amb = ggplot(man_amb, aes(x = dists, y = ambientdists)) +
  geom_point(size = .5, color = "grey4") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))   +
  xlab("Latent distance") + ylab("Ambient distance") + theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14)) + theme(legend.position = "none")
manifold_vs_ambient = grid.arrange(man, amb, ncol = 2)


## Alternative methods
## all to be run in the sparse regime


## t-SNE
tsne <- Rtsne(X, dims = 2, perplexity=50, verbose=TRUE, max_iter = 1000)
Xhat = tsne$Y
# procrustes
s=svd(t(Xhat)%*%Z)
c=sum(s$d)/sum(diag(Xhat%*%t(Xhat)))
T = s$u %*% t(s$v)
meanXhat=as.matrix(c(mean(Xhat[,1]),mean(Xhat[,2])))
meanZ = as.matrix(c(mean(Z[,1]),mean(Z[,2])))
b = meanZ - T %*% meanXhat
Y = Xhat%*%T
Y = c*(Xhat%*%T)
Y = cbind(Y[,1] + b[1,], Y[,2] + b[2,])
tsneY = Y


## UMAP 
um = umap(X)
Xhat = um$layout
## procrustes
s=svd(t(Xhat)%*%Z)
c=sum(s$d)/sum(diag(Xhat%*%t(Xhat)))
T = s$u %*% t(s$v)
meanXhat=as.matrix(c(mean(Xhat[,1]),mean(Xhat[,2])))
meanZ = as.matrix(c(mean(Z[,1]),mean(Z[,2])))
b = meanZ - T %*% meanXhat
Y = Xhat%*%T
Y = c*(Xhat%*%T)
Y = cbind(Y[,1] + b[1,], Y[,2] + b[2,])
umapY = Y



## node2vec
# A=A*1
# fwrite(A, file = "sparse_adj_mat.csv")
## just into 5 dimensions - for 2d don't use isomap
n2v5 = fread('n2v_n5.csv')
DistX = as.matrix(dist(n2v5))
epsilon = quantile(DistX, 0.02)
newdist=DistX
newdist[DistX > epsilon] = 0
g1 = graph_from_adjacency_matrix(newdist, weighted=TRUE,mode="min")
print(paste0("Number of graph components: ", components(g1)$no))
## isomap
N2vDist=distances(g1)
InnerProductMatrix = doublecentre(N2vDist^2)
s = svds(InnerProductMatrix,k=2)
Xhat = s$u[,1:2]%*% diag(sqrt(s$d[1:2]))
## procrustes
s=svd(t(Xhat)%*%Z)
c=sum(s$d)/sum(diag(Xhat%*%t(Xhat)))
T = s$u %*% t(s$v)
meanXhat=as.matrix(c(mean(Xhat[,1]),mean(Xhat[,2])))
meanZ = as.matrix(c(mean(Z[,1]),mean(Z[,2])))
b = meanZ - T %*% meanXhat
Y = Xhat%*%T
Y = c*(Xhat%*%T)
Y = cbind(Y[,1] + b[1,], Y[,2] + b[2,])
n2vY = Y


## graph distance method 
## input sparse adjacency matrix
A = fread('sparse_adj_mat.csv')
g2 = graph_from_adjacency_matrix(as.matrix(A))
## cmds
GraphDist=distances(g2)
InnerProductMatrix = doublecentre(GraphDist^2)
s = svds(InnerProductMatrix,k=2)
Xhat = s$u[,1:2]%*% diag(sqrt(s$d[1:2]))
## procrustes
s=svd(t(Xhat)%*%Z)
c=sum(s$d)/sum(diag(Xhat%*%t(Xhat)))
T = s$u %*% t(s$v)
meanXhat=as.matrix(c(mean(Xhat[,1]),mean(Xhat[,2])))
meanZ = as.matrix(c(mean(Z[,1]),mean(Z[,2])))
b = meanZ - T %*% meanXhat
Y = Xhat%*%T
Y = c*(Xhat%*%T)
Y = cbind(Y[,1] + b[1,], Y[,2] + b[2,])
adjgraphdistsY = Y


## plotting alternative methods
## take subset and plot on a grid
n2v_subset = as.data.frame(cbind(n2vY[index,], lvl))
tsne_subset = as.data.frame(cbind(tsneY[index,], lvl))
umap_subset = as.data.frame(cbind(umapY[index,], lvl))
graphdists_subset = as.data.frame(cbind(adjgraphdistsY[index,], lvl))


fig1 = ggplot(n2v_subset, aes(x=V1, y=V2,colour = factor(V3))) + 
  geom_point() + theme_bw() + ggtitle("node2vec (10 dimensions), followed by Isomap")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  ) +
  scale_color_manual(values=pal)
fig2 = ggplot(tsne_subset, aes(x=V1, y=V2,colour = factor(V3))) + 
  geom_point() + theme_bw() + ggtitle("Spectral embedding, followed by t-SNE")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  ) +
  scale_color_manual(values=pal)
fig3 = ggplot(umap_subset, aes(x=V1, y=V2,colour = factor(V3))) + 
  geom_point() + theme_bw() + ggtitle("Spectral embedding, followed by UMAP")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  ) +
  scale_color_manual(values=pal)
fig4 = ggplot(graphdists_subset, aes(x=V1, y= V2,colour = factor(V3))) + 
  geom_point() + theme_bw() + ggtitle("Graph distance, followed by CMDS")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  ) +
  scale_color_manual(values=pal)
grid.arrange(fig1, fig2, fig3, fig4, nrow=2)