library(data.table)
library(igraph)
library(Matrix)
library(irlba)
library(RSpectra)

doublecentre <- function(B){
  colmeans = matrix(rep(colMeans(B),each=nrow(B)),nrow=nrow(B),ncol=ncol(B), byrow=FALSE)
  rowmeans = matrix(rep(rowMeans(B),each=nrow(B)),nrow=nrow(B),ncol=ncol(B), byrow=TRUE)
  grandmean = matrix(mean(B), nrow=nrow(B), ncol=nrow(B))
  B - colmeans - rowmeans + grandmean
}

###download airports.csv here: https://ourairports.com/data/
rawairports = fread("airports.csv")
airports = rawairports[,c("iso_country", "gps_code","continent")]
airports = airports[gps_code!="",]
airports=airports[!duplicated(gps_code),]
setkey(airports, gps_code)

D=10
qepsilon = 0.05

rawdata=fread("flightlist_20200101_20200131.csv.gz")
data=rawdata[(rawdata$origin!="") & (rawdata$destination!=""),]
with_origin_info = merge(data, airports, by.x = "origin", by.y = "gps_code"); setnames(with_origin_info,"iso_country","origin_country");setnames(with_origin_info,"continent","origin_continent")
datax = merge(with_origin_info, airports, by.x = "destination", by.y = "gps_code"); setnames(datax,"iso_country","destination_country"); setnames(datax,"continent","destination_continent")
dataC = datax[,.(count = .N, origin = unique(origin), destination = unique(destination)), by = paste(pmin(origin, destination), pmax(origin, destination))]
all = unique(c(dataC$origin,dataC$destination))
ids=tapply(all,all)
names(ids)=all
id2names = names(sort(ids[!duplicated(ids)]))
E1 = ids[dataC$origin]
E2 = ids[dataC$destination]

A = sparseMatrix(i=c(E1,E2), j=c(E2,E1))
rownames(A) = id2names
e <- irlba::irlba(A, D)
X <- e$u %*% diag(sqrt(e$d))
X = X/sqrt(rowSums(X^2))

DistX = as.matrix(dist(X))
epsilon = quantile(DistX, qepsilon)
newdist=DistX
newdist[DistX > epsilon] = 0

epsilongraph = graph_from_adjacency_matrix(newdist, weighted=TRUE,mode="min")
components=clusters(epsilongraph)
largest=which.max(components$csize)
id2names[which(components(epsilongraph)$membership!=largest)]
egraph=induced.subgraph(epsilongraph, V(epsilongraph)[components$membership==largest])

##approximate geodesic distances (this could take ~30 mins)
ManifoldDist=distances(egraph)    

retained = as.numeric(names(V(egraph)))
rownames(ManifoldDist) = id2names[retained]; colnames(ManifoldDist) = id2names[retained]
Xdenoised = X[retained,]; rownames(Xdenoised)= id2names[retained]

##CMDS on graph
InnerProductMatrix = doublecentre(ManifoldDist^2)
s = svds(InnerProductMatrix,k=2)
Xhat = s$u[,1:2]%*% diag(sqrt(s$d[1:2]))
rownames(Xhat) = rownames(ManifoldDist)

##plotting
continent_colours=c("#FDAF91FF" ,"#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#00468BFF")
names(continent_colours) = c("NA","EU","AS","SA","AF","OC")
selected_cities=c("Los Angeles"="KLAX", "Chicago"="KORD", "Dallas"="KDFW", "New York"="KJFK","Dubai"="OMDB", "Paris"="LFPG", "London"="EGLL", "Rome"="LIRF", "Madrid"="LEMD", "Lisbon"="LPPT", "Athens"="LGAV", "Tokyo"="RJTT", "Sao Paolo"="SBGR","Kempton Park"="FAOR","Mexico City"="MMMX","Auckland"="NZAA", "Buenos Aires"="SAEZ")
continents=airports[rownames(Xhat), continent]

plot(Xhat[,1:2], col=continent_colours[continents], cex=.3,pch=16, xlab="", ylab="", axes=FALSE, main="Spectral embedding, followed by Isomap (Jan 20)")
selected_cities_date = selected_cities[which(selected_cities %in% rownames(Xhat))]
text(Xhat[selected_cities_date,1],Xhat[selected_cities_date,2], labels=names(selected_cities_date), cex=.7)
legend("topright", pch=16, cex=1, col=continent_colours, legend=names(continent_colours))
box()
