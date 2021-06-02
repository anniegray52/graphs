###This code is not intended to be to be run as is, because the relevant datasets need to be downloaded (links provided below), computations take a few hours, and print output to disk. The purpose of distributing this code is to keep a record of calculations.


library(data.table)
library(RSpectra)
library(igraph)
library(Matrix)
library(MASS)
library(irlba)
library(umap)
library(Rtsne)
library(lpSolve)
library(spatstat)


WassersteinD <- function(CD, N, M) {
    ##With thanks to Louis Gammelgaard Jensen
    f.obj = c(rep(1/N, N), -rep(1/M, M))    
    pairs = t(combn(1:(M+N), 2))
    
    MA = matrix(0, nrow = (M+N)*(M+N-1)/2, ncol = M+N)
    V = numeric(nrow(pairs))

    for(r in 1:nrow(pairs)) {
        rr = pairs[r,]
        MA[r,rr[1]] = 1
        MA[r,rr[2]] = -1

        V[r] = CD[rr[1],rr[2]]
    }

    VF = c(V, -V)
    MF = rbind(MA,MA)
    
    cdir = c(rep("<=", nrow(pairs)), rep(">=", nrow(pairs)))
    
    sol = lp(direction = "max", f.obj, MF, cdir, VF)
    return(sol$objval)
}




doublecentre <- function(B){
  colmeans = matrix(rep(colMeans(B),each=nrow(B)),nrow=nrow(B),ncol=ncol(B), byrow=FALSE)
  rowmeans = matrix(rep(rowMeans(B),each=nrow(B)),nrow=nrow(B),ncol=ncol(B), byrow=TRUE)
  grandmean = matrix(mean(B), nrow=nrow(B), ncol=nrow(B))
  B - colmeans - rowmeans + grandmean
}
sym <- function(s){s[lower.tri(s)] = t(s)[lower.tri(s)]; s}

##How the pandemic spread
###https://zenodo.org/record/4737390#.YJemE2ZKhBw 
###https://ourairports.com/data/
rawairports = fread("airports.csv")
airports = rawairports[,c("iso_country", "gps_code","continent")]
airports = airports[gps_code!="",]
airports=airports[!duplicated(gps_code),]
setkey(airports, gps_code)



files=c("flightlist_20191101_20191130.csv.gz","flightlist_20191201_20191231.csv.gz","flightlist_20200101_20200131.csv.gz","flightlist_20200201_20200229.csv.gz","flightlist_20200301_20200331.csv.gz","flightlist_20200401_20200430.csv.gz")

D=10
qepsilon = 0.05
remove_continents=c()


Xhats = list()
for (i in 1:length(files)){
    print(i)
    filename = files[i]
    date=strsplit(filename,"_")[[1]][2]
    rawdata=fread(filename)
    data=rawdata[(rawdata$origin!="") & (rawdata$destination!=""),]
    with_origin_info = merge(data, airports, by.x = "origin", by.y = "gps_code"); setnames(with_origin_info,"iso_country","origin_country");setnames(with_origin_info,"continent","origin_continent")
    datax = merge(with_origin_info, airports, by.x = "destination", by.y = "gps_code"); setnames(datax,"iso_country","destination_country"); setnames(datax,"continent","destination_continent")
    dataC = datax[,.(count = .N, origin = unique(origin), destination = unique(destination)), by = paste(pmin(origin, destination), pmax(origin, destination))]
    dataC=dataC[!(airports[origin,continent]%in% remove_continents) & !(airports[destination,continent]%in% remove_continents),]
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

    continents=airports[rownames(A), continent]

    DistX = as.matrix(dist(X))
    epsilon = quantile(DistX, qepsilon)
    newdist=DistX
    newdist[DistX > epsilon] = 0

    epsilongraph = graph_from_adjacency_matrix(newdist, weighted=TRUE,mode="min")

    print(paste0("Number of graph components: ", components(epsilongraph)$no))
    components=clusters(epsilongraph)
    largest=which.max(components$csize)
    id2names[which(components(epsilongraph)$membership!=largest)]

    egraph=induced.subgraph(epsilongraph, V(epsilongraph)[components$membership==largest])
    ##this could take a while
    ManifoldDist=distances(egraph)    
    retained = as.numeric(names(V(egraph)))
    rownames(ManifoldDist) = id2names[retained]; colnames(ManifoldDist) = id2names[retained]
    Xdenoised = X[retained,]; rownames(Xdenoised)= id2names[retained]
    {if (length(remove_continents)==0){savefile=paste0("ManifoldDist_",date,".Rout")}
    else savefile = paste0("ManifoldDist_",date,"_", paste(remove_continents, collapse="-"),".Rout")}
##    save(ManifoldDist, Xdenoised, epsilon, file=savefile)
##    load(savefile)

    ## ##printing out graph for node2vec
    ## rows=which((names(E1) %in% rownames(ManifoldDist)) & (names(E2) %in% rownames(ManifoldDist)))
    ## edges=cbind(names(E1),names(E2))[rows,]
    ## nd=!duplicated(apply(edges, 1, function(v) paste(sort(v), collapse=" ")))    
    ## edges=edges[nd,]
    ## selfloops=apply(edges, 1, function(v) v[1] == v[2])
    ## edges=edges[!selfloops,]
    ## G = make_graph(c(t(edges)), directed=FALSE)
    ## components=clusters(G)
    ## largest = which.max(components$csize)
    ## subgraph = apply(edges, 1, function(v) (components$membership[v[1]]==largest) & (components$membership[v[2]]==largest))    
    ## edgessub = edges[subgraph,]
    ## ## ##check
    ## ## Gsub = make_graph(c(t(edgessub)))
    ## ## clusters(Gsub)$csize    
    ## write.table(edgessub, file=paste0(date, "_graph.csv"), row.names=FALSE, col.names=FALSE,sep=",")
        
    InnerProductMatrix = doublecentre(ManifoldDist^2)
    s = svds(InnerProductMatrix,k=2)
    Xhat = s$u[,1:2]%*% diag(sqrt(s$d[1:2]))
    rownames(Xhat) = rownames(ManifoldDist)
    Xhats[[length(Xhats)+1]]=Xhat
    selected_cities_date = selected_cities[which(selected_cities %in% rownames(Xhat))]
    {if (length(remove_continents)==0){pdffile=paste0("city_embedding_",date,"isomap.pdf")}
    else pdffile=paste0("city_embedding_",date,"_", remove_continents,"_isomap.pdf")}

}


## ##test the graphs are good:
## edges = fread("20191201_graph.csv")
## g = make_graph(c(t(edges)), directed=FALSE)
## (clusters(g))$csize

##save(Xhats,file="Xhats.Rout")
##load("Xhats.Rout")


##library(ggsci)
continent_colours=c("#FDAF91FF" ,"#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#00468BFF")
names(continent_colours) = c("NA","EU","AS","SA","AF","OC")
selected_cities=c("Los Angeles"="KLAX", "Chicago"="KORD", "Dallas"="KDFW", "New York"="KJFK","Dubai"="OMDB", "Paris"="LFPG", "London"="EGLL", "Rome"="LIRF", "Madrid"="LEMD", "Lisbon"="LPPT", "Athens"="LGAV", "Tokyo"="RJTT", "Sao Paolo"="SBGR","Kempton Park"="FAOR","Mexico City"="MMMX","Auckland"="NZAA", "Buenos Aires"="SAEZ")

## picked=3
## filename = files[picked]
## date=strsplit(filename,"_")[[1]][2]
## savefile=paste0("ManifoldDist_",date,".Rout")
## load(savefile)
## load("Xhats.Rout")

pdf("Jan_isomap.pdf", height=6, width=12)
par(mfrow=c(1,2))
par(mar=c(.5, .5, 1.5, .5))
continents=airports[rownames(Xdenoised), continent]
plot(Xdenoised[,1:2], col=continent_colours[continents], cex=.3,pch=16, xlab="", ylab="", axes=FALSE, main="a) Spectral embedding, first two dimensions (Jan 20)")
selected_cities_date = selected_cities[which(selected_cities %in% rownames(Xhat))]
text(Xdenoised[selected_cities_date,1],Xdenoised[selected_cities_date,2], labels=names(selected_cities_date), cex=.7)
box()
Xhat = Xhats[[picked]]
continents=airports[rownames(Xhat), continent]
plot(Xhat[,1:2], col=continent_colours[continents], cex=.3,pch=16, xlab="", ylab="", axes=FALSE, main="b) Spectral embedding, followed by Isomap (Jan 20)")
selected_cities_date = selected_cities[which(selected_cities %in% rownames(Xhat))]
text(Xhat[selected_cities_date,1],Xhat[selected_cities_date,2], labels=names(selected_cities_date), cex=.7)
legend("topright", pch=16, cex=1, col=continent_colours, legend=names(continent_colours))
box()
dev.off()



## ###UMAP, t-SNE, node2vec2d, node2vec10d computations
## picked=3
## filename = files[picked]
## date=strsplit(filename,"_")[[1]][2]
## savefile=paste0("ManifoldDist_",date,".Rout")
## load(savefile)
## u = umap(Xdenoised)
## Xhatu = u$layout
## rownames(Xhatu) = rownames(Xdenoised)

## tsne_out <- Rtsne(Xdenoised, check_duplicates=FALSE)
## Xhatt = tsne_out$Y
## rownames(Xhatt) = rownames(Xdenoised)

## node2vec2=fread(paste0(date,"_n2v2d.csv"))
## graph_data = as.matrix(fread(paste0(date,"_graph.csv"),header=FALSE))
## names = unique(c(t(graph_data)))
## Xhatn2v2=data.matrix(node2vec2)
## rownames(Xhatn2v2) = names

## savefile=paste0("ManifoldDist_",date,"_node2vec.Rout")
## load(savefile)
## InnerProductMatrix = doublecentre(ManifoldDist^2)
## s = svds(InnerProductMatrix,k=2) ####six continents    
## Xhatn2viso = s$u[,1:2]%*% diag(sqrt(s$d[1:2]))
## rownames(Xhatn2viso) = rownames(ManifoldDist)


## pdf("Jan_comparison_etc.pdf", width=12, height=12)
## par(mfrow=c(2,2))
## par(mar=c(.5, .5, 1.5, .5))

## selected_cities=c("Los Angeles"="KLAX", "Chicago"="KORD", "Dallas"="KDFW", "New York"="KJFK","Dubai"="OMDB", "London"="EGLL", "Tokyo"="RJTT", "Sao Paolo"="SBGR")

## continents=airports[rownames(Xhatu), continent]
## xlim = c(quantile(Xhatu[,1], .01), quantile(Xhatu[,1], .99))
## ylim = c(quantile(Xhatu[,2], .01), quantile(Xhatu[,2], .99)) 
## plot(Xhatu, col=continent_colours[continents], cex=.3,pch=16, xlab="", ylab="", axes=FALSE, main="a) Spectral embedding, followed by UMAP", xlim=xlim,ylim=ylim)
## selected_cities_date = selected_cities[which(selected_cities %in% rownames(Xhatu))]
## text(Xhatu[selected_cities_date,1],Xhatu[selected_cities_date,2], labels=names(selected_cities_date), cex=1.2)
## box()

## selected_cities=c("Los Angeles"="KLAX", "Chicago"="KORD", "Dallas"="KDFW", "New York"="KJFK","Dubai"="OMDB", "London"="EGLL", "Tokyo"="RJTT", "Sao Paolo"="SBGR")

## continents=airports[rownames(Xhatt), continent]
## xlim = c(quantile(Xhatt[,1], .01), quantile(Xhatt[,1], .99))
## ylim = c(quantile(Xhatt[,2], .01), quantile(Xhatt[,2], .99)) 
## plot(Xhatt, col=continent_colours[continents], cex=.3,pch=16, xlab="", ylab="", axes=FALSE, main="a) Spectral embedding, followed by t-SNE", xlim=xlim,ylim=ylim)
## selected_cities_date = selected_cities[which(selected_cities %in% rownames(Xhatt))]
## text(Xhatt[selected_cities_date,1],Xhatt[selected_cities_date,2], labels=names(selected_cities_date), cex=1.2)
## box()

## selected_cities=c("Los Angeles"="KLAX", "Chicago"="KORD", "Dallas"="KDFW", "Dubai"="OMDB", "London"="EGLL")
## continents=airports[rownames(Xhatn2v2), continent]
## plot(Xhatn2v2, col=continent_colours[continents], cex=.3,pch=16, xlab="", ylab="", axes=FALSE, main="c) node2vec (2 dimensions)")
## selected_cities_date = selected_cities[which(selected_cities %in% rownames(Xhatn2v2))]
## text(Xhatn2v2[selected_cities_date,1],Xhatn2v2[selected_cities_date,2], labels=names(selected_cities_date), cex=1.2)
## box()

## selected_cities=c("Los Angeles"="KLAX", "Chicago"="KORD", "Dallas"="KDFW", "New York"="KJFK","Dubai"="OMDB", "London"="EGLL", "Lisbon"="LPPT", "Athens"="LGAV", "Tokyo"="RJTT", "Sao Paolo"="SBGR","Kempton Park"="FAOR","Mexico City"="MMMX","Auckland"="NZAA", "Buenos Aires"="SAEZ")
## continents=airports[rownames(Xhatn2viso), continent]
## plot(Xhatn2viso, col=continent_colours[continents], cex=.3,pch=16, xlab="", ylab="", axes=FALSE, main="d) node2vec (10 dimensions), followed by isomap")
## selected_cities_date = selected_cities[which(selected_cities %in% rownames(Xhatn2viso))]
## text(Xhatn2viso[selected_cities_date,1],Xhatn2viso[selected_cities_date,2], labels=names(selected_cities_date), cex=1.2)
## box()
## dev.off()



cities = sapply(1:length(Xhats), function(i) rownames(Xhats[[i]]))
shared_cities=Reduce(intersect, cities)

Procrustes = function(Xhatf, Zf, shared){
    Xhat = Xhatf[shared,]; Z = Zf[shared,]
    s=svd(t(Xhat)%*%Z)
    c=sum(s$d)/sum(diag(Xhat%*%t(Xhat)))
    T = s$u %*% t(s$v)
    meanXhat=as.matrix(c(mean(Xhat[,1]),mean(Xhat[,2])))
    meanZ = as.matrix(c(mean(Z[,1]),mean(Z[,2])))
    b = meanZ - T %*% meanXhat
    Y = Xhatf%*%T
    Y = c*(Xhatf%*%T)
    cbind(Y[,1] + b[1,], Y[,2] + b[2,])
}

main=c("Nov 19","Dec 19","Jan 20", "Feb 20", "Mar 20", "Apr 20")
pdf("sixmonth.pdf", width=12, height=8)
par(mfrow=c(2,3))
par(mar=c(.5, .5, 1.5, .5))
par(oma=rep(.5,4))
for (i in 1:length(Xhats)){
    Xhat = Xhats[[i]]
    continents=airports[rownames(Xhat), continent]
    plot(Procrustes(Xhat[,1:2],Xhats[[3]], shared_cities), col=continent_colours[continents], cex=.4,pch=16, xlab="", ylab="", main=main[i], axes=FALSE)
    if (i == 1) mtext("c)", side=3, at=-3.5, line=.5, font=2, cex=1.3) #legend("topright", pch=16, cex=1, col=continent_colours, legend=names(continent_colours))
    box()
}
dev.off()

## main=c("Nov 19","Dec 19","Jan 20", "Feb 20", "Mar 20", "Apr 20")
## pdf("threemonth.pdf", width=12, height=4)
## par(mfrow=c(1,3))
## par(mar=c(.5, .5, 1.5, .5))
## par(oma=rep(.5,4))
## for (i in 4:length(Xhats)){
##     Xhat = Xhats[[i]]
##     continents=airports[rownames(Xhat), continent]
##     plot(Procrustes(Xhat[,1:2],Xhats[[3]], shared_cities), col=continent_colours[continents], cex=.4,pch=16, xlab="", ylab="", main=main[i], cex.main=1.7, axes=FALSE)
##     if (i == 4) mtext("c)", side=3, at=-3.5, line=.5, font=2, cex=1.2) #legend("topright", pch=16, cex=1, col=continent_colours, legend=names(continent_colours))
##     box()
## }
## dev.off()




## ##All continents from the rest of the world
## load("Xhats.Rout")
## resultsbycontinent=list()
## rep=100
## ns=100
## for (continent in names(continent_colours)){
##     print(continent)
##     results=c()
##     for (i in 1:length(Xhats)){
##         ##Xhat = Xhats[[i]]
##         print(i)
##         filename = files[i]
##         date=strsplit(filename,"_")[[1]][2]
##         savefile=paste0("ManifoldDist_",date,".Rout")
##         load(savefile)
##         continents=airports[rownames(ManifoldDist), continent]
##         res=replicate(rep, {
##             names1=rownames(ManifoldDist[continents==continent,]); names1=sample(names1,ns, replace=TRUE)
##             names2=rownames(ManifoldDist[!continents==continent,]); names2=sample(names2,ns, replace=TRUE)
##             WassersteinD(ManifoldDist[c(names1,names2),c(names1,names2)],ns,ns)
##         })
##         results=rbind(results, c(mean(res), sd(res)/sqrt(length(res))))
##     }
##     resultsbycontinent[[length(resultsbycontinent)+1]] = results    
## }
## ##save(resultsbycontinent, file="resultsbycontinent.Rout")


## pdf("Wasserstein_all.pdf", width=6, height=6)
## labels=c("Nov 19","Dec 19","Jan 20", "Feb 20", "Mar 20", "Apr 20")
## load("resultsbycontinent.Rout")
## plot(c(), axes=FALSE, xlab="", ylab="Distance", xlim=c(.5, nrow(resultsbycontinent[[1]])), ylim=c(1.5, 3.5))
## for (j in 1:length(resultsbycontinent)){
##     results=resultsbycontinent[[j]]
##     lines(1:nrow(results), results[,1], type="l", col=continent_colours[j], lwd=2)
##     for (i in 1:nrow(results)){
##         lines(c(i,i),c(results[i,1]+2*results[i,2],results[i,1]-2*results[i,2]), col=continent_colours[j])
##     }
## }
## ##points(1:6, results[,1], pch=16)
## axis(1, at = 1:nrow(results)-0.5, labels=labels)
## axis(2, at = seq(1.5, 3.5, length=5), seq(1.5, 3.5, length=5))
## legend("topleft", pch=16, cex=1, col=continent_colours, legend=names(continent_colours))
## box()
## dev.off()





