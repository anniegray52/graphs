library(tidyverse)
library(corrr)
#library(corrplot)
library(data.table)
library(RSpectra)
library(igraph)
library(plotly)
library(Matrix)
#library(ks)

figfont = 20 #font size for plotting


# Functions:
doublecentre <- function(B){
  colmeans = matrix(rep(colMeans(B),each=nrow(B)),nrow=nrow(B),ncol=ncol(B), byrow=FALSE)
  rowmeans = matrix(rep(rowMeans(B),each=nrow(B)),nrow=nrow(B),ncol=ncol(B), byrow=TRUE)
  grandmean = matrix(mean(B), nrow=nrow(B), ncol=nrow(B))
  B - colmeans - rowmeans + grandmean
}
sym <- function(s){s[lower.tri(s)] = t(s)[lower.tri(s)]; s}

##data loading and wrangling
############################
#the temperature data is available from:
#https://www.kaggle.com/berkeleyearth/climate-change-earth-surface-temperature-data
temps = read_csv("GlobalLandTemperaturesByCity.csv")

just_temps = select(temps, -c(Latitude, Longitude, AverageTemperatureUncertainty))
##remove any records if more than one temp recorded in one city (NB can be cities in different countries with same names!!) per day
just_temps <- just_temps %>% distinct(.,dt,City,Country, .keep_all=TRUE)
just_temps <- just_temps %>% unite(CityCountry, City, Country, sep = ",", remove = TRUE)
just_temps_by_city <- just_temps %>% pivot_wider(., names_from = CityCountry, values_from = AverageTemperature )

Y = just_temps_by_city %>% select(-dt) %>% drop_na()
Y <- Y[ , order(names(Y))]

#some data wrangling to get the latitude and longitude
CityCountry = names(Y)
geog = data.frame(CityCountry)
geog['lat'] <- NA
geog['lng'] <- NA

#load a file that has correct lat/long for a variety of places from
#https://simplemaps.com/resources/free-country-cities
cities = read_csv("worldcities.csv")
cities<- cities%>% unite(CityCountry, city, country, sep = ",", remove = TRUE)
cities <- cities[order(cities$CityCountry),]

#find the lat/long for the cities that match
city_idx = match(geog$CityCountry,cities$CityCountry)
geog=geog[!is.na(city_idx),]  #drop those from geog who are not listed in cities
cities=cities[city_idx[!is.na(city_idx)],] #...and then also drop those ones from cities
geog$lat <- cities$lat
geog$lng <- cities$lng

#drop the data for which we don't know city's lat/long
Y=Y[,!is.na(city_idx)]
#....and finally compute the A matrix
A <-cor(Y, method = "pearson")

#colour by latitude
fig1=plot_ly(data=geog,x=~lng, y=~lat, marker = list(color = ~geog$lat, 
                                                     colorscale = "Portland",size=3, showscale = FALSE,showlegend =FALSE)) 
#...or by longitude
fig1=plot_ly(data=geog,x=~lng, y=~lat, marker = list(color = ~geog$lng, 
                                                     colorscale = "Portland",size=3, showscale = FALSE,showlegend =FALSE)) 

fig1       

fig1 = fig1 %>% layout(
  #title = list(text="<b>Latitude and longitude of temperature recordings</b>",y=0.99),
  #titlefont=list(size=22),
  xaxis = list(  title="",tickfont=list(size=figfont )),
  yaxis = list(  title="",tickfont=list(size=figfont ))
  )
fig1

## now get on with the analysis
##############################

## do the SVD of A
n = nrow(A)
p = 3; #choose rank for spectral embedding
s=svds(A,p)
X = s$u %*% diag(sqrt(s$d))
DistX = as.matrix(dist(X))

## dataframe for use in manifold plot: 
df <- as.data.frame(X)

## plot the spectral embedding
# plotting when p=2
fig2 = plot_ly(x = ~df$V1, y = ~df$V2, marker = list(color = ~geog$lat, colorscale = "Portland",size=4, showscale = FALSE,showlegend =FALSE)) 
fig2
fig2 = fig2 %>% layout(
  xaxis = list(range = c(-1, 1), title="", tickfont=list(size=figfont )),
  yaxis = list(range = c(-1, 1), title="",tickfont=list(size=figfont ))
  )
fig2

# plotting when p=3
fig2 = plot_ly(x = ~df$V1, y = ~df$V2,z =~df$V3,  marker = list(color = ~geog$lat, colorscale = "Portland",size=4, showscale = FALSE,showlegend =FALSE)) 
fig2 = fig2 %>% layout(
  scene=list(
  xaxis = list(range = c(-1, 1), title="", tickfont=list(size=figfont-6 )),
  yaxis = list(range = c(-1, 1), title="",tickfont=list(size=figfont-6 )),
  zaxis = list(range = c(-1, 1), title="",tickfont=list(size=figfont-6 ))
  )
)
fig2

## k-nearest neighbour graph
library(cccd)
k = 200 
knn = nng(X, k=k, mutual=TRUE)
components(knn)$no
adj = as_adjacency_matrix(knn)
newdist = as.matrix(DistX)*adj
kgraph = graph_from_adjacency_matrix(newdist, weighted=TRUE, mode="min")
print(paste0("Number of graph components: ", components(kgraph)$no))

## epsilon graph
a = 0.5
epsilon = quantile(DistX, a)
newdist=DistX
newdist[DistX > epsilon] = 0
epsilongraph = graph_from_adjacency_matrix(newdist, weighted=TRUE,mode="min")
##want eps to be as low as possible for number of components to be 1
print(paste0("Number of graph components: ", components(epsilongraph)$no))

##use kgraph or epsilongraph
ManifoldDist=distances(kgraph)

##CMDS if d=1
d<-1 
#apply it to the manifold distances or ambient distances
InnerProductMatrix = doublecentre(ManifoldDist^2)
#InnerProductMatrix = doublecentre(DistX^2)
s = svds(InnerProductMatrix,k=d)
Zhat = s$u*sqrt(s$d[1:d])

#if d = 1
##centre and rescale Zhat for plotting
Zhat <- Zhat-min(Zhat)
Zhat <- Zhat/max(Zhat)
Zhat <- Zhat*2 -1

fig3=plot_ly(data= cbind(as.data.frame(Zhat),Latitude=geog$lat,Longitude=geog$lng),
              y=~Longitude, x=~Zhat[,1], marker = list(color = ~geog$lng, colorscale = "Portland",
                                                      size=4, showscale = FALSE)) 
fig3 <- fig3 %>% layout(
  xaxis = list(  title="",tickfont=list(size=figfont )),
  yaxis = list(  title="",tickfont=list(size=figfont ))
)

fig3

#make the plot for the main part of the paper...
figs = subplot(fig1,fig2,fig3, widths=c(1/2,1/4,1/4),titleX = FALSE ,titleY = FALSE, margin=0.035) %>% layout(annotations = list(
    list(x = -0.05 , y = 1.1, text = "<b>a)</b>", showarrow = F, xref='paper', yref='paper', font = list(size = figfont+2)),
    list(x = 0.5 , y = 1.1, text = "<b>b)</b>", showarrow = F, xref='paper', yref='paper', font = list(size = figfont+2)),
    list(x = 0.77 , y = 1.1, text = "<b>c)</b>", showarrow = F, xref='paper', yref='paper', font = list(size = figfont+2))))
figs <- figs %>% layout(showlegend = FALSE)
figs

#make the longitude plot for the appendix
figs = subplot(fig2,fig3, widths=c(1/2,1/2),titleX = FALSE ,titleY = FALSE, margin=0.035) %>% layout(annotations = list(
  list(x = -0.05 , y = 1.1, text = "<b>a)</b>", showarrow = F, xref='paper', yref='paper', font = list(size = figfont+2)),
  list(x = 0.5 , y = 1.1, text = "<b>b)</b>", showarrow = F, xref='paper', yref='paper', font = list(size = figfont+2))))
figs <- figs %>% layout(showlegend = FALSE)
figs

###CMDS if d >1
d<-2
InnerProductMatrix = doublecentre(ManifoldDist^2)
#InnerProductMatrix = doublecentre(DistX^2)
s = svds(InnerProductMatrix,k=d)
Zhat = s$u[,1:d]%*% diag(sqrt(s$d[1:d]))

## centre and rescale Zhat for plotting
Zhat[,1] <- Zhat[,1]-min(Zhat[,1])
Zhat[,1] <- Zhat[,1]/max(Zhat[,1])
Zhat[,1] <- Zhat[,1]*2 -1
Zhat[,2] <- Zhat[,2]-min(Zhat[,2])
Zhat[,2] <- Zhat[,2]/max(Zhat[,2])
Zhat[,2] <- Zhat[,2]*2 -1

## plotting for d>1: 
fig5 = plot_ly(data= cbind(as.data.frame(Zhat),Latitude=geog$lat,Longitude=geog$lng),y=~Longitude, x=~Zhat[,1], marker = list(color = ~geog$lng, colorscale = "Portland",size=3, showscale = FALSE)) 
fig4 = plot_ly(data= cbind(as.data.frame(Zhat),Latitude=geog$lat,Longitude=geog$lng),y=~Latitude, x=~Zhat[,2], marker = list(color = ~geog$lat, colorscale = "Portland",size=3, showscale = FALSE)) 
fig4 = fig4 %>% layout(
  xaxis = list(  title="",tickfont=list(size=figfont )),
  yaxis = list(  title="",tickfont=list(size=figfont ))
)
fig5 = fig5 %>% layout(
  xaxis = list(  title="",tickfont=list(size=figfont )),
  yaxis = list(  title="",tickfont=list(size=figfont ))
)

figs = subplot(fig4,fig5, widths=c(1/2,1/2),titleX = FALSE ,titleY = FALSE, margin=0.035) %>% layout(annotations = list(
  list(x = -0.05 , y = 1.1, text = "<b>a)</b>", showarrow = F, xref='paper', yref='paper', font = list(size = figfont+2)),
  list(x = 0.5 , y = 1.1, text = "<b>b)</b>", showarrow = F, xref='paper', yref='paper', font = list(size = figfont+2))))
figs <- figs %>% layout(showlegend = FALSE)
figs



