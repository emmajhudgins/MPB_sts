### data generation to test MPB interprovincial cooperation model
library(sf)
library(pdist)
library(ggplot2)
library(terra)
data_mpb<-read_sf('./Data/Spatial/grid_sum_5k.shp')
data_mpb<-data_mpb[,-c(60:62)]
spread_area_obs<-c(252,157,71,133,141,252,202,144,55,92) # observed spread area in the presence of management from Goodsman et al.
host_data<-read.table('./Data/Spatial/pine_density.txt')
host_coords<-read.table('./Data/Spatial/xycoords.txt')
host_data<-cbind(host_coords[,2:3], host_data[,2])

host_data2<-vect('./Data/Spatial/full_host.shp')
host_data2<-as.data.frame(cbind(host_data2$`_mean`/100, st_coordinates(st_centroid(st_as_sf(host_data2)))))
colnames(host_data)<-colnames(host_data2[c(2,3,1)])
host_data2<-rbind(host_data, host_data2)
host_data2<-unique.data.frame(host_data2)


empty<-which(rowSums(as.data.frame(data_mpb[,c(7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67)])[,1:21])==0)
full<-which(rowSums(as.data.frame(data_mpb[,c(7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67)])[,1:21])!=0)
grid_old<-data.frame(st_coordinates(st_centroid(data_mpb)$geometry)) 
dists<-as.matrix(pdist(grid_old[empty,], grid_old[full,]))
distmin<-apply(dists,1,min) #minimum distance between invaded and uninvaded cells
data_mpb<-data_mpb[which(distmin<300000),] # remove faraway cells
grid<-as.data.frame(cbind(grid_old[which(distmin<300000),],data_mpb))
colnames(grid)[1:2]<-c("x_coord", "y_coord")
dist<-(dist(grid[,1:2]))
dist<-dist/5000 # scale to grid size



w_det=0.008 # min density for management
w_min=0.004
thresh=1
tmax=16
w_allee<-0
e<-0.99
b<-2
w_spr<-0.02
w_max<-rep(1, nrow(grid))
mpb_time<-matrix(0, nrow(grid), 10)
mpb_time[,1]<-grid$`_06fcount_`+grid$`_06redsum_` # initial infestation

for (t in 2:10)
{
  mpb_time[,t]<-mpb_time[,t-1]+grid[, (-3*t)+41]+grid[,(3*t)+39]
}
front_time<-matrix(0, nrow(grid), 10)

front_time[,1]<-grid$`_06fsum_su`
for (t in 2:10)
{
  front_time[,t]<-mpb_time[,t-1]+grid[, (-3*t)+41]
}

host_data<-st_as_sf(host_data2, coords = c("X","Y"))
#convert host_data to an sf grid
st_crs(host_data)<-st_crs(data_mpb)
#take value of host_data at grid coordinates
grid$host<-host_data$V1[st_nearest_feature(data_mpb, host_data)]
grid$host<-scale(grid$host, center = T)
heur_fit<<-matrix(0,nrow(grid2),6)
tmax=17
set.seed(1)
coin<-runif(nrow(grid))

mpb_fit<-function(par)
{
  neighbourhood<-8.1
    D<-exp(-1*par[1]*dist)
    D[which(dist>neighbourhood)]<-0
    D<-as.numeric(D>coin)
    D<-lapply(1:nrow(grid), function(i){which(D[((i-1)*nrow(grid)+1):(i*nrow(grid))]==1)})
    w_spr<<-0.02
    w_det<-w_min*2
    w_treat<-1
    y<<-rep(1.45, nrow(grid))
  
  w_spr <<- rep(0.02, nrow(grid))+par[2]*grid$host #negative impact on spread when higher


  y <<-rep(1.45, nrow(grid))+par[3]*grid$host # positive impact on spread when higher

  
  w_time<<-w_prime<-w2prime<-w3prime<-w4prime<-q<-matrix(0,nrow(grid),tmax)
w_time[which(grid$`_06fsum_su`>0),1]<<-w_spr[which(grid$`_06fsum_su`>0)]
w_time[which(grid$`_06redcoun`>0& front_time[,1]==0),1]<<-w_spr[which(grid$`_06redcoun`>0& front_time[,1]==0)]*1.45^3
w_time[which(grid$`_07fsum_su`>0),1]<<-w_det/y[which(grid$`_07fsum_su`>0)]
w_time[which(grid$`_08fsum_su`>0),1]<<-w_det/y[which(grid$`_08fsum_su`>0)]^2

for (time in 1:(tmax))
{

   w_prime[,time]<-w_time[, time]*y*(1-q[,time]*e)
   w_prime[which(w_prime[,time]>w_max),time]<-w_max[which(w_prime[,time]>w_max)]
  w2prime[,time]<-w_prime[,time]
  w2prime[which((w_prime[,time]-w_max)>0),time]<-w_max[which((w_prime[,time]-w_max)>0)]
  w3prime[,time]<-w2prime[,time]
   # only allow emigration from cells with w2prime>w_spr within D==1
  for (i in which(w2prime[,time]<w_min)){
        if(length(D[[i]][which(w2prime[D[[i]],time]>=w_spr[D[[i]]]& w2prime[D[[i]],time]<w_treat)])>0)
        {w3prime[i,time]<-w_min}
  }
  
 w4prime[,time]<-w3prime[,time]
 w4prime[which(w4prime[,time]>=w_max),time]<-w_max[which(w4prime[,time]>=w_max)]
if (time!=tmax){
 w_time[,time+1]<<-w4prime[,time]}
}

accuracy<-0
for (i in 7:16){
accuracy[i-5]<-spread_area_obs[i-6]*(1/0.6)-(length(which(w_time[,i]>w_det))-length(which(w_time[,i-1]>w_det)))
}
spread=0
for (i in 7:16){
  spread[i-6]<-sum(w_time[,i]>w_det)-sum(w_time[,i-1]>w_det)
}
if (length(which(spread==0))>4){
  return(1e+11)
}
print(colSums(w_time>w_det))
return(sum(accuracy^2))
}

m<-optim(par=c( 0.9, -0.04332864,  0.33481716), fn=mpb_fit, control=list(trace=100, parscale=c(0.1,0.01,0.1)))
xx<-m$value
yy<-1e+11
while(xx!=yy)
{
  yy<-xx
  m<-optim(par=c(m$par), fn=mpb_fit, control=list(trace=100))
  xx<-m$value

}
saveRDS(m, file="./Results/mpb_unmanage_pars_clean.RDS")

m<-readRDS('./Results/mpb_unmanage_pars_clean.RDS')
mpb_fit(m$par)

#RMSE
sqrt(m$value/10)

write.csv(w_time[,17], file="./Input/real_modelled_startingpoint_newthresh.csv", row.names=F)
