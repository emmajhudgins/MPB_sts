# ### data generation to test MPB interprovincial cooperation model
library(sf)
library(tidyr)
library(basemaps)
library(pdist)
library(ggplot2)
library(terra)
library(scico)
library(RColorBrewer)
library(ggpubr)
library(viridisLite)
library(canadianmaps)
data_mpb<-read_sf('./Data/Spatial/grid_sum_5k.shp')
data_mpb<-data_mpb[,-c(60:62)]

scenarios<-c("1-Base", "2-High_eff","3-Low_det","4-High_det", "5-Low_min","6-High_min", "7-Two_bud", "8-C_fix", "9-C_both", "10-Low_bud", "11-Low_eff","12-Lowest_bud", "13-Low_bud_det")


scenario<-scenarios[1] #choose from 1 to 13 (using a loop won't print to pdf)
##load individual sensitivity analysis

empty<-which(rowSums(as.data.frame(data_mpb[,c(7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67)])[,1:21])==0)
full<-which(rowSums(as.data.frame(data_mpb[,c(7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67)])[,1:21])!=0)
grid_old<-data.frame(st_coordinates(st_centroid(data_mpb)$geometry))
dists<-as.matrix(pdist(grid_old[empty,], grid_old[full,]))
distmin<-apply(dists,1,min)
data_mpb<-data_mpb[which(distmin<300000),]
grid<-as.data.frame(cbind(grid_old[which(distmin<300000),],data_mpb))
colnames(grid)[1:2]<-c("x_coord", "y_coord")
w_det=0.008 # min density for management
w_min=0.004

results_mid<-read.table(paste0('./Results/',scenario,'/popul_nt_final', '.txt'), sep="", header=F,skip=1, na.strings="")
colnames(results_mid)[1]<-"ID"
#infestation below minimum
results_mid_u0<-results_mid[,1:12]
colnames(results_mid_u0)[2:12]<-c(1:11)
#popden
results_mid_w<-abs(results_mid[,c(1,14:24)])
colnames(results_mid_w)[2:12]<-c(1:11)
#expected popden * 1-q
results_mid_w1q<-results_mid[,c(1,26:36)]
colnames(results_mid_w1q)[2:12]<-c(1:11)
#detectable and treatable
results_mid_v<-results_mid[,c(1,38:48)]
colnames(results_mid_v)[2:12]<-c(1:11)
#capable of spread
results_mid_z<-results_mid[,c(1,50:60)]
colnames(results_mid_z)[2:12]<-c(1:11)
#first infestation
results_mid_x<-results_mid[,c(1,62:72)]
colnames(results_mid_x)[2:12]<-c(1:11)
#management actions
results_mid_q<-results_mid[,c(1,74:84)]
colnames(results_mid_q)[2:12]<-c(1:11)


results_heur<-read.table(paste0('./Results/',scenario,'/popul_nt_1.txt'), sep="", header=F,skip=1, na.strings="")
colnames(results_heur)[1]<-"ID"
#infestation below minimum
results_heur_u0<-results_heur[,1:12]
colnames(results_heur_u0)[2:12]<-c(1:11)
#popden
results_heur_w<-abs(results_heur[,c(1,14:24)])
colnames(results_heur_w)[2:12]<-c(1:11)
#expected popden * 1-q
results_heur_w1q<-results_heur[,c(1,26:36)]
colnames(results_heur_w1q)[2:12]<-c(1:11)
#detectable and treatable
results_heur_v<-results_heur[,c(1,38:48)]
colnames(results_heur_v)[2:12]<-c(1:11)
#capable of spread
results_heur_z<-results_heur[,c(1,50:60)]
colnames(results_heur_z)[2:12]<-c(1:11)
#first infestation
results_heur_x<-results_heur[,c(1,62:72)]
colnames(results_heur_x)[2:12]<-c(1:11)
#management actions
results_heur_q<-results_heur[,c(1,74:84)]
colnames(results_heur_q)[2:12]<-c(1:11)


results_init<-read.table(paste0('./Results/',scenario,'/popul_nt_ini.txt'), sep="", header=F,skip=1, na.strings="")
colnames(results_init)[1]<-"ID"
#infestation below minimum
results_init_u0<-results_init[,1:12]
colnames(results_init_u0)[2:12]<-c(1:11)
#popden
results_init_w<-abs(results_init[,c(1,14:24)])
colnames(results_init_w)[2:12]<-c(1:11)
#expected popden * 1-q
results_init_w1q<-results_init[,c(1,26:36)]
colnames(results_init_w1q)[2:12]<-c(1:11)
#detectable and treatable
results_init_v<-results_init[,c(1,38:48)]
colnames(results_init_v)[2:12]<-c(1:11)
#capable of spread
results_init_z<-results_init[,c(1,50:60)]
colnames(results_init_z)[2:12]<-c(1:11)
#first infestation
results_init_x<-results_init[,c(1,62:72)]
colnames(results_init_x)[2:12]<-c(1:11)
#management actions
results_init_q<-results_init[,c(1,74:84)]
colnames(results_init_q)[2:12]<-c(1:11)

#variability in site management  rates and densities at management
mean(results_mid_w[, 2:12][which(results_mid_q[,2:12]>0, arr.ind=T)], na.rm=T)
sd(rowSums(results_mid_w[, 2:12])[which(rowSums(results_mid_q[,2:12])>0, arr.ind=T)], na.rm=T)
mean(rowSums(results_mid_q[, 2:12])[which(rowSums(results_mid_q[,2:12])>0, arr.ind=T)], na.rm=T)
sd(rowSums(results_mid_q[, 2:12])[which(rowSums(results_mid_q[,2:12])>0, arr.ind=T)], na.rm=T)

mean(results_heur_w[, 2:12][which(results_heur_q[,2:12]>0, arr.ind=T)], na.rm=T)
sd(rowSums(results_heur_w[, 2:12])[which(rowSums(results_mid_q[,2:12])>0, arr.ind=T)], na.rm=T)
mean(rowSums(results_heur_q[, 2:12])[which(rowSums(results_heur_q[,2:12])>0, arr.ind=T)], na.rm=T)
sd(rowSums(results_heur_q[, 2:12])[which(rowSums(results_heur_q[,2:12])>0, arr.ind=T)], na.rm=T)

mean(results_init_w[, 2:12][which(results_init_q[,2:12]>0, arr.ind=T)], na.rm=T)
sd(rowSums(results_init_w[, 2:12])[which(rowSums(results_mid_q[,2:12])>0, arr.ind=T)], na.rm=T)
mean(rowSums(results_init_q[, 2:12])[which(rowSums(results_init_q[,2:12])>0, arr.ind=T)], na.rm=T)
sd(rowSums(results_init_q[, 2:12])[which(rowSums(results_init_q[,2:12])>0, arr.ind=T)], na.rm=T)

#proportion of sites managed at different times after detection

p1<-sum(which(results_heur_w[, 2:12][which(results_heur_q[,2:12]>0, arr.ind=T)]-0.008<=0.004)%in%which(results_heur_w[, 2:12][which(results_heur_q[,2:12]>0, arr.ind=T)]-0.008>=0))
p2<-sum(which(results_heur_w[, 2:12][which(results_heur_q[,2:12]>0, arr.ind=T)]-0.012<=0.004)%in%which(results_heur_w[, 2:12][which(results_heur_q[,2:12]>0, arr.ind=T)]-0.012>=0))
p3<-sum(which(results_heur_w[, 2:12][which(results_heur_q[,2:12]>0, arr.ind=T)]-0.012>=0.004)%in%which(results_heur_w[, 2:12][which(results_heur_q[,2:12]>0, arr.ind=T)]-0.012>=0))
p1/(p1+p2+p3)
p2/(p1+p2+p3)
p3/(p1+p2+p3)

p1<-sum(which(results_init_w[, 2:12][which(results_init_q[,2:12]>0, arr.ind=T)]-0.008<=0.004)%in%which(results_init_w[, 2:12][which(results_init_q[,2:12]>0, arr.ind=T)]-0.008>=0))
p2<-sum(which(results_init_w[, 2:12][which(results_init_q[,2:12]>0, arr.ind=T)]-0.012<=0.004)%in%which(results_init_w[, 2:12][which(results_init_q[,2:12]>0, arr.ind=T)]-0.012>=0))
p3<-sum(which(results_init_w[, 2:12][which(results_init_q[,2:12]>0, arr.ind=T)]-0.012>=0.004)%in%which(results_init_w[, 2:12][which(results_init_q[,2:12]>0, arr.ind=T)]-0.012>=0))
p1/(p1+p2+p3)
p2/(p1+p2+p3)
p3/(p1+p2+p3)


p1<-sum(which(results_mid_w[, 2:12][which(results_mid_q[,2:12]>0, arr.ind=T)]-0.008<=0.004)%in%which(results_mid_w[, 2:12][which(results_mid_q[,2:12]>0, arr.ind=T)]-0.008>=0))
p2<-sum(which(results_mid_w[, 2:12][which(results_mid_q[,2:12]>0, arr.ind=T)]-0.012<=0.004)%in%which(results_mid_w[, 2:12][which(results_mid_q[,2:12]>0, arr.ind=T)]-0.012>=0))
p3<-sum(which(results_mid_w[, 2:12][which(results_mid_q[,2:12]>0, arr.ind=T)]-0.012>=0.004)%in%which(results_mid_w[, 2:12][which(results_mid_q[,2:12]>0, arr.ind=T)]-0.012>=0))
p1/(p1+p2+p3)
p2/(p1+p2+p3)
p3/(p1+p2+p3)


results_spr<-read.table(paste0('./Results/',scenario,'/popul_nt0.txt'), sep="", header=F,skip=1, na.strings="")
colnames(results_spr)[1]<-"ID"
#infestation below minimum
results_spr_u0<-results_spr[,1:12]
colnames(results_spr_u0)[2:12]<-c(1:11)
#popden
results_spr_w<-abs(results_spr[,c(1,14:24)])
colnames(results_spr_w)[2:12]<-c(1:11)
#expected popden * 1-q
results_spr_w1q<-results_spr[,c(1,26:36)]
colnames(results_spr_w1q)[2:12]<-c(1:11)
#detectable and treatable
results_spr_v<-results_spr[,c(1,38:48)]
colnames(results_spr_v)[2:12]<-c(1:11)
#capable of spread
results_spr_z<-results_spr[,c(1,50:60)]
colnames(results_spr_z)[2:12]<-c(1:11)
#first infestation
results_spr_x<-results_spr[,c(1,62:72)]
colnames(results_spr_x)[2:12]<-c(1:11)
#management actions
results_spr_q<-results_spr[,c(1,74:84)]
colnames(results_spr_q)[2:12]<-c(1:11)

#coords
coordinates<-cbind(read.table('./Input/x_coord.txt'), (read.table('./Input/y_coord.txt')[,2]))
st_bbox(w_coordinates[which(results_mid_x[,12]==1),])
st_centroid(st_combine(w_coordinates[which(results_mid_x[,12]==1),]))
q_coordinates<-st_as_sf(cbind(coordinates, results_mid_q), coords=c(2,3))
colnames(q_coordinates)[3:13]<-c(2024:2034)
q_coordinates$sum<-rowSums(st_drop_geometry(q_coordinates[,3:13]))

# mapping setup
set_defaults(map_service = "osm", map_type = "topographic")
#create new column with binned values of sum c(0,2,4,6,8,10)
q_coordinates$sum2<-cut(q_coordinates$sum, breaks=c(-3,0,3,6,9,12), labels=c(NA,"1-3","4-6","7-9","9-12"))
st_crs(q_coordinates)<-st_crs(data_mpb)
box<-st_bbox(st_transform(q_coordinates, st_crs(4326)))
box<-c(-119.5,53.2,-108,60.9)
sask_border<-data.frame(x=-110, y=56)
sask_border<-st_as_sf(sask_border, crs=st_crs(4326), coords=c(1,2))
sask_border<-data.frame(st_coordinates(sask_border))[1,1]


#how far East does spread go?
max(coordinates$V2[which((results_spr_v[,12])>0)])/1000
max(coordinates$V2[which((results_init_v[,12])>0)])/1000
max(coordinates$V2[which((results_heur_v[,12])>0)])/1000
max(coordinates$V2[which((results_mid_v[,12])>0)])/1000


pdf("./Plots/Fig1.pdf", width=8, height=6)
par(mar=c(4,4,2,2))
time<-seq(1,16)
y<-rep(0,16)
y[1]<-w_min
for(time in 2:16)
{
  y[time]<-y[time-1]*1.45
}
plot(y=c(y[1:15],0.9),x=c(seq(0,14),14.55), xlim=c(0,15), ylim=c(0,1),type="l", col="black", lwd=2, xlab="Time since invasion", ylab="MPB density relative to carrying capacity", main=NULL,xaxt='n')
points(y=c(0.9,0.004),x=c(14.55,15), ylim=c(0,1),type="l", col="black", lwd=2, xlab="Time since invasion", ylab="MPB density relative to carrying capacity", main=NULL,xaxt='n')

axis(1, at=seq(0,16, by=2), labels=seq(0,16, by=2))
points(y=w_det,x=2.5-0.7, col=brewer.pal(4, "BrBG")[2],pch=19)
points(y=w_min, x=1-1, col=brewer.pal(4, "BrBG")[1], pch=19)
points(y=0.02,x=5.2-1, col=brewer.pal(4, "BrBG")[3], pch=19)
abline(v=15.55-1, col=brewer.pal(4, "BrBG")[4], lty=3, lwd=2)
abline(v=5.2-1, col=brewer.pal(4, "BrBG")[3], lty=3, lwd=2)
abline(v=2.5-0.7, col=brewer.pal(4, "BrBG")[2], lty=3, lwd=2)
abline(v=1-1, col=brewer.pal(4, "BrBG")[1], lty=3, lwd=2)
points(y=0.9, x=15.55-1,col=brewer.pal(4, "BrBG")[4], pch=19)
text(x=15.55-1, y=0.9, labels=expression("w"["max"],),pos=2)
text(x=5.2+0.7-1, y=0.02, labels=expression("w"["spr"],),pos=3)
text(x=2.5+0.7-1, y=w_det, labels=expression("w"["det"]),pos=3)
text(x=1+0.7-1, y=w_min, labels=expression("w"["min"]),pos=3)
dev.off()

pdf(paste0("./Plots/", scenario,"/Fig2.pdf"), width=8, height=6) 
    par(mar=c(4,4,2,2))
    par(mfrow=c(2,2))
    plot(colSums(results_spr_v[,2:12]), type="l", col=scico(4, palette="roma")[1],ylab="Invaded Cells", xlab="Timestep", ylim=c(0,2200), lwd=2)
    points(colSums(results_heur_v[,2:12]), type="l",col=scico(4, palette="roma")[2], lwd=2)
    points(colSums(results_mid_v[,2:12]), type="l",col=scico(4, palette="roma")[4], lwd=2)
    points(colSums(results_init_v[,2:12]), type="l",col=scico(4, palette="roma")[3], lwd=2)
    legend('topleft', fill=scico(4, palette="roma"), legend=c("No Management", "Single Period","Single Period + Spread", "Multiperiod"))
    mtext("a.",side=3, adj=0)
    
hist(rowSums(results_heur_q[which(rowSums(results_heur_q[,2:12])>0),2:12]), xlab="Management frequency",  col=scico(4,palette="roma")[2], prob=T, breaks=c(1:10),  ylim=c(0,0.8), main="") #main="Single Period",
mtext("b.",side=3, adj=0)
hist(rowSums(results_init_q[which(rowSums(results_init_q[,2:12])>0),2:12]), xlab="Management frequency", col=scico(4,palette="roma")[3], prob=T, breaks=c(1:10), ylim=c(0,0.8), main="") #main="Single Period + Spread"
mtext("c.",side=3, adj=0)
hist(rowSums(results_mid_q[which(rowSums(results_mid_q[,2:12])>0),2:12]), xlab="Management frequency",prob=T, col=scico(4,palette="roma")[4], breaks=c(1:10),, ylim=c(0,0.8), main="")# main="Multiperiod"
mtext("d.", side=3, adj=0)

dev.off()

#time of first infestation
first_inf_heur<-apply(results_heur_w[,2:12],1,function(x){min(which(x>=w_det))})
first_inf_heur[which(first_inf_heur==Inf)]<-NA
unique(first_inf_heur)
w_coordinates$first_inf_heur<-first_inf_heur
w_coordinates$first_inf_heur<-first_inf_heur-first_inf
w_coordinates$first_inf_heur[which(is.na(first_inf))]<-200
w_coordinates$first_inf_heur[which(is.na(first_inf_heur))]<-100
w_coordinates$first_inf_heur[which(is.na(first_inf_heur)&is.na(first_inf))]<-NA

w_coordinates$first_inf_heur[which(results_heur_w[,2]>w_det)]<-300
infested_cells<-w_coordinates[which(!is.na(w_coordinates$first_inf)),3:14]
w_coordinates$first_inf_heur2<-cut(w_coordinates$first_inf_heur, breaks=c(-6:4, 101,201,301), labels=c(-5:4, "Never infested in myopic", "Never infested in full", "Initial distribution"))


st_crs(w_coordinates)<-st_crs(data_mpb)
coords<-st_coordinates(st_transform(w_coordinates, st_crs(3857)))
CD_AL<-subset(CD, PRNAME%in%c("Alberta", "Saskatchewan"))
st_bbox(CD_AL)



## Frequency of management figure
pdf(paste0('./Plots/',scenario,'/Fig3.pdf'), width=4, height=7)
g1<-ggplot()+geom_cd(data = CD_AL, colour = "black", size = 0.2, fill=NA)+scale_colour_manual(values=c("NA",rep("NA",5)), na.value=NA)+
  geom_sf(data =st_transform(q_coordinates, st_crs(4326)), aes(colour=sum2,fill = sum2, alpha=sum2>0),pch=22,size=0.6,inherit.aes = FALSE)+
  scale_x_continuous(limits = c(-118.00128,-108.5), expand=c(0,0))+scale_y_continuous(limits = c(53, 58.0006), expand=c(0,0))+
  guides(alpha="none", colour='none')+
  geom_text()+
  annotate("text",label="ALBERTA", colour="darkgrey", x=-117.0, y=57.7, size=2)+
  annotate("text",label="SASKATCHEWAN", colour="darkgrey", x=-109.5, y=56.7, size=2, angle=270)+
    theme(legend.title = element_text(size=10))+
  #guides(fill=guide_legend(title="Management Frequency", position="bottom", direction="horizontal", nrow=1, override.aes = list(size = 5)))+
  #ggtitle("Multiperiod Model")+
  scale_fill_manual(values=c("NA",brewer.pal(4,'YlOrBr')),breaks=as.factor(c(NA,"1-3","4-6","7-9","9-12")), na.value = NA, drop=F)+
  theme_minimal()+theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),
    title=element_text(size=10))+
  geom_vline(colour='darkred', xintercept=sask_border, show.legend=FALSE)+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0, "pt"), #length of tick marks
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(), 
        plot.margin = unit(c(0,0,0,0),"mm"))+labs(tag="c.")

q_coordinates$sum_init<-rowSums(results_init_q[,2:12])
q_coordinates$sum_init2<-cut(q_coordinates$sum_init, breaks=c(-3,0,3,6,9,12), labels=c(NA,"1-3","4-6","7-9","9-12"))

g2<-ggplot()+geom_cd(data = CD_AL, colour = "black", size = 0.2, fill=NA)+scale_colour_manual(values=c("NA",rep("NA",5)), na.value=NA)+
  geom_sf(data =st_transform(q_coordinates, st_crs(4326)), aes(colour=sum_init2,fill = sum_init2, alpha=sum_init2>0),pch=22,size=0.6,inherit.aes = FALSE)+
  scale_x_continuous(limits = c(-118.00128,-108.5), expand=c(0,0))+scale_y_continuous(limits = c(53, 58.0006), expand=c(0,0))+
  guides(alpha="none", colour='none')+
  geom_text()+
  annotate("text",label="ALBERTA", colour="darkgrey", x=-117.0, y=57.7, size=2)+
  annotate("text",label="SASKATCHEWAN", colour="darkgrey", x=-109.5, y=56.7, size=2, angle=270)+
  theme(legend.title = element_text(size=10))+
  guides(fill=guide_legend(title="Management Frequency", position="bottom", direction="horizontal", nrow=1, override.aes = list(size = 5)))+
  #ggtitle("Single Period + Spread Model")+
  scale_fill_manual(values=c("NA",brewer.pal(3,'YlOrBr')),breaks=as.factor(c(NA,"1-3","4-6","7-9","9-12")), na.value = NA, drop=F)+
  theme_minimal()+theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),
    title=element_text(size=10))+
  geom_vline(colour='darkred', xintercept=sask_border, show.legend=FALSE)+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0, "pt"), #length of tick marks
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(), 
        plot.margin = unit(c(0,0,0,0),"mm"))+labs(tag="b.")

q_coordinates$sum_heur<-rowSums(results_heur_q[,2:12])
q_coordinates$sum_heur2<-cut(q_coordinates$sum_heur, breaks=c(-3,0,3,6,9,12), labels=c(NA,"1-3","4-6","7-9","9-12"))

g3<-ggplot()+geom_cd(data = CD_AL, colour = "black", size = 0.2, fill=NA)+scale_colour_manual(values=c("NA",rep("NA",5)), na.value=NA)+
  geom_sf(data =st_transform(q_coordinates, st_crs(4326)), aes(colour=sum_heur2,fill = sum_heur2, alpha=sum_heur2>0),pch=22,size=0.6,inherit.aes = FALSE)+
  scale_x_continuous(limits = c(-118.00128,-108.5), expand=c(0,0))+scale_y_continuous(limits = c(53, 58.0006), expand=c(0,0))+
  guides(alpha="none", colour='none')+
  geom_text()+
  annotate("text",label="ALBERTA", colour="darkgrey", x=-117.0, y=57.7, size=2)+
  annotate("text",label="SASKATCHEWAN", colour="darkgrey", x=-109.5, y=56.7, size=2, angle=270)+
  theme(legend.title = element_text(size=10))+
  guides(fill=guide_legend(title="Management Frequency", position="bottom", direction="horizontal", nrow=1, override.aes = list(size = 5)))+
  #ggtitle("Single Period Model")+
  scale_fill_manual(values=c("NA",brewer.pal(3,'YlOrBr')),breaks=as.factor(c(NA,"1-3","4-6","7-9","9-12")), na.value = NA, drop=F)+
  theme_minimal()+theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),
    title=element_text(size=10))+
  geom_vline(colour='darkred', xintercept=sask_border, show.legend=FALSE)+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0, "pt"), #length of tick marks
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(), 
        plot.margin = unit(c(0,0,0,0),"mm"))+labs(tag="a.")
g4<-ggplot(data =data.frame(sum_heur2=as.factor(c("1-3","4-6","7-9", "10")), x=c(1:4), y=c(1:4)), aes(x=x,y=y,fill = sum_heur2))+
  geom_point(pch=22)+
  theme_minimal()+
  guides(fill=guide_legend(title="Management\n Frequency", position="bottom", direction="horizontal", nrow=2, byrow=T,override.aes = list(size = 5), theme=theme(legend.title=element_text(size=10))))+
  scale_fill_manual(values=c(brewer.pal(4,'YlOrBr')),breaks=as.factor(c("1-3","4-6","7-9","10")), na.value = NA, drop=F)
  
ggarrange(g3, g2, g1, ncol=1, nrow=3, legend.grob=get_legend(g4), legend="bottom", common.legend=T)
dev.off()

pdf(paste0("./Plots/",scenario,"/Fig4.pdf"), width=7, height=5.5)
ggplot()+
  geom_cd(data = CD_AL, colour = "black", size = 0.2, fill=NA)+scale_fill_identity()+
  geom_sf(data =st_transform(w_coordinates, st_crs(4326)), aes(colour = first_inf_heur2),pch=15,size=0.7,inherit.aes = FALSE)+
  #scale_x_continuous(limits = c(-120.00128,-108.5), expand=c(0,0))+scale_y_continuous(limits = c(48.9975, 60.0006), expand=c(0,0))+
  scale_x_continuous(limits = c(-118.00128,-108.5), expand=c(0,0))+scale_y_continuous(limits = c(53, 58.0006), expand=c(0,0))+
  guides(alpha="none")+
  geom_text()+
  annotate("text",label="ALBERTA", colour="darkgrey", x=-117.0, y=57.7, size=5)+
  annotate("text",label="SASKATCHEWAN", colour="darkgrey", x=-109.5, y=56.7, size=5, angle=270)+
  theme(legend.title = element_text(size=10))+
  guides(colour=guide_legend(title=NULL, position="inside", direction="vertical", ncol=1, override.aes = list(size = 3)))+
  ggtitle("Relative time of infestation")+
  scale_colour_manual(values=c(viridis(5), "darkgrey"), labels=c("Never infested in full","Infested earlier in myopic", "Infested at the same time", "Infested earlier in full","Never infested in myopic", "Initial distribution" ), limits=c("Never infested in full",-5,0,4, "Never infested in myopic", "Initial distribution"), na.value = NA)+
  theme_minimal()+theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),
    title=element_text(size=10))+
  geom_vline(colour='darkred', xintercept=sask_border, show.legend=FALSE)+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(), 
        plot.margin = unit(c(0,0,0,0),"mm"))+
  guides(colour=guide_legend(title=NULL, position="right", direction="vertical", ncol=1, override.aes = list(size = 3)))
dev.off()


#example binary spread pattern

dist<-read.table('./Input/spread_2024.txt')
dist<-subset(dist, V3==1 & V1==255)
data2<-st_transform(q_coordinates, st_crs(4326))
data2$x<-st_coordinates(data2)[,1]
data2$y<-st_coordinates(data2)[,2]
data2<-as.data.frame(data2)


# if(scenario=="Base"){
# pdf("./Plots/FigS1.pdf", width=8, height=6)
#   ggplot()+geom_cd(data = CD_AL, colour = "black", size = 0.2, fill=NA)+scale_fill_identity()+
#     geom_sf(data =st_transform(st_as_sf(data2), st_crs(4326)), pch=0, alpha=0.5, cex=0.8,inherit.aes = FALSE)+
#     scale_x_continuous(limits = c(-118.00128,-108.5), expand=c(0,0))+scale_y_continuous(limits = c(53, 58.0006), expand=c(0,0))+
#     guides(alpha="none")+
#   geom_sf(data =st_transform(st_as_sf(data2[which(data2$V1%in%dist$V2),]), st_crs(4326)), pch=15, size=0.8, colour="black")+
#   geom_sf(data =st_transform(st_as_sf(data2[which(data2$V1==255),]), st_crs(4326)), pch=15, size=0.8, colour="yellow")+
#     theme(legend.title = element_text(size=10))+
#   theme_minimal()+theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     title=element_text(size=10))+
#   geom_vline(colour='darkred', xintercept=sask_border, show.legend=FALSE)+
#   theme(axis.line = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks.length = unit(0, "pt"), #length of tick marks
#         panel.background = element_blank(),
#         panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.background = element_blank(),
#         plot.margin = unit(c(0,0,0,0),"mm"))+
#   guides(colour=guide_legend(title=NULL, position="right", direction="vertical", ncol=1, override.aes = list(size = 3)))+
#   geom_text()+
#     annotate("text",label="ALBERTA", colour="darkgrey", x=-117, y=57.7, size=4)+
#     annotate("text",label="SASKATCHEWAN", colour="darkgrey", x=-109, y=56.7, size=4, angle=270)
#   dev.off()
# }

# Spread progression over time for supplement

w_coordinates<-st_as_sf(cbind(coordinates, results_heur_w), coords=c(2,3))
colnames(w_coordinates)[3:13]<-c(2024:2034)
w_long<-pivot_longer(w_coordinates, cols=3:13, names_to="year", values_to="density")
st_crs(w_long)<-st_crs(data_mpb)
w_long$density2<-cut(w_long$density, breaks=c(-0.02,0.004,0.008,0.2,1,18), labels=c("Absent", "Undetectable", "Detectable", "Unmanageable", "Collapsed"))

pdf(paste0("./Plots/", scenario, "/FigS2a.pdf"), width=6, height=7)
ggplot()+
  geom_cd(data = CD_AL, colour = "black", size = 0.2, fill=NA)+scale_fill_identity()+
  geom_sf(data =st_transform(w_long, st_crs(4326)), aes(colour = density2, alpha=density>0),pch=15,size=0.5,inherit.aes = FALSE)+
  scale_x_continuous(limits = c(-120.00128,-108.5), expand=c(0,0))+scale_y_continuous(limits = c(48.9975, 60.0006), expand=c(0,0))+
  guides(alpha="none")+
  geom_text()+
  annotate("text",label="ALBERTA", colour="darkgrey", x=-117.0, y=59.2, size=2)+
  annotate("text",label="SASKATCHEWAN", colour="darkgrey", x=-109.5, y=55.7, size=2, angle=270)+
 guides(alpha="none")+
  facet_wrap(~year, ncol=3)+
  theme(legend.title = element_text(size=10))+
  guides(colour=guide_legend(title="Management Frequency", position="bottom", direction="horizontal", nrow=1, override.aes = list(size = 5)))+
  ggtitle("Single-Period Model")+
  scale_colour_manual(values=c(magma(4,direction=-1)), limits=c("Undetectable", "Detectable", "Unmanageable", "Collapsed"), na.value = NA)+
  theme_minimal()+theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),
    title=element_text(size=10))+
  geom_vline(colour='darkred', xintercept=sask_border, show.legend=FALSE)+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(), 
        plot.margin = unit(c(0,0,0,0),"mm"))+
  guides(colour=guide_legend(title=NULL, position="right", direction="vertical", ncol=1, override.aes = list(size = 3)))
dev.off()


w_coordinates<-st_as_sf(cbind(coordinates, results_init_w), coords=c(2,3))
colnames(w_coordinates)[3:13]<-c(2024:2034)
w_long<-pivot_longer(w_coordinates, cols=3:13, names_to="year", values_to="density")
st_crs(w_long)<-st_crs(data_mpb)
w_long$density2<-cut(w_long$density, breaks=c(-0.02,0.004,0.008,0.2,1,18), labels=c("Absent", "Undetectable", "Detectable", "Unmanageable", "Collapsed"))

pdf(paste0("./Plots/", scenario, "/FigS2b.pdf"), width=6, height=7)
ggplot()+
  geom_cd(data = CD_AL, colour = "black", size = 0.2, fill=NA)+scale_fill_identity()+
  geom_sf(data =st_transform(w_long, st_crs(4326)), aes(colour = density2, alpha=density>0),pch=15,size=0.5,inherit.aes = FALSE)+
  scale_x_continuous(limits = c(-120.00128,-108.5), expand=c(0,0))+scale_y_continuous(limits = c(48.9975, 60.0006), expand=c(0,0))+
  guides(alpha="none")+
  geom_text()+
  annotate("text",label="ALBERTA", colour="darkgrey", x=-117.0, y=59.2, size=2)+
  annotate("text",label="SASKATCHEWAN", colour="darkgrey", x=-109.5, y=55.7, size=2, angle=270)+
  facet_wrap(~year, ncol=3)+
  theme(legend.title = element_text(size=10))+
  guides(colour=guide_legend(title="Management Frequency", position="bottom", direction="horizontal", nrow=1, override.aes = list(size = 5)))+
  ggtitle("Single-Period + Spread Model")+
  scale_colour_manual(values=c(magma(4,direction=-1)), limits=c("Undetectable", "Detectable", "Unmanageable", "Collapsed"), na.value = NA)+
  theme_minimal()+theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),
    title=element_text(size=10))+
  geom_vline(colour='darkred', xintercept=sask_border, show.legend=FALSE)+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(), 
        plot.margin = unit(c(0,0,0,0),"mm"))+
  guides(colour=guide_legend(title=NULL, position="right", direction="vertical", ncol=1, override.aes = list(size = 3)))
dev.off()


w_coordinates<-st_as_sf(cbind(coordinates, results_mid_w), coords=c(2,3))
colnames(w_coordinates)[3:13]<-c(2024:2034)
w_long<-pivot_longer(w_coordinates, cols=3:13, names_to="year", values_to="density")
st_crs(w_long)<-st_crs(data_mpb)
w_long$density2<-cut(w_long$density, breaks=c(-0.02,0.004,0.008,0.2,1,18), labels=c("Absent", "Undetectable", "Detectable", "Unmanageable", "Collapsed"))

pdf(paste0("./Plots/", scenario, "/FigS2c.pdf"), width=6, height=7)
  ggplot()+
  geom_cd(data = CD_AL, colour = "black", size = 0.2, fill=NA)+scale_fill_identity()+
  geom_sf(data =st_transform(w_long, st_crs(4326)), aes(colour = density2, alpha=density>0),pch=15,size=0.5,inherit.aes = FALSE)+
  scale_x_continuous(limits = c(-120.00128,-108.5), expand=c(0,0))+scale_y_continuous(limits = c(48.9975, 60.0006), expand=c(0,0))+
  guides(alpha="none")+
  geom_text()+
  annotate("text",label="ALBERTA", colour="darkgrey", x=-117.0, y=59.2, size=2)+
  annotate("text",label="SASKATCHEWAN", colour="darkgrey", x=-109.5, y=55.7, size=2, angle=270)+
  facet_wrap(~year, ncol=3)+
  theme(legend.title = element_text(size=10))+
  guides(colour=guide_legend(title="Management Frequency", position="bottom", direction="horizontal", nrow=1, override.aes = list(size = 5)))+
  ggtitle("Multiperiod Model")+
  scale_colour_manual(values=c(magma(4,direction=-1)), limits=c("Undetectable", "Detectable", "Unmanageable", "Collapsed"), na.value = NA)+
  theme_minimal()+theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),
    title=element_text(size=10))+
  geom_vline(colour='darkred', xintercept=sask_border, show.legend=FALSE)+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0, "pt"), #length of tick marks
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(), 
        plot.margin = unit(c(0,0,0,0),"mm"))+
  guides(colour=guide_legend(title=NULL, position="right", direction="vertical", ncol=1, override.aes = list(size = 3)))
dev.off()
