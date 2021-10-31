
## trying to plot collagen and laminin gradients
library(splines)

#take each embryo and do a sliding window to get averages for each embryo

WT_embryos<-read.csv('WT_collagen_normalized.csv',stringsAsFactors = F)
KO_embryos<-read.csv('KO_collagen_normalized.csv',stringsAsFactors = F)

WT_embryos<-read.csv('WT_laminin_normalized.csv',stringsAsFactors = F)
KO_embryos<-read.csv('KO_laminin_normalized.csv',stringsAsFactors = F)

##sliding window of 0.3
#WT analysis
windows<-1:1:100/100
windows<-c(0,windows)
embryos<-unique(WT_embryos$Embryo)
window_data_wt<-data.frame(Distance=as.numeric(),Intensity=as.numeric(),Embryo=as.character(),Condition=as.character())

for (g in 1:length(embryos)){
  emb_data<-WT_embryos[WT_embryos$Embryo==embryos[g],]
  for (i in 1:(length(windows)-1)){
    temp_values<-emb_data[emb_data$Distance>windows[i] & emb_data$Distance<windows[i+1],]
    temp_intensity<-mean(temp_values$Intensity)
    
    window_data_wt<-rbind(window_data_wt,data.frame(windows[i],temp_intensity,embryos[g],"WT"))
  }
}

colnames(window_data_wt)<-c("Distance","Intensity","Embryo","Condition")

window_wt_lines<-matrix(,nrow=length(windows)-1,ncol = 3)

for (i in 1:(length(windows)-1)){
  temp_data<-window_data_wt[window_data_wt$Distance==windows[i],]
  temp_mean<-mean(temp_data$Intensity)
  temp_SD<-sd(temp_data$Intensity)
  window_wt_lines[i,]<-c(windows[i],temp_mean,temp_SD)
}
colnames(window_wt_lines)<-c("Distance","Mean","SD")
window_wt_lines<-as.data.frame(window_wt_lines)

wt_mean<-as.data.frame(spline(window_wt_lines$Distance,window_wt_lines$Mean))
wt_SDup<-as.data.frame(spline(window_wt_lines$Distance,(window_wt_lines$Mean+window_wt_lines$SD)))
wt_SDdown<- as.data.frame(spline(window_wt_lines$Distance,(window_wt_lines$Mean-window_wt_lines$SD)))


ggplot(window_data_wt, aes(x=Distance ,y=Intensity )) + geom_point() + geom_smooth(se=F)+
  geom_smooth(data = wt_SDup, aes(x=x, y=y), color = "darkred",se=F)+
  geom_smooth(data = wt_SDdown, aes(x=x, y=y), color = "darkred",se=F)+
  theme(panel.background = element_blank())
ggsave('WT_collagen_scatter.pdf', height = 10, width=12)

##KO analysis
windows<-1:1:100/100
windows<-c(0,windows)
embryos<-unique(KO_embryos$Embryo)
window_data_KO<-data.frame(Distance=as.numeric(),Intensity=as.numeric(),Embryo=as.character(),Condition=as.character())

for (g in 1:length(embryos)){
  emb_data<-KO_embryos[KO_embryos$Embryo==embryos[g],]
  for (i in 1:(length(windows)-1)){
    temp_values<-emb_data[emb_data$Distance>windows[i] & emb_data$Distance<windows[i+1],]
    temp_intensity<-mean(temp_values$Intensity)
    
    window_data_KO<-rbind(window_data_KO,data.frame(windows[i],temp_intensity,embryos[g],"KO"))
  }
}

colnames(window_data_KO)<-c("Distance","Intensity","Embryo","Condition")

window_KO_lines<-matrix(,nrow=length(windows)-1,ncol = 3)

for (i in 1:(length(windows)-1)){
  temp_data<-window_data_KO[window_data_KO$Distance==windows[i],]
  temp_mean<-mean(temp_data$Intensity)
  temp_SD<-sd(temp_data$Intensity)
  window_KO_lines[i,]<-c(windows[i],temp_mean,temp_SD)
}
colnames(window_KO_lines)<-c("Distance","Mean","SD")
window_KO_lines<-as.data.frame(window_KO_lines)

KO_mean<-as.data.frame(spline(window_KO_lines$Distance,window_KO_lines$Mean))
KO_SDup<-as.data.frame(spline(window_KO_lines$Distance,(window_KO_lines$Mean+window_KO_lines$SD)))
KO_SDdown<- as.data.frame(spline(window_KO_lines$Distance,(window_KO_lines$Mean-window_KO_lines$SD)))


ggplot(window_data_KO, aes(x=Distance ,y=Intensity )) + geom_point() + geom_smooth(se=F)+
  geom_smooth(data = KO_SDup, aes(x=x, y=y), color = "darkred",se=F)+
  geom_smooth(data = KO_SDdown, aes(x=x, y=y), color = "darkred",se=F)+
  theme(panel.background = element_blank())
ggsave('KO_collagen_scatter.pdf', height = 10, width=12)


ggplot(window_data_KO, aes(x=Distance ,y=Intensity ,group = Embryo))  + 
  geom_smooth(data=window_data_KO,aes(x=Distance, y=Intensity,group = NA),se=F)+
  geom_smooth(data = KO_SDup, aes(x=x, y=y,group = NA), color = "darkred",se=F)+
  geom_smooth(data = KO_SDdown, aes(x=x, y=y,group = NA), color = "darkred",se=F)+
  geom_line(data = window_data_KO, aes(x=Distance ,y=Intensity ,group = Embryo), color = "grey")+
  theme(panel.background = element_blank())
ggsave('KO_collagen_line_smoothed.pdf', height = 10, width=12)


ggplot(window_data_KO, aes(x=Distance ,y=Intensity ,group = Embryo))  + 
  geom_line(data = KO_SDup, aes(x=x, y=y,group = NA), color = "darkred",se=F)+
  geom_line(data = KO_SDdown, aes(x=x, y=y,group = NA), color = "darkred",se=F)+
  geom_line(data = window_data_KO, aes(x=Distance ,y=Intensity ,group = Embryo), color = "grey")+
  theme(panel.background = element_blank())
ggsave('KO_collagen_line_raw.pdf', height = 10, width=12)


ggplot(window_data_wt, aes(x=Distance ,y=Intensity ,group = Embryo))  + 
  geom_line(data = wt_SDup, aes(x=x, y=y,group = NA), color = "darkred",se=F)+
  geom_line(data = wt_SDdown, aes(x=x, y=y,group = NA), color = "darkred",se=F)+
  geom_line(data = window_data_wt, aes(x=Distance ,y=Intensity ,group = Embryo), color = "grey")+
  theme(panel.background = element_blank())
ggsave('wt_collagen_line_raw.pdf', height = 10, width=12)


