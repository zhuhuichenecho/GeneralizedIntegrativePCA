library(RColorBrewer)
library(ggplot2)
library(gplots)
library(graphics)

rm(list = ls())


load("SwitzerlandItaly.RData")
load("SwitzerlandItalysize.RData")

data = list()
data[[1]] = Death_data[,1:91]
data[[2]] = Death_data[,112:202]


N.size = list()
N.size[[1]] = size_data[,1:91]
N.size[[2]] = size_data[,112:202]

n.size = do.call(cbind,N.size)

DeathRate = (Death_data/size_data)[,c(1:91,112:202)]

col = rep("Regular",nrow(DeathRate))
col[47] = "Spanish Flu"
col[43:46]  = "WWI"
col[68:74] = "WWII"

col1 = col
col1[43:46] = "WWI"
col1[68:74] = "WWII"

lwd = rep(1,nrow(DeathRate))
lwd[47] = 2

lwd1 = lwd
lwd1[43:46] = 2
lwd1[68:74] = 2

DeathRateItaly = as.vector(t(DeathRate[,1:91]))
year = rep(rownames(DeathRate),each = ncol(DeathRate)/2)
age = rep(as.numeric(colnames(DeathRate)[1:91]),times = nrow(DeathRate))
COL1 = rep(col1,each = ncol(DeathRate)/2)
DeathRateItaly = data.frame(DeathRateItaly,year,age)
size =  rep(1,length(COL1))
size[COL1=="Regular"] = 1
size[COL1!="Regular"] = 1.05
linetype = rep("a",length(COL1))
linetype[COL1=="WWI"] = "b"
linetype[COL1=="WWII"] = "d"
linetype[COL1=="Spanish Flu"] = "c"



######## Figure 2a ########
geom_line(data = DeathRateItaly[which(COL1=="Regular"),],aes(x = age,y = DeathRateItaly,
                                                             group = year,colour = COL1[which(COL1=="Regular")]),size = 0.5)+
  geom_line(data = DeathRateItaly[-which(COL1=="Regular"),],aes(x = age,y = DeathRateItaly,
                                                                group = year,colour = COL1[which(COL1!="Regular")], linetype =  
                                                                  linetype[which(COL1!="Regular")]),size = 1.05)

p1 = ggplot()+ geom_line(data = DeathRateItaly,aes(x = age,y = DeathRateItaly,
                                                   group = year,colour = COL1, linetype =  
                                                     linetype),size = 1)+
  labs(x = "Age",y = "Italy Mortality Rate")+
  scale_y_continuous(limits = range(DeathRate,na.rm =T))+
  scale_colour_manual(name = "", labels = c("Regular","Spanish flu","WWI","WWII"),
                      values=c("black","red","blue","green"))+
  scale_linetype_manual(name = "",labels = c("Regular","Spanish flu","WWI","WWII"),
                        values = c("solid","dotdash","dashed","dotted"))+theme(legend.text=element_text(size=30,face="bold"),
                                                                               axis.text=element_text(size=30),
                                                                               axis.title=element_text(size=30,face="bold"), 
                                                                               legend.key = element_rect(size = 5),
                                                                               legend.key.size = unit(2, 'lines'),
                                                                               legend.position = c(0.15,0.9))+
  geom_line(data = DeathRateItaly[-which(COL1=="Regular"),],aes(x = age,y = DeathRateItaly,
                                                                group = year,colour = COL1[which(COL1!="Regular")], linetype =  
                                                                  linetype[which(COL1!="Regular")]),size = 1.2)


pdf("MortRateVSAgeItaly.pdf",width = 11, height = 8.5)
p1
dev.off()


######## Figure 2b ########
DeathRateSwitz = as.vector(t(DeathRate[,-(1:91)]))
year = rep(rownames(DeathRate),each = ncol(DeathRate)/2)
age = rep(as.numeric(colnames(DeathRate)[1:91]),times = nrow(DeathRate))
COL2 = rep(col,each = ncol(DeathRate)/2)
DeathRateSwitz = data.frame(DeathRateSwitz,year,age)
size =  rep(1,length(COL2))
size[COL2=="Regular"] = 1
size[COL2!="Regular"] = 1.05

p2 = ggplot()+geom_line(data = DeathRateSwitz ,aes(x = age,y = DeathRateSwitz,
                                                   group = year,colour = COL2 ,linetype =  
                                                     linetype ),size = 1)+
  labs(x = "Age",y = "Switzerland Mortality Rate")+
  scale_y_continuous(limits = range(DeathRate,na.rm =T))+
  scale_colour_manual(name = "", labels = c("Regular","Spanish flu","WWI","WWII"),
                      values=c("black","red","blue","green"))+
  scale_linetype_manual(name = "",labels = c("Regular","Spanish flu","WWI","WWII"),
                        values = c("solid","dotdash","dashed","dotted"))+theme(legend.text=element_text(size=30,face="bold"),
                                                                               axis.text=element_text(size=30),
                                                                               axis.title=element_text(size=30,face="bold"), 
                                                                               legend.key = element_rect(size = 5),
                                                                               legend.key.size = unit(2, 'lines'),
                                                                               legend.position = c(0.15,0.9))

pdf("MortRateVSAgeSwitz.pdf",width = 11, height = 8.5)
p2
dev.off()




family = c("binomial","binomial")

tol = 0.1
P1 = ncol(data[[1]])
P2 = ncol(data[[2]])
D = c(P1,P2)
P = sum(D)
N = nrow(data[[1]])
X = do.call(cbind,data)

rank = c(0,0,0)
criteria_cur =  EPCAJIVEinterPoissonMissBIC(rank[1],rank[-1],0.1); criteria_pre = Inf

library(doParallel)
registerDoParallel(cores = 3)

Criteria = NULL
while(criteria_pre-criteria_cur>=0){
  
  print(rank)
  criteria_pre = criteria_cur
  
  calCri = function(i,rank){
    
    source("functions.R")
    X = do.call(cbind,data)
    rank.new = rank
    
    if(i>3){
      
      rank.new[i-3] = ifelse(rank.new[i-3]>0,rank.new[i-3]-1,0)
      
    }else{
      
      rank.new[i] = rank.new[i]+1
      
    }
    
    stepwise(rank.new[1],rank.new[-1],0.1)
    
  }
  
  criteria = foreach(i = 1:(length(rank)),.combine = c,
                     .packages="Matrix")%dopar%calCri(i,rank)
  Criteria = c(Criteria,criteria)
  
  if(which.min(criteria)>3){
    
    rank[which.min(criteria)-3] = ifelse(rank[which.min(criteria)-3]>0,rank[which.min(criteria)-3]-1,0)
    
  }else{
    
    rank[which.min(criteria)] = rank[which.min(criteria)]+1
    
  }
  
  criteria_cur = min(criteria)
  print(criteria)

  
}

rank[which.min(criteria)] = rank[which.min(criteria)]-1

rank = c(1,1,0)
result = EPCAJIVEMissbio(data,rank[1],rank[-1],D = D,family = family, tol = tol,max.iter = 500,n.size = n.size)
Theta.est = result$Theta.joint+result$Theta.inds
Theta.est[which(is.na(result$Theta.inds),arr.ind = T)] = result$Theta.joint[which(is.na(result$Theta.inds),arr.ind = T)]
colnames(Theta.est) = colnames(DeathRate)
rownames(Theta.est) = rownames(DeathRate)

colnames(result$U.joint) = c("1","U.join1")
rownames(result$U.joint) = rownames(DeathRate)
rownames(result$V.joint) = c("Int","V.joint1")
colnames(result$V.joint) = colnames(DeathRate)
rownames(result$U.ind[[1]]) = rownames(DeathRate)
colnames(result$V.ind[[1]]) = colnames(DeathRate)[1:91]
colnames(result$Theta.inds) = colnames(DeathRate)
rownames(result$Theta.inds) = rownames(DeathRate)
colnames(result$Theta.joint) = colnames(DeathRate)
rownames(result$Theta.joint) = rownames(DeathRate)

MU = matrix(result$mu,nrow = nrow(DeathRate),ncol = ncol(DeathRate),byrow = T)
rownames(MU) = rownames(DeathRate)
colnames(MU) = colnames(DeathRate)

EstProp = exp(Theta.est)/(1+exp(Theta.est))

########## Heatmap ##########

Deathrate = data.frame(year = rownames(DeathRate),DeathRate[,1:91])
names(Deathrate)[-1] = colnames(DeathRate[,1:91])
heatrateItaly = melt(Deathrate,id.vars = "year")
names(heatrateItaly) = c("Year","Age","Death_Rate")

dat = data.frame(x=runif(1), y=runif(1))


######## Figure 2c ########
pdf("TrueProp1.pdf",width = 10)
ggplot(heatrateItaly, aes(factor(Age),factor(Year) )) +
  geom_tile(aes(fill = Death_Rate)) +
  # scale_fill_gradient2()+
  scale_fill_gradient(low = "lightblue",high = "darkblue") +
  ylab("Year ") +
  xlab("Age") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Death Rate")+
  geom_point(aes(x=24, y=45), data=dat, size=30, shape=1, color="red")+
  annotate("text", x = 24, y = 50, label = "WWI",color = "red")+
  geom_point(aes(x=24, y=73), data=dat, size=30, shape=1, color="red")+
  annotate("text", x = 24, y = 80, label = "WWII",color = "red")+
  geom_segment(aes(x = 50, xend = 50, y = 55,yend = 48),size = 0.5,color = "red",arrow = arrow(length = unit(0.05,"cm")))+
  annotate("text",x = 60,y = 52,label = "Flu Pademic",color = "red")+
  scale_x_discrete(breaks = seq(0,90, by = 10)) +
  scale_y_discrete(breaks = seq(1872,2014, by = 10))
dev.off()

Deathrate = data.frame(year = rownames(DeathRate),DeathRate[,92:182])
names(Deathrate)[-1] = colnames(DeathRate[,92:182])
heatrateSwiz = melt(Deathrate,id.vars = "year")
names(heatrateSwiz) = c("Year","Age","Death_Rate")


######## Figure 2d ########
pdf("TrueProp2.pdf",width = 10)
ggplot(heatrateSwiz, aes(factor(Age),factor(Year) )) +
  geom_tile(aes(fill = Death_Rate)) +
  scale_fill_gradient(low = "lightblue",high = "darkblue")+
  ylab("Year ") +
  xlab("Age") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Death Rate")+
  geom_segment(aes(x = 25, xend = 25, y = 55,yend = 48),size = 0.5,color = "red",arrow = arrow(length = unit(0.05,"cm")))+
  annotate("text",x = 35,y = 52,label = "Flu Pademic",color = "red")+
  scale_x_discrete(breaks = seq(0,90, by = 10)) +
  scale_y_discrete(breaks = seq(1872,2014, by = 10))
dev.off()


######## Figure 3a ########
Italymean = data.frame(age = as.numeric(colnames(Death_data)[1:91]),mean = result$mu[1:91])
p1 = ggplot(data = Italymean)+geom_line(aes(x = age,y  = mean,group = 1),size = 1.5)+
  labs(x = "Age",y = "Italy estimated mean")+theme(legend.text=element_text(size=16),
                                                   axis.text=element_text(size=30),
                                                   axis.title=element_text(size=30,face="bold"))

pdf("Italymean.pdf",width = 11, height = 8.5)
p1
dev.off()


######## Figure 3b ########
Switzmean = data.frame(age = as.numeric(colnames(Death_data)[1:91]),mean = result$mu[92:182])
p2 = ggplot(data = Switzmean)+geom_line(aes(x = age,y  = mean,group = 1),size = 1.5)+
  labs(x = "Age",y = "Switzerland estimated mean")+theme(legend.text=element_text(size=30),
                                                         axis.text=element_text(size=30),
                                                         axis.title=element_text(size=30,face="bold"))

pdf("Switzmean.pdf",width = 11, height = 8.5)
p2
dev.off()



######## Figure 3c ########
U.joint = result$U.joint[,2]*svd(scale(result$Theta.joint,scale = F),nu = 1,nv = 1)$d[1]

jointscore = data.frame(U.joint = U.joint,year = as.numeric(rownames(Death_data)))
rec = data.frame(x = c(1917,1917,1919,1919),y = c(min(U.joint[46:48])-1,max(U.joint[46:48])+1,
                                                  max(U.joint[46:48])+1,min(U.joint[46:48])-1))
p3 = ggplot(data = jointscore)+geom_line(aes(x = year,y = U.joint,group = 1),size = 1.5)+
  labs(x = "Year",y = "Joint score")+geom_line(data = rec[1:2,],aes(x = x,y = y,group = 1),size = 1.2,color = "red")+
  scale_x_continuous(breaks = seq(1885,2014,25))+
  geom_line(data = rec[3:4,],aes(x = x,y = y,group = 1),size = 1.2,color = "red")+
  annotate("text",x = 1918,y = -3,label = "Spanish Flu",color = "black",size = 8,fontface = 2)+theme(legend.text=element_text(size=16),
                                                                                                     axis.text=element_text(size=30),
                                                                                                     axis.title=element_text(size=30,face="bold"))
geom_polygon(data = rec,aes(x = x,y = y,group = 1),
             color="grey20",
             alpha=0.3,
             inherit.aes = FALSE)

pdf("Ujoint.pdf",width = 11, height = 8.5) 
p3
dev.off()



######## Figure 3d ########
jointloadItaly = data.frame(age = as.numeric(colnames(Death_data)[1:91]),load = result$V.joint[2,1:91])
p4 = ggplot(data = jointloadItaly)+geom_line(aes(x = age,y = load,group = 1),size = 1.5)+
  labs(x = "Age",y = "Joint loading for Italy")+theme(legend.text=element_text(size=16),
                                                      axis.text=element_text(size=30),
                                                      axis.title=element_text(size=30,face="bold"))

pdf("ItalyVjoint.pdf",width = 11, height = 8.5) 
p4
dev.off()


######## Figure 3e ########
jointloadSwitz = data.frame(age = as.numeric(colnames(Death_data)[1:91]),load = result$V.joint[2,92:182])
p5 = ggplot(data = jointloadSwitz)+geom_line(aes(x = age,y = load,group = 1),size = 1.5)+
  labs(x = "Age",y = "Joint loading for Switzerland")+theme(legend.text=element_text(size=16),
                                                            axis.text=element_text(size=30),
                                                            axis.title=element_text(size=30,face="bold"))
pdf("SwitzVjoint.pdf",width = 11, height = 8.5) 
p5
dev.off()


######## Figure 3f ########
indscoreItaly = data.frame(year =as.numeric(rownames(Death_data)),score = result$U.ind[[1]])
rec1 = data.frame(x = c(1914,1914,1918,1918),y = c(min(result$U.ind[[1]][43:47])-1,max(result$U.ind[[1]][43:47])+1,
                                                   max(result$U.ind[[1]][43:47])+1,min(result$U.ind[[1]][43:47])-1))
rec2 = data.frame(x = c(1939,1939,1945,1945),y = c(min(result$U.ind[[1]][68:74])-1,max(result$U.ind[[1]][68:74])+1,
                                                   max(result$U.ind[[1]][68:74])+1,min(result$U.ind[[1]][68:74])-1))
p6 = ggplot(data = indscoreItaly)+geom_line(aes(x = year,y = score,group = 1),size = 1.5)+
  geom_line(data = rec1[1:2,],aes(x = x,y = y),size = 1.2,color = "red")+
  geom_line(data = rec1[3:4,],aes(x = x,y = y),size = 1.2,color = "red")+
  geom_line(data = rec2[1:2,],aes(x = x,y = y),size = 1.2,color = "red")+
  geom_line(data = rec2[3:4,],aes(x = x,y = y),size = 1.2,color = "red")+
  scale_x_continuous(breaks = seq(1885,2014,25))+
  labs(x = "Year",y = "Individual score for Italy")+
  annotate("text",x = 1915,y = 1,label = "WWI",color = "black",size = 8,fontface = 2)+
  annotate("text",x = 1942,y = 2,label = "WWII",color = "black",size = 8,fontface = 2)+theme(legend.text=element_text(size=16),
                                                                                             axis.text=element_text(size=30),
                                                                                             axis.title=element_text(size=30,face="bold"))

+geom_polygon(data = rec1,aes(x = x,y = y,group = 1),
              color="grey20",
              alpha=0.3,
              inherit.aes = FALSE)+
  +geom_polygon(data = rec2,aes(x = x,y = y,group = 1),
                color="grey20",
                alpha=0.3,
                inherit.aes = FALSE)

pdf("ItalyUind.pdf",width = 11, height = 8.5)
p6
dev.off()



######## Figure 3g ########
indloadSwitz = data.frame(age = as.numeric(colnames(Death_data)[1:91]),load = as.vector(result$V.ind[[1]]))
p7 = ggplot(indloadSwitz)+geom_line(aes(x = age,y = load,group = 1),size = 1.5)+
  labs(x = "Age",y = "Indvidual loading for Italy")+theme(legend.text=element_text(size=16),
                                                          axis.text=element_text(size=30),
                                                          axis.title=element_text(size=30,face="bold"))

pdf("ItalyVind.pdf",width = 11, height = 8.5)
p7
dev.off()