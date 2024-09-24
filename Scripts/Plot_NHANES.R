

setwd("~/Desktop/Rugao/Longevity//SWEET/NHANES/Individual/")

library(readr)
lf=list.files()
dl=matrix(0,nrow = 253,ncol=length(lf))

k=1
for(f in lf[grep("txt",lf)]){
  d=read.table(f,header = T)
  dl[,k]=as.numeric(d[[3]])
  k=k+1
}

rownames(dl)=paste(d$gene1,d$gene2,sep="--")
colnames(dl)=gsub(".txt","",lf)

info=readRDS("../NHANES.info.rds")
table(info$status)
table(info$year)
info$id=sapply(info$sampleID,function(x) stringr::str_split(x,"_")[[1]][2])
info$newid=paste(as.numeric(as.factor(info$year)),info$id,sep="_")
nhanes_merged1=readRDS("../NHANES.DATA.rds")
nhanes_merged1$id=paste(nhanes_merged1$SDDSRVYR,nhanes_merged1$SEQN,sep="_")
info$status1=nhanes_merged1$MORTSTAT[match(info$newid,nhanes_merged1$id)]
info$time1=nhanes_merged1$PERMTH_INT[match(info$newid,nhanes_merged1$id)]
info$cause=nhanes_merged1$UCOD_LEADING[match(info$newid,nhanes_merged1$id)]
info$RDW=nhanes_merged1$LBXRDW[match(info$newid,nhanes_merged1$id)]


info=info[match(colnames(dl),info$sampleID),]
dim(info)
dim(info[info$age<30,])

dl1=dl
for(i in 1:nrow(dl)){
  intl=quantile(dl[i,][which(info$age<30)], probs = c(0,0.05,0.5,0.95,1)) #[which(info$trueage<80)]
  tt=rep(0,ncol(dl))
  tt[dl[i,]<intl[2]]=1
  tt[dl[i,]>intl[4]]=1
  dl1[i,]=tt
}

Biomarker1=sapply(rownames(dl1),function(x) strsplit(x,"--")[[1]][1])
Biomarker2=sapply(rownames(dl1),function(x) strsplit(x,"--")[[1]][2])

library(igraph)
###################################################################
cal_deg <- function(kk){
  dg=data.frame(Biomarker1,Biomarker2,
                w=kk)
  dg=dg[dg$w>0,]
  g=graph.data.frame(dg, directed = F,
                     vertices = unique(c(Biomarker1,Biomarker2)))
  degree(g)
}

net=dl1
net.deg=apply(net, 2, cal_deg)


phe=c("albumin","alp","bun","creat","glucose","ttbl",
      "uap","basopa","eosnpa","lymph","monopa","neut","wbc","crp",
      "rbc","mcv","rdw","ggt","dbp","sbp","hdl","totchol","hba1c")
length(phe)

apply(net.deg, 1, function(x) cor.test(x,info$age,method = "spearman")$estimate)
cor.test(info$age,net[2,])

cs=apply(info[,colnames(info) %in% phe], 2, function(x) cor.test(x,info$age,method = "spearman")$estimate)

names(cs)[abs(cs)>0.2]

rownames(net.deg)

##################################
cal_topo <- function(kk){
  dg=data.frame(Biomarker1,Biomarker2,
                w=kk)
  dg=dg[dg$w>0,]
  g=graph.data.frame(dg, directed = F,
                     vertices = unique(c(Biomarker1,Biomarker2)))
  g1=graph.data.frame(dg)
  g=igraph::simplify(g)
  #c1 = cluster_leading_eigen(g)
  c1 = cluster_fast_greedy(g)
  ll=c( average.degree = mean(igraph::degree(g)),
        average.path.length = average.path.length(g), 
        #edge_num = ecount(g),
        vertex = vcount(g1),
        connectance=edge_density(g1,loops=FALSE),
        #average_density=graph.density(g1),
        #diameter = diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL),
        #connectivity = edge_connectivity(g),
        #
        #edge.connectivity = edge_connectivity(g),
        #
        
        centralization.betweenness = centralization.betweenness(g)$centralization,
        #centralization.eigen=centr_eigen(g, directed = FALSE)$centralization,
        #centralization.pg=mean(page_rank(g)$vector),
        #centr_clo(g, mode = "all")$centralization
        centralization.degree = centralization.degree(g)$centralization,
        clustering.coefficient = transitivity(g,type = "global"),
        Q=modularity(c1)
  )
  return(ll)
}

net.stat=apply(dl1, 2, cal_topo)

edges=apply(dl1, 2, sum)

cor.test(info$age,edges,method = "spearman")
#plot(info$age,edges)
plot(density(info$age,na.rm = T))

setwd("~/Desktop/Rugao/Longevity//SWEET/NHANES/")
library(haven)

##################################
library(readr)

pd.vd=data.frame(age=info$age,
                 wave=info$wave,
                 id=info$sampleID,
                 #id1=sapply(info$sampleID,function(x) strsplit(x,"_")[[1]][2]),
                 sex=info$gender,
                 edu=info$edu,
                 race=info$race,
                 time=info$time1,
                 status=info$status1,
                 adl=info$adl,
                 walk=info$lnwalk,
                 grip=info$grip_scaled,
                 health=info$health,
                 cause=info$ucod_leading,
                 hyperten=info$hyperten,
                 diabetes=info$diabetes)

pd.vd=cbind(pd.vd,t(net.stat))
pd.vd=pd.vd[which(pd.vd$status>=0),]
dim(pd.vd)
#write.csv(pd.vd,"NHANES.NB.csv",row.names = F)

dim(pd.vd[!is.na(pd.vd$walk),])
mean(pd.vd$age[!is.na(pd.vd$walk)])
sd(pd.vd$age[!is.na(pd.vd$walk)])

dim(pd.vd[!is.na(pd.vd$grip),])
mean(pd.vd$age[!is.na(pd.vd$grip)])
sd(pd.vd$age[!is.na(pd.vd$grip)])

pd.vd$ad=cut(pd.vd$average.degree,
             breaks=quantile(pd.vd$average.degree,probs = c(0,0.5,1)),
             include.lowest = T,right = T)

pd.vd$adly=cut(pd.vd$adl,breaks = c(0,1,20),include.lowest = T,right = F)
l1=glm(adly~average.degree,data = pd.vd,family = "binomial")
s1=summary(l1)
s1

l1=glm(adly~age+sex+average.degree,data = pd.vd,family = "binomial")
s1=summary(l1)
s1

l1=lm(health~age+sex+average.degree,data = pd.vd)
s1=summary(l1)
s1

l1=lm(walk~age+sex+average.degree,data = pd.vd)
s1=summary(l1)
s1

l1=lm(grip~age+sex+average.degree,data = pd.vd)
s1=summary(l1)
s1

#################
## Age stratification
pd.vd$ad=cut(pd.vd$average.degree,
             breaks=quantile(pd.vd$average.degree,probs = c(0,0.5,1)),
             include.lowest = T,right = T)

pd.vd$cc=cut(pd.vd$clustering.coefficient,
             breaks=quantile(pd.vd$clustering.coefficient,probs = c(0,0.5,1),na.rm = T),
             include.lowest = T,right = T)

##############
pd.vd$time1=pd.vd$time/12
km.to <- survfit(Surv(time1,status)~ad,data = pd.vd )
ggsurvplot(km.to,data=pd.vd,pval = TRUE)

km.to1 <- survfit(Surv(time1,status)~cc,data = pd.vd )
ggsurvplot(km.to1,data=pd.vd,pval = TRUE)

p1=ggsurvplot(km.to, # 
              data = pd.vd,  
              conf.int = TRUE, # 显示置信区间
              pval = TRUE, # 添加P值
              #risk.table = TRUE, # 绘制累计风险曲线
              #surv.median.line = "hv", # 添加中位生存时间线
              #add.all = TRUE, # 添加总患者生存曲线
              censor = FALSE,
              ylim = c(0.6, 1),
              #xlim = c(15, 96),
              font.x=12,font.y=12,
              xlab="Years",
              legend.title="NHANES",
              legend.labs=c("Low AD","High AD"),
              pval.coord = c(1,0.7),
              palette = "hue")  #+ # 自定义调色板 +

p2=ggsurvplot(km.to1, # 
              data = pd.vd,  
              conf.int = TRUE, # 显示置信区间
              pval = TRUE, # 添加P值
              #risk.table = TRUE, # 绘制累计风险曲线
              #surv.median.line = "hv", # 添加中位生存时间线
              #add.all = TRUE, # 添加总患者生存曲线
              censor = FALSE,
              ylim = c(0.6, 1),
              #xlim = c(15, 96),
              font.x=12,font.y=12,
              xlab="Years",
              legend.title="NHANES",
              legend.labs=c("Low CC","High CC"),
              pval.coord = c(1,0.7),
              palette = "hue")  #+ # 自定义调色板 +

pdf("NHANES.Km.curve.Years.AD.pdf",height = 3,width = 4,onefile = F)
p1$plot
dev.off()

pdf("NHANES.Km.curve.Years.CC.pdf",height = 3,width = 4,onefile = F)
p2$plot
dev.off()

##############


#form = formula(paste("survival::Surv(time,status)~", i, "+", covars, sep = ""))
#cox[[i]] = survival::coxph(form, data = dat)

#cox1 <- coxph(Surv(time,status)~age+sex+edu+race+average.degree,data = pd.vd)
#cox1

table(pd.vd$cause)

covars="age+sex+edu+race"
cl=list(0,1,2,3,5,6,7,8,9)
cl[[1]]=seq(1,10)
cl

ml=c("All-cause","Heart disease","Malignant neoplasms","Respiratory diseases",
     "Cerebrovascular diseases","Alzheimer’s disease","Diabetes mellitus",
     "Influenza and pneumonia","Nephritis/nephrosis")
k=1
pl=list()
pl1=list()
for(j in cl){
  cl1=which(pd.vd$status==0)
  cl2=which(pd.vd$cause %in% c(j))
  pd.vd1=pd.vd[c(cl1,cl2),]
  #print(dim(pd.vd1))
  print(table(pd.vd1$status))
  
  km.to1 <- survfit(Surv(time1,status)~ad,data = pd.vd1 )
  km.to2 <- survfit(Surv(time1,status)~cc,data = pd.vd1 )
  
  p1=ggsurvplot(km.to1, # 
                data = pd.vd1,  
                conf.int = TRUE, # 显示置信区间
                pval = TRUE, # 添加P值
                #risk.table = TRUE, # 绘制累计风险曲线
                #surv.median.line = "hv", # 添加中位生存时间线
                #add.all = TRUE, # 添加总患者生存曲线
                censor = FALSE,
                ylim = c(0.6, 1),
                #xlim = c(15, 96),
                font.x=12,font.y=12,
                xlab="Years",
                title=ml[k],
                legend.title="NHANES",
                legend.labs=c("Low AD","High AD"),
                pval.coord = c(1,0.7),
                palette = "hue")  #+ # 自定义调色板 +
  
  p2=ggsurvplot(km.to2, # 
                data = pd.vd1,  
                conf.int = TRUE, # 显示置信区间
                pval = TRUE, # 添加P值
                #risk.table = TRUE, # 绘制累计风险曲线
                #surv.median.line = "hv", # 添加中位生存时间线
                #add.all = TRUE, # 添加总患者生存曲线
                censor = FALSE,
                ylim = c(0.6, 1),
                #xlim = c(15, 96),
                font.x=12,font.y=12,
                xlab="Years",
                title=ml[k],
                legend.title="NHANES",
                legend.labs=c("Low AD","High AD"),
                pval.coord = c(1,0.7),
                palette = "hue")  #+ # 自定义调色板 +
  
  pl[[k]]=p1
  pl1[[k]]=p2
  k=k+1
}

library(ggpubr)

#1: diseases of heart
#2: malignant neoplasms
#3: chronic lower respiratory diseases
#4: accidents
#5: cerebrovascular diseases
#6: Alzheimer’s disease
#7: diabetes mellitus
#8: influenza and pneumonia
#9: nephritis, nephrotic syndrome and nephrosis
#10: all other causes

od=data.frame()
k=1
for(j in cl){
  for(i in rownames(net.stat)){
    cl1=which(pd.vd$status==0)
    cl2=which(pd.vd$cause %in% c(j))
    pd.vd1=pd.vd[c(cl1,cl2),]
    table(pd.vd1$status)
    
    form = formula(paste("survival::Surv(time,status)~", i, "+", covars, sep = ""))
    cox = survival::coxph(form, data = pd.vd1)
    
    s=summary(cox)
    p=s$coefficients[1,5]
    or=s$conf.int[1,1]
    #print(paste(i,": ",round(or,2)," (",paste(ci,collapse = "-"),"), ",p,sep=""))
    
    beta=round(s$coefficients[1,1],2)
    se=s$coefficients[1,3]
    #print(paste(i,": ",beta,se,p))
    #print(paste(i,": ",or,ci,p))
    res=c(i,beta,se,p,ml[k])
    od=rbind(od,res)
  }
  k=k+1
}

colnames(od)=c("name","beta","se","p","cause")
od$beta=as.numeric(od$beta)
od$se=as.numeric(od$se)
od$name[od$name=="Q"]="modularity"
od$name=gsub("[.]"," ",od$name)

write.csv(od,"NHANES.Mortality.Cause.csv",quote = F)

od$lb=cut(as.numeric(od$p),breaks = c(0,0.001,0.01,0.05,1),right = T)
table(od$lb)
levels(od$lb)=c("P < 0.001","P < 0.01","P < 0.05","NS")

od$name=factor(od$name,levels=unique(od$name)[c(1,3,4,2,8,7,5,6)])
od$cause=factor(od$cause,levels=rev(unique(od$cause)))
library(ggforestplot)

pdf("Mortality_Net_topo.pdf",height = 6,width = 10)
ggplot(od,aes(beta,cause,color=lb)) +
  geom_point() + 
  geom_errorbar(aes(xmin=beta-se,xmax=beta+se), width = 0.5)+
  geom_vline(xintercept = 0, color = "black",linetype="dashed",alpha=0.6)+
  facet_wrap(name~.,scales = "free_x",nrow = 2) +
  theme_bw() + labs(color="Significance") +
  scale_x_continuous(n.breaks = 4) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        panel.spacing.x = unit(0.3,"lines"),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=12),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")) +
  scale_color_npg() + ylab("")
dev.off()  

###########
pd.vd$ageg=cut(pd.vd$age,breaks = c(18,39,59,90),right = T)
table(pd.vd$ageg)
#levels(pd.vd$ageg)=c("Young (18-44 years)","Middle-Aged (45-64 years)","Old-aged (65-85 years)") 
levels(pd.vd$ageg)=c("Young","Middle-Aged","Old-aged") 


pd.vd$sexg=as.factor(pd.vd$sex)
levels(pd.vd$sexg)=c("Male","Female")

pd.vd$edug=cut(pd.vd$edu,breaks=c(0,1,2,4),right = T)
table(pd.vd$edug)
levels(pd.vd$edug)=c("Less than HS","HS/GED","Colloge")

pd.vd$raceg=as.factor(pd.vd$race)
levels(pd.vd$raceg)=c("White","Black","Hispanic") 

dl=list(pd.vd$ageg,pd.vd$sexg,pd.vd$edug,pd.vd$raceg)
names(dl)=c("Age","Gender","Education","Race")
od=data.frame()
k=1
for(j in dl){
  for(l in levels(j)){
    pd.vd1=pd.vd[which(j==l),]
    table(pd.vd1$status)
    print(nrow(pd.vd1))
    for(i in rownames(net.stat)){
      form = formula(paste("survival::Surv(time,status)~", i, "+", covars, sep = ""))
      cox = survival::coxph(form, data = pd.vd1)
      
      s=summary(cox)
      p=s$coefficients[1,5]
      or=s$conf.int[1,1]
      #print(paste(i,": ",round(or,2)," (",paste(ci,collapse = "-"),"), ",p,sep=""))
      
      beta=round(s$coefficients[1,1],2)
      se=s$coefficients[1,3]
      #print(paste(i,": ",beta,se,p))
      #print(paste(i,": ",or,ci,p))
      res=c(i,beta,se,p,l,names(dl)[k])
      od=rbind(od,res)
    }
  }
  k=k+1
}

colnames(od)=c("name","beta","se","p","Groups","Category")
od$beta=as.numeric(od$beta)
od$se=as.numeric(od$se)
od$name[od$name=="Q"]="Modularity"
od$name=gsub("[.]"," ",od$name)

write.csv(od,"NHANES.Mortality.Subgroups.csv",quote = F)


od$lb=cut(as.numeric(od$p),breaks = c(0,0.001,0.01,0.05,1),right = T)
table(od$lb)
levels(od$lb)=c("P < 0.001","P < 0.01","P < 0.05","NS")

od$name=factor(od$name,levels=unique(od$name)[c(1,3,4,2,8,7,5,6)])
od$Groups=factor(od$Groups,levels=rev(unique(od$Groups)))
library(ggforestplot)
library(ggh4x)

p1=ggplot(od[od$name %in% levels(od$name)[1:4],],aes(beta,Groups,color=lb)) +
  geom_point() + 
  geom_errorbar(aes(xmin=beta-se,xmax=beta+se), width = 0.5)+
  geom_vline(xintercept = 0, color = "black",linetype="dashed",alpha=0.6)+
  facet_grid(Category~name,scales = "free",space = "free_y") +
  theme_bw() + labs(color="Significance") +
  scale_x_continuous(n.breaks = 4) +
  #scale_y_discrete(position = "right") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        panel.spacing.x = unit(0.3,"lines"),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=12),
        strip.text.y = element_blank(),
        #strip.text.y.left  = element_text(size=12,angle = 0),
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "lines")) +
  scale_color_npg() + ylab("") + xlab("")

p2=ggplot(od[od$name %in% levels(od$name)[5:8],],aes(beta,Groups,color=lb)) +
  geom_point() + 
  geom_errorbar(aes(xmin=beta-se,xmax=beta+se), width = 0.5)+
  geom_vline(xintercept = 0, color = "black",linetype="dashed",alpha=0.6)+
  facet_grid(Category~name,scales = "free",space = "free_y") +
  theme_bw() + labs(color="Significance") +
  scale_x_continuous(n.breaks = 4) +
  #scale_y_discrete(position = "right") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        panel.spacing.x = unit(0.3,"lines"),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=12),
        strip.text.y = element_blank(),
        #strip.text.y.left  = element_text(size=12,angle = 0),
        plot.margin = unit(c(-0.8, 0.5, 0.5, 0.5), "lines")) +
  scale_color_npg() + ylab("") + xlab("Beta")

pdf("Subgroups_Net_topo.pdf",height = 6,width = 10)

ggarrange(p1,p2,nrow = 2,ncol = 1,
          common.legend = T,
          legend = "right")
dev.off()  

################
pd.vd$walky=cut(pd.vd$walk,breaks=quantile(pd.vd$walk,probs = c(0,0.8,1),na.rm = T),
                right = T,include.lowest = T)
table(pd.vd$walky)

th1=quantile(pd.vd$grip[pd.vd$sex==1],probs = c(0,0.2,1),na.rm = T)[2]
th2=quantile(pd.vd$grip[pd.vd$sex==2],probs = c(0,0.2,1),na.rm = T)[2]  

pd.vd$gripy=rep(0,nrow(pd.vd))
pd.vd$gripy[which(pd.vd$sex==1 & pd.vd$grip <th1)]=1
pd.vd$gripy[which(pd.vd$sex==2 & pd.vd$grip <th2)]=1
pd.vd$gripy[is.na(pd.vd$grip)]=NA
table(pd.vd$gripy)

pd.vd$healthy=cut(pd.vd$health,
                  breaks=c(0,3,5),
                  right = T,include.lowest = T)

pd.vd$adly=cut(pd.vd$adl,
               breaks=c(-1,0,20),
               right = T,include.lowest = T)
table(pd.vd$adly)

covars="age+sex+edu+race"
res=data.frame()
for(i in rownames(net.stat)){
  tt=c()
  for(yy in c("adly","walky","gripy")){
    form = formula(paste(yy,"~", i, "+", covars, sep = ""))
    l=glm(form,data=pd.vd,family = "binomial")
    s=summary(l)
    p=s$coefficients[2,4]
    beta=round(s$coefficients[2,1],2)
    se=s$coefficients[2,2]
    tt=c(tt,c(beta,se,p))
  }
  #print(paste(i,": ",beta,se,p))
  res=rbind(res,tt)
}
rownames(res)=rownames(net.stat)
colnames(res)=c("ADL beta","ADL se","ADL p",
                "Walk beta","Walk se","Walk p",
                "Grip beta","Grip se","Grip p")
res=res[c(1,3,4,2,8,7,5,6),]
res

res2=data.frame(b1=paste(sprintf("%.2f",as.numeric(res[,1]))," (",
                         sprintf("%.2f",as.numeric(res[,2])),")",sep=""),
                p1=sprintf("%.3f",as.numeric(res[,3])),
                b2=paste(sprintf("%.2f",as.numeric(res[,4]))," (",
                         sprintf("%.2f",as.numeric(res[,5])),")",sep=""),
                p2=sprintf("%.3f",as.numeric(res[,6])),
                b3=paste(sprintf("%.2f",as.numeric(res[,7]))," (",
                          sprintf("%.2f",as.numeric(res[,8])),")",sep=""),
                p3=sprintf("%.3f",as.numeric(res[,9])))
rownames(res2)=rownames(res)  
res2

write.csv(res,"NHANES.ADL.stat.csv",quote = F)
write.csv(res,"NHANES.ADL.stat2.csv",quote = F)
