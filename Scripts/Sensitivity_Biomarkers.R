

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

thres=c(0.05)
dl1=dl
for(i in 1:nrow(dl)){
  q=0.05
  intl=quantile(dl[i,][which(info$age<30)], probs = c(0,q,0.5,1-q,1)) #[which(info$trueage<80)]
  tt=rep(0,ncol(dl))
  tt[dl[i,]<intl[2]]=1
  tt[dl[i,]>intl[4]]=1
  dl1[i,]=tt
}

library(igraph)
###################################################################
## Sampling the nodes
phe=c("albumin","alp","bun","creat","glucose","ttbl",
      "uap","basopa","eosnpa","lymph","monopa","neut","wbc","crp",
      "rbc","mcv","rdw","ggt","dbp","sbp","hdl","totchol","hba1c")
length(phe)
cc=apply(info[,phe], 2, function(x) cor(info$age,x,use="pairwise.complete.obs",method = "spearman"))
phe1=phe[order(abs(cc),decreasing = T)]
phe1

f1=lf[grep("txt",lf)]
d=read.table(f,header = T)

dlk=list()
l=seq(5,23,1)
k=1
for(nl in l){
  for(j in 1:10){
    bio=phe1[sample(length(phe1),nl)]
    sl=which(d$gene1 %in% bio & d$gene2 %in% bio)
    dlk[[k]]=dl1[sl,]
    k=k+1
  }
}

length(dlk)

##################################
cal_topo <- function(kk){
  Biomarker1=sapply(names(kk),function(x) strsplit(x,"--")[[1]][1])
  Biomarker2=sapply(names(kk),function(x) strsplit(x,"--")[[1]][2])
  dg=data.frame(Biomarker1,Biomarker2,
                w=kk)

  dg=dg[dg$w>0,]
  library(igraph)
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

library(parallel)
# 用system.time来返回计算所需时间

se=sample(nrow(info),0.2*nrow(info))
length(se)

net.stat=list()
for(i in 1:length(dlk)){
  print(paste("####: ",i))
  #net.stat[[i]]=apply(dlk[[i]], 2, cal_topo)
  system.time({
    cl <- makeCluster(16) # 初始化四核心集群
    net.stat[[i]] <- parApply(cl,dlk[[i]][,se], 2, cal_topo) 
    stopCluster(cl)
  })
}

#saveRDS(net.stat,"Biomarkers.robust.rds")
setwd("~/Desktop/Rugao/Longevity//SWEET/NHANES/")
library(haven)
net.stat=readRDS("Biomarkers.robust.rds")

plot(info$age[se],net.stat[[1]][1,])
cor.test(info$age[se],net.stat[[1]][1,])
cor.test(info$age[se],net.stat[[2]][1,])
cor.test(info$age[se],net.stat[[3]][1,])

##################################
library(readr)

cal_cindex <- function(net.stat){
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
  pd.vd=pd.vd[se,]
  pd.vd=cbind(pd.vd,t(net.stat))
  pd.vd=pd.vd[which(pd.vd$status>=0),]
  dim(pd.vd)

  ss=c()
  for(i in rownames(net.stat)){
    #i="average.degree"
    form = formula(paste("survival::Surv(time,status)~", i, sep = ""))
    cox = survival::coxph(form, data = pd.vd)
    
    s=summary(cox)
    ss=c(ss,s$concordance)
  }
  ss
}

cc=matrix(0,ncol = length(net.stat),nrow = 8*2)
colnames(cc)=rep(seq(5,23),each=10)
row.names(cc)=rep(rownames(net.stat[[1]]),2)
ce=cc
for(i in 1:length(net.stat)){
  ci=cal_cindex(net.stat[[i]])  
  cc[,i]=ci[seq(1,32,2)]
  ce[,i]=ci[seq(2,32,2)]
  #print(cc[,i])
}
cc


library(reshape2)
dc=melt(cc[1:8,])
dc1=melt(ce[1:8,])
dc = dc %>% group_by(Var1,Var2) %>%
  summarise(q1=quantile(value)[2],
            q3=quantile(value)[4],
            value=median(value))

dc$Var1=gsub("[.]"," ",as.character(dc$Var1))
dc$Var1[grep("cent",dc$Var1)]=gsub(" ","\n",dc$Var1[grep("cent",dc$Var1)])
dc$Var1[grep("clu",dc$Var1)]=gsub(" ","\n",dc$Var1[grep("clu",dc$Var1)])
dc$Var1[dc$Var1=="average path length"]="average path\nlength" 

dc$Var1[dc$Var1=="Q"]="modularity"
unique(dc$Var1)
dc$Var1=factor(dc$Var1,unique(dc$Var1)[c(1,3,4,2,8,7,5,6)])
pal <- hcl.colors(50,"Zissou1")

pdf("Biomarkers.sel.cindex.pdf",height = 5,width = 10)
ggplot(dc,aes(Var2,value)) +
  geom_line(aes(group=Var1)) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin = q1, ymax = q3),width=0.3) +
  geom_hline(yintercept = 0.5,linetype="dashed") +
  facet_wrap(Var1~.,nrow = 2,scales = "free_x") +
  xlab("Number of biomarkers") + ylab("Survival analysis C-index") + #labs(fill="C-index") +
  theme_bw() +
  scale_y_continuous(limits = c(0.49,0.6),breaks = seq(0.5,0.65,0.05)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=16,face = "bold"),
        panel.spacing.x = unit(0.8,"line"),
        axis.text.x  = element_text(size=14),
        axis.title = element_text(size=16))

dev.off()