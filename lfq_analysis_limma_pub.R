proteomicstable <- read.delim("proteomicstable.txt",stringsAsFactors = F)

entrezrat=read.delim("gene_result.txt",header=T, stringsAsFactors = F)
ensemblrat=read.delim("mart_rat_export.txt",header=T,stringsAsFactors=F)
unirat=read.delim("uniprot_rat.tsv",header=T,stringsAsFactors = F)

unirat$ensembl=""
unirat$genename_ensembl=""
unirat$entrez_ensembl=""
for (i in 1:nrow(ensemblrat)){
  if (ensemblrat$UniProtKB.Swiss.Prot.ID[i]!=""){
    index=which(unirat$Entry==ensemblrat$UniProtKB.Swiss.Prot.ID[i])
    if (length(index)==1){
      unirat$ensembl[index]=ensemblrat$Gene.stable.ID[i]
      unirat$genename_ensembl[index]=ensemblrat$Gene.name[i]
      unirat$entrez_ensembl[index]=ensemblrat$NCBI.gene..formerly.Entrezgene..ID[i]
    }
  }
  if (ensemblrat$UniProtKB.TrEMBL.ID[i]!=""){
    index=which(unirat$Entry==ensemblrat$UniProtKB.TrEMBL.ID[i])
    if (length(index)==1){
      unirat$ensembl[index]=ensemblrat$Gene.stable.ID[i]
      unirat$genename_ensembl[index]=ensemblrat$Gene.name[i]
      unirat$entrez_ensembl[index]=ensemblrat$NCBI.gene..formerly.Entrezgene..ID[i]
    }
  }
}

library(org.Rn.eg.db)
entrezvec=mapIds(org.Rn.eg.db,keys=unirat$Entry,column="ENTREZID",keytype="UNIPROT",multiVals = "first")
unirat$entrez2=as.character(entrezvec)

entrezvec=mapIds(org.Rn.eg.db,keys=unirat$Gene.Names..primary.,column="ENTREZID",keytype="SYMBOL",multiVals = "first")
unirat$entrez3=as.character(entrezvec)


unirat$entrez=""
for (i in 1:nrow(unirat)){
  if (is.element(unirat$Gene.Names..primary.[i],entrezrat$Symbol)){
    index=which(entrezrat$Symbol==unirat$Gene.Names..primary.[i])
    unirat$entrez[i]=entrezrat$GeneID[index[1]]
  } else {
      has_alias=grep(unirat$Gene.Names..primary.[i],entrezrat$Aliases,fixed=T)
      if (length(has_alias)==1){
        unirat$entrez[i]=entrezrat$GeneID[has_alias]
      }  
  }
}


allentrez=unirat[,c(13,18,19,20,21)]
allentrez[allentrez=="NULL"]=""
allentrez[is.na(allentrez)]=""
colSums(allentrez=="")
aaa=rowSums(allentrez=="")
sum(aaa==5)

for (i in 1:nrow(allentrez)){
  entrezstr=unlist(strsplit(allentrez[i,1],split=";",fixed=T))
  allentrez[i,1]=entrezstr[1]
}
allentrez[is.na(allentrez)]=""

allentrez$join=allentrez[,5]
allentrez$join[allentrez$join==""]=allentrez[allentrez$join=="",4]
allentrez$join[allentrez$join==""]=allentrez[allentrez$join=="",3]
allentrez$join[allentrez$join==""]=allentrez[allentrez$join=="",2]
allentrez$join[allentrez$join==""]=allentrez[allentrez$join=="",1]
colSums(allentrez=="")

unirat$entrezjoin=allentrez$join


proteomicstable$gene_name="-"
genenames=list()
for (i in 1:nrow(proteomicstable)){
  uniprotid=unlist(strsplit(proteomicstable$T..Majority.protein.IDs[i],split=";",fixed=T))
  gname=rep(NA,length(uniprotid))
  for (j in 1:length(uniprotid)){
    index=which(unirat$Entry==uniprotid[j])
    if (length(index)==1){
      gname[j]=unirat$Gene.Names..primary.[index]  
    }
  }
  genenames[[i]]=gname
  gname=gname[!is.na(gname)]
  if (length(gname)>0){
    uniname=paste(unique(gname),sep=";")    
    proteomicstable$gene_name[i]=uniname
  }
  
}

for (i in 1:nrow(proteomicstable)){
  Aempty=(proteomicstable$T..Gene.names[i]=="")
  Bempty=(proteomicstable$gene_name[i]=="")|(proteomicstable$gene_name[i]=="-")
  if (Aempty & Bempty){
    proteomicstable$gene_name[i]=paste(i)
  } else {
    if (Bempty){
      proteomicstable$gene_name[i]=proteomicstable$T..Gene.names[i]
    }  else {
      if (!Aempty){
        proteomicstable$gene_name[i]=paste(c(proteomicstable$T..Gene.names[i],proteomicstable$gene_name[i]),sep=";")
      }
    }
  }
}

ptab=proteomicstable[proteomicstable$C..Reverse!="+",]
ptab=ptab[ptab$C..Potential.contaminant!="+",]
ptab=ptab[ptab$N..Unique.peptides>1,]
lfq=ptab[,1:6]

colSums(lfq)
colSums(lfq==0)

colSums(lfq)/colSums(lfq!=0)

lfqna=lfq
lfqna[lfqna==0]=NA

summary(lfqna)

loglfqna=log2(lfqna)
summary(loglfqna)

nacount=data.frame(ko=rowSums(is.na(lfqna[,1:3])),wt=rowSums(is.na(lfqna[,4:6])))

lfqmedians=apply(loglfqna,2,median,na.rm=T)
lfqmediansrobust=apply(loglfqna[rowSums(nacount)==0,],2,median,na.rm=T)

apply(loglfqna[rowSums(nacount)==0,],2,sd,na.rm=T)


lfqnorm=loglfqna
for (i in 1:6){
  lfqnorm[,i]=loglfqna[,i]-lfqmediansrobust[i]  
}
summary(lfqnorm[rowSums(nacount)==0,])


m1=rowMeans(lfqnorm[,1:3],na.rm=T)
m2=rowMeans(lfqnorm[,4:6],na.rm=T)

mean_na=data.frame(means=c(m1,m2),nas=c(nacount[,1],nacount[,2]))

mean(mean_na$means[mean_na$nas==0],na.rm=T)
mean(mean_na$means[mean_na$nas==1],na.rm=T)
mean(mean_na$means[mean_na$nas==2],na.rm=T)

min(mean_na$means[mean_na$nas==0],na.rm=T)
min(mean_na$means[mean_na$nas==1],na.rm=T)
min(mean_na$means[mean_na$nas==2],na.rm=T)

boxplot(m~na,data=data.frame(m=rowMeans(lfqnorm[,1:3],na.rm=T),na=nacount$ko))
boxplot(m~na,data=data.frame(na=nacount$wt,m=rowMeans(lfqnorm[,4:6],na.rm=T)))


tokeep=!(nacount$ko>1 & nacount$wt>1)                  

lfqnafilt=lfqnorm[tokeep,]                   
ptabfilt=ptab[tokeep,]
nacountfilt=nacount[tokeep,]


summary(lfqnafilt)

min(lfqnafilt,na.rm=T)

#study deviation from mean across lfq scale

deviations=vector()
meanvec=vector()
k=1
for (i in 1:nrow(lfqnafilt)){
  if (nacountfilt$ko[i]==0){
    meanvalue=mean(as.numeric(lfqnafilt[i,1:3]),na.rm=T)
    meanvec[k]=meanvalue
    deviations[k]=lfqnafilt[i,1]-meanvalue
    k=k+1
    meanvec[k]=meanvalue
    deviations[k]=lfqnafilt[i,2]-meanvalue
    k=k+1
    meanvec[k]=meanvalue
    deviations[k]=lfqnafilt[i,3]-meanvalue
    k=k+1
  } else {
    if (nacountfilt$ko[i]==1){
      valvec=lfqnafilt[i,1:3]
      valvec=valvec[!is.na(valvec)]
      meanvalue=mean(as.numeric(valvec),na.rm=T)
      meanvec[k]=meanvalue
      deviations[k]=valvec[1]-meanvalue
      k=k+1
      meanvec[k]=meanvalue
      deviations[k]=valvec[2]-meanvalue
      k=k+1
    }
  }
  
  if (nacountfilt$wt[i]==0){
    meanvalue=mean(as.numeric(lfqnafilt[i,4:6]),na.rm=T)
    meanvec[k]=meanvalue
    deviations[k]=lfqnafilt[i,4]-meanvalue
    k=k+1
    meanvec[k]=meanvalue
    deviations[k]=lfqnafilt[i,5]-meanvalue
    k=k+1
    meanvec[k]=meanvalue
    deviations[k]=lfqnafilt[i,6]-meanvalue
    k=k+1
  } else {
    if (nacountfilt$wt[i]==1){
      valvec=lfqnafilt[i,4:6]
      valvec=valvec[!is.na(valvec)]
      meanvalue=mean(as.numeric(valvec),na.rm=T)
      meanvec[k]=meanvalue
      deviations[k]=valvec[1]-meanvalue
      k=k+1
      meanvec[k]=meanvalue
      deviations[k]=valvec[2]-meanvalue
      k=k+1
    }
  }
}

min(meanvec)
max(meanvec)
hist(meanvec,50)


globalsd=sqrt(sum(deviations^2)/(length(deviations)-1))
minlfq=min(lfqnafilt,na.rm=T)

deviationsord=deviations[order(meanvec)]
meanvecord=meanvec[order(meanvec)]

sdvec=vector()
for (i in 1:6834){
  sdvec[i]=sqrt(sum(deviationsord[i:(i+499)]^2)/499)
}

plot(sdvec)

#missing value inputation
# if 3 NA - input random number average min-1 and global sd
# if 2 NA - same as for 3 NA
# if 1 NA - knn

knnlookup_ko=lfqnafilt[nacountfilt$ko==0,1:3]
knnlookup_wt=lfqnafilt[nacountfilt$wt==0,4:6]
names(knnlookup_ko)=c("A","B","C")
names(knnlookup_wt)=c("A","B","C")

knnlookup=rbind(knnlookup_ko,knnlookup_wt)

findknn=function(seed,lookuptab,nneib){
  seedmat=matrix(seed,nrow=nrow(lookuptab),ncol=ncol(lookuptab),byrow = T)
  dist2seed=rowSums((seedmat-lookuptab)^2,na.rm=T)
  lookupord=lookuptab[order(dist2seed),]
  knnavg=colMeans(lookupord[1:nneib,])
  seed[is.na(seed)]=knnavg[is.na(seed)]
  seed
}

lfqpreimput=lfqnafilt


miximpute=function(lfqmat,nacountmat,minlfq,globalsd,nneib,knnlookup){
  
  for (i in 1:nrow(lfqmat)){
    if (nacountmat$ko[i]==3){
      lfqmat[i,1:3]=rnorm(3,mean=minlfq,sd=globalsd)
    } else {
      if (nacountmat$ko[i]==2){
        rval=rnorm(2,mean=minlfq,sd=globalsd)
        lfqmat[i,which(is.na(lfqmat[i,1:3]))]=rval
      } else {
        if (nacountmat$ko[i]==1){
          lfqmat[i,1:3]=findknn(as.numeric(lfqmat[i,1:3]),knnlookup,5)
        }
      }
    }
    if (nacountmat$wt[i]==3){
      lfqmat[i,4:6]=rnorm(3,mean=minlfq,sd=globalsd)
    } else {
      if (nacountmat$wt[i]==2){
        rval=rnorm(2,mean=minlfq,sd=globalsd)
        lfqmat[i,3+which(is.na(lfqmat[i,4:6]))]=rval
      } else {
        if (nacountmat$wt[i]==1){
          lfqmat[i,4:6]=findknn(as.numeric(lfqmat[i,4:6]),knnlookup,5)
        }
      }
    }
    
  }
  lfqmat
}

lfqfinal=list()

for (i in 1:30){ 
  lfqfinal[[i]]=miximpute(lfqnafilt,nacountfilt,minlfq,globalsd,5,knnlookup)
}

#limma analysis
library(limma)

toptablelist=list()

for (i in 1:30){
  protein.matrix = (as.matrix(lfqfinal[[i]]))
  row.names(protein.matrix)=ptabfilt$T..id
  
  class = as.factor(c("ko","ko","ko","wt","wt","wt"))
  design = model.matrix(~0+class) # fitting without intercept
  
  fit1 = lmFit(protein.matrix,design = design)
  cont <- makeContrasts(classko-classwt, levels = design)
  fit2 = contrasts.fit(fit1,contrasts = cont)
  fit3 <- eBayes(fit2)
  
  toptablelist[[i]]=topTable(fit3,coef=1,number=Inf)

}

toptable10=toptablelist[[1]]
toptable10$nsig=0
for (i in 1:30){
  toptable10$nsig=toptable10$nsig+(toptablelist[[i]]$adj.P.Val<0.1 & abs(toptablelist[[i]]$logFC)>1)
}

ids=row.names(toptable10)
toptable10$na_ko=0
toptable10$na_wt=0
toptable10$all_ids=""
toptable10$gene_name=""
for (j in 1:nrow(toptable10)){
  toptable10$na_ko[j]=nacountfilt$ko[ptabfilt$T..id==ids[j]]
  toptable10$na_wt[j]=nacountfilt$wt[ptabfilt$T..id==ids[j]]
  toptable10$all_ids[j]=ptabfilt$T..Majority.protein.IDs[ptabfilt$T..id==ids[j]]
  toptable10$gene_name[j]=ptabfilt$gene_name[ptabfilt$T..id==ids[j]]
}

sum(toptable10$nsig>15)

volcanoplot(fit3,coef=1, style = "p-value", highlight = 168,
            names=rownames(fit3$coefficients))


toptable10$entrez=""
toptable10$gobp=""
toptable10$gene=""
for (i in 1:nrow(toptable10)){
  gname=unlist(strsplit(toptable10$gene_name[i],split=";",fixed=T))
  if (length(gname)==1){
    toptable10$gene[i]=gname
    index=which(unirat$Gene.Names..primary.==gname)
    if (length(index)>0){
      toptable10$entrez[i]=unirat$entrez[index[1]]
      toptable10$gobp[i]=unirat$Gene.Ontology..biological.process.[index[1]]
    }
  } else {
    entrezvec=vector(length=length(gname))
    govec=vector(length=length(gname))
    ngovec=vector(length=length(gname))
    for (j in 1:length(gname)){
      index=index=which(unirat$Gene.Names..primary.==gname[j])
      if (length(index)>0){
        entrezvec[j]=unirat$entrez[index[1]]
        govec[j]=unirat$Gene.Ontology..biological.process.[index[1]]
        ngovec[j]=length(unlist(strsplit(govec[j],split=";",fixed=T)))
      }
      selectname=which.max(ngovec)
      if (length(selectname)>0){
        toptable10$gene[i]=gname[selectname]
        toptable10$entrez[i]=entrezvec[selectname]
        toptable10$gobp[i]=govec[selectname]
      }
    }
  }
}

toptable10$gene[toptable10$gene_name=="LOC108348048"]="Vwa5a"
toptable10$entrez[toptable10$gene_name=="LOC108348048"]="301097"
toptable10$gene[toptable10$all_ids=="D3ZCH7;G3V9D7"]="Add3"
toptable10$entrez[toptable10$all_ids=="D3ZCH7;G3V9D7"]="25230"
toptable10$entrez[toptable10$all_ids=="Q498U0"]="691807"
toptable10$gene[toptable10$all_ids=="Q498U0"]="C2h4orf3"
toptable10$entrez[toptable10$all_ids=="Q3KRE0;A0A8I6AD85"]="298682"
toptable10$entrez[toptable10$gene_name=="ENSRNOG00000063418"]="363248"
toptable10$gene[toptable10$gene_name=="ENSRNOG00000063418"]="Rpl37a"
toptable10$gene[toptable10$gene_name=="Txn"]="Txn1"
toptable10$entrez[toptable10$gene_name=="Txn"]="116484"
toptable10$gene[toptable10$gene_name=="Cox6c2"]="Cox6c"
toptable10$entrez[toptable10$gene_name=="Cox6c2"]="54322"
toptable10$entrez[toptable10$all_ids=="P16391"]="100188935"
toptable10$gene[toptable10$all_ids=="P16391"]="RT1-A"
toptable10$entrez[toptable10$gene_name=="Idh3B"]="94173"
toptable10$entrez[toptable10$gene_name=="Rab1A"]="81754"
toptable10$entrez[toptable10$gene_name=="ENSRNOG00000066953"]="79449"
toptable10$gene[toptable10$gene_name=="ENSRNOG00000066953"]="Rpl21"
toptable10$entrez[toptable10$gene_name=="AABR07021955.1"]="89789"
toptable10$gene[toptable10$gene_name=="AABR07021955.1"]="Lbr"
toptable10$entrez[toptable10$gene_name=="ENSRNOG00000062895"]="81765"
toptable10$gene[toptable10$gene_name=="ENSRNOG00000062895"]="Rpl13"
toptable10$entrez[toptable10$gene_name=="H4c2; Hist1h4m; H4c16"]="64627"
toptable10$gene[toptable10$gene_name=="H4c2; Hist1h4m; H4c16"]="Hist1h4b"
toptable10$entrez[toptable10$gene_name=="ENSRNOG00000065195"]="191572"
toptable10$gene[toptable10$gene_name=="ENSRNOG00000065195"]="Ahnak"


write.csv(toptable10,file="toptable10.csv",row.names = F)

