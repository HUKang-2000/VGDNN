library(caret)

load("Final_annotated.RData")

temp<-summary(as.factor(cv.data9$GENEINFO))
gene_name<-names(temp)

for(i in 1:dim(cv.data9)[1]){
  if(cv.data9$CLNSIG[i]=="Benign/Likely_benign" || cv.data9$CLNSIG[i]=="Likely_benign")
    cv.data9$CLNSIG[i]<-"Benign"
  
  if(cv.data9$CLNSIG[i]=="Likely_pathogenic" ||cv.data9$CLNSIG[i]=="Pathogenic/Likely_pathogenic")
    cv.data9$CLNSIG[i]<-"Pathogenic"
}

total<-matrix(ncol=3,nrow=30)
colnames(total)<-c("Benign","Pathogenic","Uncertain_significance")

for(i in 1:30){
  index<-which(cv.data9$GENEINFO==gene_name[i])
  total[i,]<-c(length(which(cv.data9$CLNSIG[index]=="Benign")),
               length(which(cv.data9$CLNSIG[index]=="Pathogenic")),
               length(which(cv.data9$CLNSIG[index]=="Uncertain_significance")))

}

rownames(total)<-gene_name

def<-function(self,...){
  self<-as.numeric(self)
  return(self[1]/self[2])
}
Pos_names<-c("CSQ.EXON","CSQ.cDNA_position","CSQ.CDS_position","CSQ.Protein_position")

for(i in 1:4){
  temp<-strsplit(cv.data9[,Pos_names[i]],"/")
  temp<-unlist(lapply(temp, def))
  cv.data9[,Pos_names[i]]<-temp
}

na_fill <-function(self,...){
  self[is.na(self)]<-0
  return(self)
}

## KOVA
cv.data9$KOVA_KOVA_AF<-as.numeric(cv.data9$KOVA_KOVA_AF)
cv.data9$KOVA_KOVA_AF<-na_fill(cv.data9$KOVA_KOVA_AF)

##KRGDB
cv.data9$KRGDB_KRGDB_AF<-as.numeric(cv.data9$KRGDB_KRGDB_AF)
cv.data9[,"naFlag_KRGDB_AF"]<-is.na(cv.data9$KRGDB_KRGDB_AF)
cv.data9$KRGDB_KRGDB_AF<-na_fill(cv.data9$KRGDB_KRGDB_AF)

##scSNV_ADA_SCORE
cv.data9$scSNV_ADA_SCORE<-as.numeric(cv.data9$scSNV_ADA_SCORE)

##scSNV_RF_SCORE
cv.data9$scSNV_RF_SCORE<-as.numeric(cv.data9$scSNV_RF_SCORE)

##dbNSFP_MutationAssessor_score
cv.data9$dbNSFP_MutationAssessor_score<-as.numeric(cv.data9$dbNSFP_MutationAssessor_score)

#gnomAD3_InbreedingCoeff
cv.data9$gnomAD3_InbreedingCoeff<-as.numeric(cv.data9$gnomAD3_InbreedingCoeff)

#gnomAD3_segdup
cv.data9$gnomAD3_segdup<-as.logical(cv.data9$gnomAD3_segdup)

##dbNSFP_SIFT_score
cv.data9$dbNSFP_SIFT_score<-as.numeric(cv.data9$dbNSFP_SIFT_score) # non-exist not reported 

##dbNSFP_SIFT4G_score
cv.data9$dbNSFP_SIFT4G_score<-as.numeric(cv.data9$dbNSFP_SIFT4G_score) # non-exist not reported 

##dbNSFP_APPRIS
cv.data9$dbNSFP_APPRIS<-as.factor(cv.data9$dbNSFP_APPRIS)
cv.data9$dbNSFP_APPRIS[is.na(cv.data9$dbNSFP_APPRIS)]<-"."

##dbNSFP_codon_degeneracy
cv.data9$dbNSFP_codon_degeneracy<-as.factor(cv.data9$dbNSFP_codon_degeneracy)

##dbNSFP_codonpos
cv.data9$dbNSFP_codonpos<-as.numeric(cv.data9$dbNSFP_codonpos)
cv.data9$dbNSFP_codonpos<-cv.data9$dbNSFP_codonpos-2

##dbNSFP_ALSPAC_AF
cv.data9$dbNSFP_ALSPAC_AF<-as.numeric(cv.data9$dbNSFP_ALSPAC_AF)
cv.data9[,"naFlag_dbNSFP_ALSPAC_AF"]<-is.na(cv.data9$dbNSFP_ALSPAC_AF)
cv.data9$dbNSFP_ALSPAC_AF<-na_fill(cv.data9$dbNSFP_ALSPAC_AF)

##dbNSFP_TWINSUK_AF
cv.data9$dbNSFP_TWINSUK_AF<-as.numeric(cv.data9$dbNSFP_TWINSUK_AF)
cv.data9[,"naFlag_dbNSFP_TWINSUK_AF"]<-is.na(cv.data9$dbNSFP_TWINSUK_AF)
cv.data9$dbNSFP_TWINSUK_AF<-na_fill(cv.data9$dbNSFP_TWINSUK_AF)

##dbNSFP_UK10K_AF
cv.data9$dbNSFP_UK10K_AF<-as.numeric(cv.data9$dbNSFP_UK10K_AF)
cv.data9$dbNSFP_UK10K_AF<-na_fill(cv.data9$dbNSFP_UK10K_AF)

##dbNSFP_bStatistic
cv.data9$dbNSFP_bStatistic<-as.numeric(cv.data9$dbNSFP_bStatistic)

##dbNSFP_GERP___RS
cv.data9$dbNSFP_GERP___RS<-as.numeric(cv.data9$dbNSFP_GERP___RS)

##dbNSFP_phastCons100way_vertebrate
cv.data9$dbNSFP_phastCons100way_vertebrate<-as.numeric(cv.data9$dbNSFP_phastCons100way_vertebrate)

##dbNSFP_phastCons30way_mammalian
cv.data9$dbNSFP_phastCons30way_mammalian<-as.numeric(cv.data9$dbNSFP_phastCons30way_mammalian)

##dbNSFP_phastCons17way_primate
cv.data9$dbNSFP_phastCons17way_primate<-as.numeric(cv.data9$dbNSFP_phastCons17way_primate)

##dbNSFP_phyloP100way_vertebrate
cv.data9$dbNSFP_phyloP100way_vertebrate<-as.numeric(cv.data9$dbNSFP_phyloP100way_vertebrate)

##dbNSFP_phyloP30way_mammalian
cv.data9$dbNSFP_phyloP30way_mammalian<-as.numeric(cv.data9$dbNSFP_phyloP30way_mammalian)

##dbNSFP_phyloP17way_primate
cv.data9$dbNSFP_phyloP17way_primate<-as.numeric(cv.data9$dbNSFP_phyloP17way_primate)

##dbNSFP_SiPhy_29way_logOdds
cv.data9$dbNSFP_SiPhy_29way_logOdds<-as.numeric(cv.data9$dbNSFP_SiPhy_29way_logOdds)

##dbNSFP_integrated_fitCons_score
cv.data9$dbNSFP_integrated_fitCons_score<-as.numeric(cv.data9$dbNSFP_integrated_fitCons_score)

##dbNSFP_Reliability_index
cv.data9$dbNSFP_Reliability_index<-as.numeric(cv.data9$dbNSFP_Reliability_index)
cv.data9$dbNSFP_Reliability_index[is.na(cv.data9$dbNSFP_Reliability_index)]<-10

##dbNSFP_LRT_Omega
cv.data9$dbNSFP_LRT_Omega<-as.numeric(cv.data9$dbNSFP_LRT_Omega)

##dbNSFP_LRT_score
cv.data9$dbNSFP_LRT_score<-as.numeric(cv.data9$dbNSFP_LRT_score)

##dbNSFP_Interpro_domain
cv.data9$dbNSFP_Interpro_domain<-as.logical(cv.data9$dbNSFP_Interpro_domain==".")

##gnomAD3_AF
cv.data9$gnomAD3_AF<-as.numeric(cv.data9$gnomAD3_AF)
cv.data9[,"naFlag_gnomAD3_AF"]<-is.na(cv.data9$gnomAD3_AF)
cv.data9<-cv.data9[,-270]
cv.data9$gnomAD3_AF<-na_fill(cv.data9$gnomAD3_AF)
gnomAD_names<-colnames(cv.data9)[63:266]

for(i in 1:length(gnomAD_names)){
  cv.data9[,gnomAD_names[i]]<-as.numeric(cv.data9[,gnomAD_names[i]])
  cv.data9[,gnomAD_names[i]]<-na_fill(cv.data9[,gnomAD_names[i]])
  if(sum(is.na(cv.data9[,gnomAD_names[i]]))!=0)
    print(paste(i,"omg"))

}

numeric_names<-c()
categoric_names<-c()
na_names<-c()

for(i in 28:269){
  if(is.numeric(cv.data9[,i]))
    numeric_names[i]<-colnames(cv.data9)[i]
  else
    categoric_names[i]<-colnames(cv.data9)[i]
  if(sum(is.na(cv.data9[,i]))!=0)
    na_names[i]<-colnames(cv.data9)[i]
}

numeric_names<-numeric_names[!is.na(numeric_names)]
categoric_names<-categoric_names[!is.na(categoric_names)]
na_names<-na_names[!is.na(na_names)]

modelSet_commTrt<-cv.data9[,c("GENEINFO","CLNSIG",numeric_names,categoric_names)]

dup<-duplicated(modelSet_commTrt)
modelSet_commTrt<-unique(modelSet_commTrt)

for(i in 1:30){
  index<-which(modelSet_commTrt$GENEINFO==gene_name[i])
  total[i,]<-c(length(which(modelSet_commTrt$CLNSIG[index]=="Benign")),
               length(which(modelSet_commTrt$CLNSIG[index]=="Pathogenic")),
               length(which(modelSet_commTrt$CLNSIG[index]=="Uncertain_significance")))
  
}

rownames(total)<-gene_name

length(which(modelSet_commTrt$GENEINFO=="MLH1"))
length(which(modelSet_commTrt$GENEINFO=="MSH2"))

MLH1_sam<-matrix(ncol=30,nrow=34)#9:25(P:B), 7:20
MSH2_sam<-matrix(ncol=30,nrow=31)#6:21(P:B), 4:14

set.seed(12)

for(i in 1:30){
  temp_MLH1<-modelSet_commTrt[which(modelSet_commTrt$GENEINFO=="MLH1"),]
  index1<-which(temp_MLH1$CLNSIG=="Benign")
  index2<-sample(1:length(index1),17)
  
  temp_MLH1_train_B<-temp_MLH1[index1[-index2],]
  temp_MLH1_test_B<-temp_MLH1[index1[index2],]
  
  index1<-which(temp_MLH1$CLNSIG=="Pathogenic")
  index2<-sample(1:length(index1),6)
  
  temp_MLH1_train_P<-temp_MLH1[index1[-index2],]
  temp_MLH1_test_P<-temp_MLH1[index1[index2],]
  
  assign(paste("train_",i,"_MLH1",sep=""),rbind(temp_MLH1_train_B,temp_MLH1_train_P))
  assign(paste("test_",i,"_MLH1",sep=""),rbind(temp_MLH1_test_B,temp_MLH1_test_P))
  
  temp_MSH2<-modelSet_commTrt[which(modelSet_commTrt$GENEINFO=="MSH2"),]
  
  index1<-which(temp_MSH2$CLNSIG=="Benign")
  index2<-sample(1:length(index1),14)
  temp_MSH2_train_B<-temp_MSH2[index1[-index2],]
  temp_MSH2_test_B<-temp_MSH2[index1[index2],]
  
  index1<-which(temp_MSH2$CLNSIG=="Pathogenic")
  
  index2<-sample(1:length(index1),4)
  
  temp_MSH2_train_P<-temp_MSH2[index1[-index2],]
  temp_MSH2_test_P<-temp_MSH2[index1[index2],]
  
  assign(paste("train_",i,"_MSH2",sep=""),rbind(temp_MSH2_train_B,temp_MSH2_train_P))
  assign(paste("test_",i,"_MSH2",sep=""),rbind(temp_MSH2_test_B,temp_MSH2_test_P))
}

train_for_MLH1<-modelSet_commTrt[which(modelSet_commTrt$GENEINFO!="MLH1"),]
train_for_MSH2<-modelSet_commTrt[which(modelSet_commTrt$GENEINFO!="MSH2"),]

VUS_MLH1<-modelSet_commTrt[which(modelSet_commTrt$GENEINFO=="MLH1"& modelSet_commTrt$CLNSIG=="Uncertain_significance" ),]
VUS_MSH2<-modelSet_commTrt[which(modelSet_commTrt$GENEINFO=="MSH2"& modelSet_commTrt$CLNSIG=="Uncertain_significance" ),]
train_for_MLH1<-rbind(train_for_MLH1,VUS_MLH1)
train_for_MSH2<-rbind(train_for_MSH2,VUS_MSH2)

for (i in 1:30) {
  d<-paste("train_",i,"_MLH1",sep="")
  t<-paste("test_",i,"_MLH1",sep="")
  
  d<-get(d)
  t<-get(t)
  
  d<-paste("train_",i,"_MSH2",sep="")
  t<-paste("test_",i,"_MSH2",sep="")
  
  d<-get(d)
  t<-get(t)
}

FIPS_median<-matrix(0,nrow = 30,ncol=length(na_names))
colnames(FIPS_median)<-na_names
FIPS_names<-na_names

## MLH1
for(i in 1:30){
  d<-paste("train_",i,"_MLH1",sep="")
  t<-paste("test_",i,"_MLH1",sep="")
  d<-get(d)
  t<-get(t)
  d<-rbind(d,train_for_MLH1)
  
  for(j in 1:length(FIPS_names)){
    train_median<-median(d[,FIPS_names[j]],na.rm=TRUE)
    d[,FIPS_names[j]]<-ifelse(is.na(d[,FIPS_names[j]]),train_median,d[,FIPS_names[j]])
    t[,FIPS_names[j]]<-ifelse(is.na(t[,FIPS_names[j]]),train_median,t[,FIPS_names[j]])
    FIPS_median[i,j]<-train_median
  }
  
  assign(paste("train_",i,"_MLH1",sep=""),d)
  assign(paste("test_",i,"_MLH1",sep=""),t)
}
## MSH2
for(i in 1:30){
  d<-paste("train_",i,"_MSH2",sep="")
  t<-paste("test_",i,"_MSH2",sep="")
  d<-get(d)
  t<-get(t)
  
  d<-rbind(d,train_for_MSH2)
  for(j in 1:length(FIPS_names)){
    FIPS_names[j]
    train_median<-median(d[,FIPS_names[j]],na.rm=TRUE)
    d[,FIPS_names[j]]<-ifelse(is.na(d[,FIPS_names[j]]),train_median,d[,FIPS_names[j]])
    t[,FIPS_names[j]]<-ifelse(is.na(t[,FIPS_names[j]]),train_median,t[,FIPS_names[j]])
    FIPS_median[i,j]<-train_median
  }
  assign(paste("train_",i,"_MSH2",sep=""),d)
  assign(paste("test_",i,"_MSH2",sep=""),t)
}

##MLH1
for (i in 1:30) {
  d<-paste("train_",i,"_MLH1",sep="")
  t<-paste("test_",i,"_MLH1",sep="")
  d<-get(d)
  t<-get(t)
  
  d[,"CLNSIG"]<-as.character(d[,"CLNSIG"])
  t[,"CLNSIG"]<-as.character(t[,"CLNSIG"])
  d[,"CLNSIG"]<-as.factor(d[,"CLNSIG"])
  t[,"CLNSIG"]<-as.factor(t[,"CLNSIG"])
  
  assign(paste("train_",i,"_MLH1",sep=""),d)
  assign(paste("test_",i,"_MLH1",sep=""),t)
}

##MSH2
for (i in 1:30) {
  d<-paste("train_",i,"_MSH2",sep="")
  t<-paste("test_",i,"_MSH2",sep="")
  d<-get(d)
  t<-get(t)
  
  d[,"CLNSIG"]<-as.character(d[,"CLNSIG"])
  t[,"CLNSIG"]<-as.character(t[,"CLNSIG"])
  d[,"CLNSIG"]<-as.factor(d[,"CLNSIG"])
  t[,"CLNSIG"]<-as.factor(t[,"CLNSIG"])
  
  assign(paste("train_",i,"_MSH2",sep=""),d)
  assign(paste("test_",i,"_MSH2",sep=""),t)
}

Mean_list_MLH1<-list()
Sd_list_MLH1<-list()
i<-1

##MLH1
for(i in 1:30){
  d<-paste("train_",i,"_MLH1",sep="")
  t<-paste("test_",i,"_MLH1",sep="")
  
  temp_train<-get(d)
  temp_test<-get(t)

  train_GI<-temp_train[,"GENEINFO"]
  test_GI<-temp_test[,"GENEINFO"]
  
  temp_train<-temp_train[,-c(1)]
  temp_test<-temp_test[,-c(1)]
  
  dmy<-dummyVars(~.,data = temp_train)
  new_train_AllGene_2cls<-predict(dmy,newdata = temp_train)
  new_test_AllGene_2cls<-predict(dmy,newdata = temp_test)
  
  y_train_temp<-new_train_AllGene_2cls[,"CLNSIG.Pathogenic"]
  y_test_temp<-new_test_AllGene_2cls[,"CLNSIG.Pathogenic"]
  
  VUS_tag<-new_train_AllGene_2cls[,"CLNSIG.Uncertain_significance"]
  
  t<-which(colnames(new_train_AllGene_2cls)=="dbSNP_ASPTRUE")
  t<-c(t,which(colnames(new_train_AllGene_2cls)=="gnomAD3_segdupTRUE"))
  
  s<-which(colnames(new_train_AllGene_2cls)=="dbNSFP_Interpro_domainTRUE")
  e<-which(colnames(new_train_AllGene_2cls)=="naFlag_dbNSFP_TWINSUK_AFTRUE")
  
  temp<-seq(from=s,to=e,by=2)
  
  t<-c(t,temp)
  
  binary_factor_train_AllGene_2cls<-new_train_AllGene_2cls[,t]
  binary_factor_test_AllGene_2cls<-new_test_AllGene_2cls[,t]
  
  s<-which(colnames(new_train_AllGene_2cls)=="dbNSFP_APPRIS..")
  e<-which(colnames(new_train_AllGene_2cls)=="dbNSFP_codon_degeneracy.3")
  
  mul_factor_train_AllGene_2cls<-new_train_AllGene_2cls[,s:e]
  mul_factor_test_AllGene_2cls<-new_test_AllGene_2cls[,s:e]
  
  s<-which(colnames(new_train_AllGene_2cls)=="KOVA_KOVA_AF")
  e<-which(colnames(new_train_AllGene_2cls)=="gnomAD3_AF_sas_XY")
  
  numeric_train_AllGene_2cls<-new_train_AllGene_2cls[,s:e]
  numeric_test_AllGene_2cls<-new_test_AllGene_2cls[,s:e]
  
  numeric_train_AllGene_2cls<-scale(numeric_train_AllGene_2cls)
  Mean_list_MLH1[[i]]<-attr(numeric_train_AllGene_2cls,"scaled:center")
  Sd_list_MLH1[[i]]<-attr(numeric_train_AllGene_2cls,"scaled:scale")
  numeric_test_AllGene_2cls<-scale(numeric_test_AllGene_2cls,
                                          center = attr(numeric_train_AllGene_2cls,"scaled:center"),
                                          scale = attr(numeric_train_AllGene_2cls,"scaled:scale")
                                          )
  print(i)
  
  assign(paste("x_train_",i,"_MLH1",sep=""),cbind(train_GI,binary_factor_train_AllGene_2cls,mul_factor_train_AllGene_2cls,numeric_train_AllGene_2cls)[!as.logical(VUS_tag),])
  assign(paste("x_test_",i,"_MLH1",sep=""),cbind(test_GI,binary_factor_test_AllGene_2cls,mul_factor_test_AllGene_2cls,numeric_test_AllGene_2cls))
  assign(paste("x_pre_train_",i,"_MLH1",sep=""),cbind(train_GI,binary_factor_train_AllGene_2cls,mul_factor_train_AllGene_2cls,numeric_train_AllGene_2cls)[as.logical(VUS_tag),])
  
  assign(paste("y_train_",i,"_MLH1",sep=""),y_train_temp[!as.logical(VUS_tag)])
  assign(paste("y_test_",i,"_MLH1",sep=""),y_test_temp)
}


for (i in 1:30) {
  d<-paste("x_train_",i,"_MLH1",sep="")
  t<-paste("x_test_",i,"_MLH1",sep="")
  v<-paste("x_pre_train_",i,"_MLH1",sep="")
  
  d<-get(d)
  t<-get(t)
  v<-get(v)
  abs_col<-c()
  
  for(j in 1:249){
    if(length(unique(d[,j]))==1)
      abs_col<-c(abs_col,j)
  }
  
  if(length(abs_col)!=0){
    d<-d[,-abs_col]
    t<-t[,-abs_col]
    v<-v[,-abs_col]
  }
  assign(paste("x_train_",i,"_MLH1",sep=""),d)
  assign(paste("x_test_",i,"_MLH1",sep=""),t)
  assign(paste("x_pre_train_",i,"_MLH1",sep=""),v)
}

##MSH2
Mean_list_MSH2<-list()
Sd_list_MSH2<-list()

for(i in 1:30){
  d<-paste("train_",i,"_MSH2",sep="")
  t<-paste("test_",i,"_MSH2",sep="")
  
  temp_train<-get(d)
  temp_test<-get(t)
  
  train_GI<-temp_train[,"GENEINFO"]
  test_GI<-temp_test[,"GENEINFO"]
  
  temp_train<-temp_train[,-c(1)]
  temp_test<-temp_test[,-c(1)]
  
  dmy<-dummyVars(~.,data = temp_train)
  new_train_AllGene_2cls<-predict(dmy,newdata = temp_train)
  new_test_AllGene_2cls<-predict(dmy,newdata = temp_test)

  y_train_temp<-new_train_AllGene_2cls[,"CLNSIG.Pathogenic"]
  y_test_temp<-new_test_AllGene_2cls[,"CLNSIG.Pathogenic"]
  
  VUS_tag<-new_train_AllGene_2cls[,"CLNSIG.Uncertain_significance"]
  
  t<-which(colnames(new_train_AllGene_2cls)=="dbSNP_ASPTRUE")
  t<-c(t,which(colnames(new_train_AllGene_2cls)=="gnomAD3_segdupTRUE"))
  
  s<-which(colnames(new_train_AllGene_2cls)=="dbNSFP_Interpro_domainTRUE")
  e<-which(colnames(new_train_AllGene_2cls)=="naFlag_dbNSFP_TWINSUK_AFTRUE")
  
  temp<-seq(from=s,to=e,by=2)
  
  t<-c(t,temp)
  
  binary_factor_train_AllGene_2cls<-new_train_AllGene_2cls[,t]
  binary_factor_test_AllGene_2cls<-new_test_AllGene_2cls[,t]
  
  s<-which(colnames(new_train_AllGene_2cls)=="dbNSFP_APPRIS..")
  e<-which(colnames(new_train_AllGene_2cls)=="dbNSFP_codon_degeneracy.3")
  
  mul_factor_train_AllGene_2cls<-new_train_AllGene_2cls[,s:e]
  mul_factor_test_AllGene_2cls<-new_test_AllGene_2cls[,s:e]
  
  s<-which(colnames(new_train_AllGene_2cls)=="KOVA_KOVA_AF")
  e<-which(colnames(new_train_AllGene_2cls)=="gnomAD3_AF_sas_XY")
  
  numeric_train_AllGene_2cls<-new_train_AllGene_2cls[,s:e]
  numeric_test_AllGene_2cls<-new_test_AllGene_2cls[,s:e]
  
  numeric_train_AllGene_2cls<-scale(numeric_train_AllGene_2cls)
  Mean_list_MSH2[[i]]<-attr(numeric_train_AllGene_2cls,"scaled:center")
  Sd_list_MSH2[[i]]<-attr(numeric_train_AllGene_2cls,"scaled:scale")
  
  numeric_test_AllGene_2cls<-scale(numeric_test_AllGene_2cls,
                                   center = attr(numeric_train_AllGene_2cls,"scaled:center"),
                                   scale = attr(numeric_train_AllGene_2cls,"scaled:scale")
  )
  print(i)
  
  assign(paste("x_train_",i,"_MSH2",sep=""),cbind(train_GI,binary_factor_train_AllGene_2cls,mul_factor_train_AllGene_2cls,numeric_train_AllGene_2cls)[!as.logical(VUS_tag),])
  assign(paste("x_test_",i,"_MSH2",sep=""),cbind(test_GI,binary_factor_test_AllGene_2cls,mul_factor_test_AllGene_2cls,numeric_test_AllGene_2cls))
  assign(paste("x_pre_train_",i,"_MSH2",sep=""),cbind(train_GI,binary_factor_train_AllGene_2cls,mul_factor_train_AllGene_2cls,numeric_train_AllGene_2cls)[as.logical(VUS_tag),])
  
  assign(paste("y_train_",i,"_MSH2",sep=""),y_train_temp[!as.logical(VUS_tag)])
  assign(paste("y_test_",i,"_MSH2",sep=""),y_test_temp)
}

for (i in 1:30) {
  d<-paste("x_train_",i,"_MSH2",sep="")
  t<-paste("x_test_",i,"_MSH2",sep="")
  v<-paste("x_pre_train_",i,"_MSH2",sep="")
  
  d<-get(d)
  t<-get(t)
  v<-get(v)
  
  abs_col<-c()
  
  for(j in 1:249){
    if(length(unique(d[,j]))==1)
      abs_col<-c(abs_col,j)
  }
  
  if(length(abs_col)!=0){
    d<-d[,-abs_col]
    t<-t[,-abs_col]
    v<-v[,-abs_col]
  }
  assign(paste("x_train_",i,"_MSH2",sep=""),d)
  assign(paste("x_test_",i,"_MSH2",sep=""),t)
  assign(paste("x_pre_train_",i,"_MSH2",sep=""),v)
}

for(i in 1:30){
  d<-paste("x_train_",i,"_MLH1",sep="")
  t<-paste("y_train_",i,"_MLH1",sep="")
  v<-paste("x_pre_train_",i,"_MLH1",sep="")
  
  dd<-paste("x_test_",i,"_MLH1",sep="")
  tt<-paste("y_test_",i,"_MLH1",sep="")
  
  save(list = c(d,t,v,dd,tt),
       file = file.path("Data",
                        paste0(i, "-th_Transfer_setting_MLH1_DiseaseSpecific_forNN_30sets_Nafill_standardized.RData")))
}

#MSH2
for(i in 1:30){
  d<-paste("x_train_",i,"_MSH2",sep="")
  t<-paste("y_train_",i,"_MSH2",sep="")
  v<-paste("x_pre_train_",i,"_MSH2",sep="")
  
  dd<-paste("x_test_",i,"_MSH2",sep="")
  tt<-paste("y_test_",i,"_MSH2",sep="")
  
  save(list = c(d,t,v,dd,tt),
       file = file.path("Data",
                        paste0(i, "-th_Transfer_setting_MSH2_DiseaseSpecific_forNN_30sets_Nafill_standardized.RData")))
}