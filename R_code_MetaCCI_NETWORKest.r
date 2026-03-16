library(HDeconometrics)


############################# Load expression levels and OBS for two cell groups (DS_OBS1 and DS_OBS2)
LR_EXP<-read.table("MetaCCI_DATA//T_LR_EXP.csv",sep=",")
ALL_OBS<-read.table("MetaCCI_DATA//T_OBS300.txt",sep="\t")

EXP_C<-LR_EXP[as.matrix(ALL_OBS["DS_OBS1",]),]
EXP_Q<-LR_EXP[as.matrix(ALL_OBS["DS_OBS2",]),]

############################# Gene network estimation by using lasso 
BETA_C<-matrix(numeric(ncol(EXP_C)*ncol(EXP_C)),ncol=ncol(EXP_C))
colnames(BETA_C)<-colnames(EXP_C);rownames(BETA_C)<-colnames(EXP_C)
BETA_Q<-matrix(numeric(ncol(EXP_Q)*ncol(EXP_Q)),ncol=ncol(EXP_Q))
colnames(BETA_Q)<-colnames(EXP_Q);rownames(BETA_Q)<-colnames(EXP_Q)
for (c in 1:ncol(LR_EXP)){
if (c<=ncol(EXP_C) && sum(EXP_C[,c]!=0)>0){
lasso=ic.glmnet(scale(EXP_C[,-c],FALSE,FALSE),scale(EXP_C[,c],FALSE,FALSE),crit = "bic")
BETA_C[names(coef(lasso)[-1]),c]<-coef(lasso)[-1]
}
if (c<=ncol(EXP_Q) && sum(EXP_Q[,c]!=0)>0){
lasso=ic.glmnet(scale(EXP_Q[,-c],FALSE,FALSE),scale(EXP_Q[,c],FALSE,FALSE),crit = "bic")
BETA_Q[names(coef(lasso)[-1]),c]<-coef(lasso)[-1]
}
}

############################# Constrcuted Network structure: Regulator genes, Target Genes, Edge weights
NP_C<-c(0,0,0)
NP_Q<-c(0,0,0)
T_COL<-max(c(ncol(BETA_C),ncol(BETA_Q)))
for (c in 1:T_COL){
if (ncol(BETA_C)>=c && sum(BETA_C[,c]!=0)>0){
NP_C<-rbind(NP_C,cbind(rownames(BETA_C)[BETA_C[,c]!=0],colnames(BETA_C)[c],BETA_C[,c][BETA_C[,c]!=0]))
}
if (ncol(BETA_Q)>=c && sum(BETA_Q[,c]!=0)>0){
NP_Q<-rbind(NP_Q,cbind(rownames(BETA_Q)[BETA_Q[,c]!=0],colnames(BETA_Q)[c],BETA_Q[,c][BETA_Q[,c]!=0]))
}
}
NP_C<-NP_C[-1,]
colnames(NP_C)<-c("RG","TG","COEF")
rownames(NP_C)<-paste(NP_C[,1],"_",NP_C[,2],sep="")

NP_Q<-NP_Q[-1,]
colnames(NP_Q)<-c("RG","TG","COEF")
rownames(NP_Q)<-paste(NP_Q[,1],"_",NP_Q[,2],sep="")


NP_C<-cbind(NP_C,abs(as.numeric(NP_C[,3])*apply(EXP_C[,NP_C[,1]],2,mean)))
NP_Q<-cbind(NP_Q,abs(as.numeric(NP_Q[,3])*apply(EXP_Q[,NP_Q[,1]],2,mean)))
colnames(NP_C)[4]<-"RE";colnames(NP_Q)[4]<-"RE";

write.table(NP_C,"MetaCCI_DATA//NP_C.csv",sep=",")
write.table(NP_Q,"MetaCCI_DATA//NP_Q.csv",sep=",")




