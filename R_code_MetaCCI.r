library(igraph)
require(dplyr)
library("Hmisc")

####################################### Screen differentially regylated genes ######################################
####################################################################################################################
####################################################################################################################

# Input 1 
# Ligand-Receptor genes
LR_GENEs<-read.table("MetaCCI_Data\\T_LR_GENEs.csv",sep=",")[,1]

# Average of expression levels of ligand-receptor genes
mEXP<-read.table("MetaCCI_DATA//mEXP.txt",sep="\t")
colnames(mEXP)<-gsub("[.]", "-", colnames(mEXP))
# Gene network of query group Q (colunes: regulator gene, target gene, edge weight, regulatory effect, absolute edge weigth)
NW_Q<-read.table("MetaCCI_DATA//NW_Q.csv",sep=",")
# Gene network of control group C (colunes: regulator gene, target gene, edge weight, regulatory effect, absolute edge weigth)
NW_C<-read.table("MetaCCI_DATA//NW_C.csv",sep=",")

############################## 1. Compute activity score of genes in query Q: A_Q
##### 1.1 Compute Beween centrality
g <- graph_from_edgelist(as.matrix(NW_Q[,1:2]), directed = TRUE)
BC_Q <- betweenness(g, normalized = TRUE)

##### 1.2 Compute regulatroy effect of genes
AT_Q<- NW_Q[,-2]%>% 
     group_by(RG) %>% 
     summarise_all(funs(sum))
AT_Q<-data.frame(AT_Q[,c("RG","absCOEF")])
AT_Q<-cbind(AT_Q,t(AT_Q[,2]*mEXP[1,AT_Q[,"RG"]]),BC_Q[AT_Q[,"RG"]],0,0,0)
colnames(AT_Q)[3:7]<-c("RE","BC","HB","JD","Act")

##### 1.3 Compute Hubness and Jaccard distance
for (r in 1:nrow(AT_Q)){
QR<-unique(c(NW_Q[NW_Q[,1]==rownames(AT_Q)[r],2],NW_Q[NW_Q[,2]==rownames(AT_Q)[r],1]))
VS<-unique(c(NW_C[NW_C[,1]==rownames(AT_Q)[r],2],NW_C[NW_C[,2]==rownames(AT_Q)[r],1]))
AT_Q[r,"HB"]<-length(QR)/length(unique(c(QR,VS)))
AT_Q[r,"JD"]<-(1-length(unique(intersect(QR,VS)))/length(unique(c(QR,VS))))
}
AT_Q[,"Act"]<-AT_Q[,"RE"]*((AT_Q[,"BC"]+AT_Q[,"HB"])/2)

############################## 2. Compute activity score of genes in control C: A_C
##### 2.1 Compute Beween centrality
g <- graph_from_edgelist(as.matrix(NW_C[,1:2]), directed = TRUE)
BC_C <- betweenness(g, normalized = TRUE)

##### 2.2 Compute regulatroy effect of genes
AT_C_Q<- NW_C[,-2]%>% 
     group_by(RG) %>% 
     summarise_all(funs(sum))
AT_C_Q<-data.frame(AT_C_Q[,c("RG","absCOEF")])
AT_C_Q<-cbind(AT_C_Q,t(AT_C_Q[,2]*mEXP[5,AT_C_Q[,"RG"]]),BC_C[AT_C_Q[,"RG"]],0,0,0)
colnames(AT_C_Q)[3:7]<-c("RE","BC","HB","JD","Act")

##### 2.3 Compute Hubness and Jaccard distance
for (r in 1:nrow(AT_C_Q)){
QR<-unique(c(NW_C[NW_C[,1]==rownames(AT_C_Q)[r],2],NW_C[NW_C[,2]==rownames(AT_C_Q)[r],1]))
VS<-unique(c(NW_Q[NW_Q[,1]==rownames(AT_Q)[r],2],NW_Q[NW_Q[,2]==rownames(AT_Q)[r],1]))
AT_C_Q[r,"HB"]<-length(QR)/length(unique(c(QR,VS)))
AT_C_Q[r,"JD"]<-(1-length(unique(intersect(QR,VS)))/length(unique(c(QR,VS))))
}
AT_C_Q[,"Act"]<-AT_C_Q[,"RE"]*((AT_C_Q[,"BC"]+AT_C_Q[,"HB"])/2)

############################## 3. Compute Difference of gene’s activities: Delta
DELTA_Q<-matrix(numeric(length(LR_GENEs)*1),nrow=1)
colnames(DELTA_Q)<-LR_GENEs

tR<-length(unique(c(rownames(AT_Q),rownames(AT_C_Q))))
AT_Q_mch<-data.frame(matrix(numeric(tR*ncol(AT_Q)),nrow=tR))
AT_C_Q_mch<-data.frame(matrix(numeric(tR*ncol(AT_C_Q)),nrow=tR))
rownames(AT_Q_mch)<-unique(c(rownames(AT_Q),rownames(AT_C_Q)))
rownames(AT_C_Q_mch)<-unique(c(rownames(AT_Q),rownames(AT_C_Q)))
colnames(AT_Q_mch)<-colnames(AT_Q)
colnames(AT_C_Q_mch)<-colnames(AT_C_Q)
AT_Q_mch[rownames(AT_Q),1:7]<-data.frame(AT_Q)
AT_C_Q_mch[rownames(AT_C_Q),]<-data.frame(AT_C_Q)
##### Delta
DELTA_Q[1,rownames(AT_Q_mch)]<-abs((AT_Q_mch[,"RE"]*AT_Q_mch[,"BC"]-AT_C_Q_mch[,"RE"]*AT_C_Q_mch[,"BC"]))*AT_Q_mch[,"JD"]

############################## 4. For the permuted gene networks of Q and C, above procedures were iterated and compute the difference of gene’s activities: Delta_omega
############################## 5. Compute permutation p.value and extract Differentailly regulated genes: Genes in V*_Q
# Output: drGENEs.txt


########################################## CCIs inference with Eigen cell ##########################################
####################################################################################################################
####################################################################################################################

# Cut of value for eigen cells numbers
PC_cut<-0.95
# Significance levels for correlation network
sgLv<-0.01

# Input 1 
# Expression levels of Ligand-Receptor genes
LR_EXP<-read.table("MetaCCI_Data\\T_LR_EXP.csv",sep=",")
# Ligand-Receptor pairs
LRpair<-readRDS("MetaCCI_Data\\human_lr_pair.rds")
mchID<-rbind(as.matrix(LRpair[,c("ligand_ensembl_gene_id","ligand_gene_symbol")]),as.matrix(LRpair[,c("receptor_ensembl_gene_id","receptor_gene_symbol")]))
mchID<-mchID[!duplicated(mchID), ]
colnames(LR_EXP)<-mchID[match(colnames(LR_EXP),mchID[,1]),2]

# Differentailly regulated genes
drGENEs<-read.table("MetaCCI_DATA//drGENEs.txt",sep="\t")
# Randomly selected 300 cells for query (1st row), target 1 (2nd row), target (3rd row), target 3 (4th row) groups.
ALL_OBS<-read.table("MetaCCI_DATA//T_OBS300.txt",sep="\t")

############################## 1. Eigen cells estimation : E_q, E_t1, E_t2, E_t3,
##### 1.1 Extract expression levels of differentially regulated genes for groups
EXP_Q<-LR_EXP[as.matrix(ALL_OBS[1,]),rownames(drGENEs)]
EXP_T1<-LR_EXP[as.matrix(ALL_OBS[2,]),rownames(drGENEs)]
EXP_T2<-LR_EXP[as.matrix(ALL_OBS[3,]),rownames(drGENEs)]
EXP_T3<-LR_EXP[as.matrix(ALL_OBS[1,]),rownames(drGENEs)]

X_Q<-as.matrix(scale(EXP_Q,TRUE,FALSE)) 
PCA<-prcomp(X_Q)
No.PC<-c(1:length(PCA$sdev))[cumsum(PCA$sdev)/sum(PCA$sdev)>=PC_cut][1]
EG_CELL_Q<-svd(X_Q)$v[,1:No.PC]

X_T1<-as.matrix(scale(EXP_T1,TRUE,FALSE)) 
PCA<-prcomp(X_T1)
No.PC<-c(1:length(PCA$sdev))[cumsum(PCA$sdev)/sum(PCA$sdev)>=PC_cut][1]
EG_CELL_T1<-svd(X_T1)$v[,1:No.PC]

X_T2<-as.matrix(scale(EXP_T2,TRUE,FALSE)) 
PCA<-prcomp(X_T2)
No.PC<-c(1:length(PCA$sdev))[cumsum(PCA$sdev)/sum(PCA$sdev)>=PC_cut][1]
EG_CELL_T2<-svd(X_T2)$v[,1:No.PC]

X_T3<-as.matrix(scale(EXP_T3,TRUE,FALSE)) 
PCA<-prcomp(X_T3)
No.PC<-c(1:length(PCA$sdev))[cumsum(PCA$sdev)/sum(PCA$sdev)>=PC_cut][1]
EG_CELL_T3<-svd(X_T3)$v[,1:No.PC]


##### 1.2 CCIs inference with the eigen cells
Pvalue<-matrix(numeric(10*1),nrow=1)
T_CELL<-cbind(EG_CELL_Q,EG_CELL_T1,EG_CELL_T2,EG_CELL_T3)
colnames(T_CELL)<-
c(paste("Q_",c(1:ncol(EG_CELL_Q)),sep=""),
paste("T1_",c(1:ncol(EG_CELL_T1)),sep=""),
paste("T2_",c(1:ncol(EG_CELL_T2)),sep=""),
paste("T3_",c(1:ncol(EG_CELL_T3)),sep="")
)
##### 1.2.1 CCIs between query and target 2 groups.
COR_P<-rcorr(as.matrix(T_CELL))$P
x<-sum(COR_P[substr(rownames(COR_P),1,1)=="Q",substr(colnames(COR_P),1,2)=="T1"]<sgLv,na.rm=TRUE)
M<-sum(COR_P[substr(rownames(COR_P),1,1)=="Q",substr(colnames(COR_P),1,1)!="Q"]<sgLv,na.rm=TRUE)
N<-sum(COR_P[substr(rownames(COR_P),1,1)=="Q",substr(colnames(COR_P),1,1)!="Q"]!=0,na.rm=TRUE)
n<-sum(COR_P[substr(rownames(COR_P),1,1)=="Q",substr(colnames(COR_P),1,2)=="T1"]!=0,na.rm=TRUE)
Pvalue[,1]<-1-phyper(q = x, m =M, n = N, k = n)

##### 1.2.2 CCIs between query and target 3 groups.
x<-sum(COR_P[substr(rownames(COR_P),1,1)=="Q",substr(colnames(COR_P),1,2)=="T2"]<sgLv,na.rm=TRUE)
M<-sum(COR_P[substr(rownames(COR_P),1,1)=="Q",substr(colnames(COR_P),1,1)!="Q"]<sgLv,na.rm=TRUE)
N<-sum(COR_P[substr(rownames(COR_P),1,1)=="Q",substr(colnames(COR_P),1,1)!="Q"]!=0,na.rm=TRUE)
n<-sum(COR_P[substr(rownames(COR_P),1,1)=="Q",substr(colnames(COR_P),1,2)=="T2"]!=0,na.rm=TRUE)
Pvalue[,2]<-1-phyper(q = x, m =M, n = N, k = n)

##### 1.2.3 CCIs between query and target 4 groups.
x<-sum(COR_P[substr(rownames(COR_P),1,1)=="Q",substr(colnames(COR_P),1,2)=="T3"]<sgLv,na.rm=TRUE)
M<-sum(COR_P[substr(rownames(COR_P),1,1)=="Q",substr(colnames(COR_P),1,1)!="Q"]<sgLv,na.rm=TRUE)
N<-sum(COR_P[substr(rownames(COR_P),1,1)=="Q",substr(colnames(COR_P),1,1)!="Q"]!=0,na.rm=TRUE)
n<-sum(COR_P[substr(rownames(COR_P),1,1)=="Q",substr(colnames(COR_P),1,2)=="T3"]!=0,na.rm=TRUE)
Pvalue[,3]<-1-phyper(q = x, m =M, n = N, k = n)

##### 1.2.4 CCIs between target 2 and target 3 groups.
x<-sum(COR_P[substr(rownames(COR_P),1,2)=="T1",substr(colnames(COR_P),1,2)=="T2"]<sgLv,na.rm=TRUE)
M<-sum(COR_P[substr(rownames(COR_P),1,2)=="T1",substr(colnames(COR_P),1,2)!="T1"]<sgLv,na.rm=TRUE)
N<-sum(COR_P[substr(rownames(COR_P),1,2)=="T1",substr(colnames(COR_P),1,2)!="T1"]!=0,na.rm=TRUE)
n<-sum(COR_P[substr(rownames(COR_P),1,2)=="T1",substr(colnames(COR_P),1,2)=="T2"]!=0,na.rm=TRUE)
Pvalue[,4]<-1-phyper(q = x, m =M, n = N, k = n)

##### 1.2.5 CCIs between target 2 and target 4 groups.
x<-sum(COR_P[substr(rownames(COR_P),1,2)=="T1",substr(colnames(COR_P),1,2)=="T3"]<sgLv,na.rm=TRUE)
M<-sum(COR_P[substr(rownames(COR_P),1,2)=="T1",substr(colnames(COR_P),1,2)!="T1"]<sgLv,na.rm=TRUE)
N<-sum(COR_P[substr(rownames(COR_P),1,2)=="T1",substr(colnames(COR_P),1,2)!="T1"]!=0,na.rm=TRUE)
n<-sum(COR_P[substr(rownames(COR_P),1,2)=="T1",substr(colnames(COR_P),1,2)=="T3"]!=0,na.rm=TRUE)
Pvalue[,5]<-1-phyper(q = x, m =M, n = N, k = n)


# P.value of CCIs
Pvalue<-round(Pvalue,5)
Pvalue

# FDR-q.value of CCIs
FDR<-round(p.adjust(Pvalue[,1:5], "BH"), 5)
FDR