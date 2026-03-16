# MetaCCI

**MetaCCI: A Gene Regulatory Network-Based Framework for Cell--Cell
Communication Inference**

------------------------------------------------------------------------

# Overview

This repository provides the implementation of **MetaCCI**, a
computational framework for inferring **cell--cell communication (CCI)**
by integrating **gene regulatory networks (GRNs)** with ligand--receptor
gene expression.

This repository provides:

• Full analysis pipeline\
• Example dataset\
• Step‑by‑step tutorial for reproducing MetaCCI results

The MetaCCI workflow consists of two major stages:

1.  **Gene regulatory network estimation**
2.  **Cell--cell communication inference**

------------------------------------------------------------------------

# Software Requirements

The implementation is written in **R**.

Required packages:

``` r
install.packages("HDeconometrics")
install.packages("igraph")
install.packages("dplyr")
install.packages("Hmisc")
```

Load packages:

``` r
library(HDeconometrics)
library(igraph)
library(dplyr)
library(Hmisc)
```

------------------------------------------------------------------------

# MetaCCI Workflow

    Gene expression data
            ↓
    Gene regulatory network estimation
            ↓
    Identification of differentially regulated genes
            ↓
    Eigen-cell estimation
            ↓
    Cell–cell communication inference

------------------------------------------------------------------------

# Step 1: Gene Regulatory Network Estimation

Gene regulatory networks are estimated using **lasso regression with BIC
model selection**.

### Input

    MetaCCI_DATA/T_LR_EXP.csv
    MetaCCI_DATA/T_OBS300.txt

-   **T_LR_EXP.csv** : ligand--receptor gene expression matrix\
-   **T_OBS300.txt** : indices of cells belonging to different cell
    groups

### Example Code

``` r
library(HDeconometrics)

LR_EXP <- read.table("MetaCCI_DATA/T_LR_EXP.csv", sep=",")
ALL_OBS <- read.table("MetaCCI_DATA/T_OBS300.txt", sep="\t")

EXP_C <- LR_EXP[as.matrix(ALL_OBS["DS_OBS1",]), ]
EXP_Q <- LR_EXP[as.matrix(ALL_OBS["DS_OBS2",]), ]

BETA_C <- matrix(0,ncol(EXP_C),ncol(EXP_C))
colnames(BETA_C) <- colnames(EXP_C)
rownames(BETA_C) <- colnames(EXP_C)

BETA_Q <- matrix(0,ncol(EXP_Q),ncol(EXP_Q))
colnames(BETA_Q) <- colnames(EXP_Q)
rownames(BETA_Q) <- colnames(EXP_Q)

for (c in 1:ncol(LR_EXP)){

 if (c <= ncol(EXP_C) && sum(EXP_C[,c] != 0) > 0){

  lasso <- ic.glmnet(
   scale(EXP_C[,-c],FALSE,FALSE),
   scale(EXP_C[,c],FALSE,FALSE),
   crit="bic"
  )

  BETA_C[names(coef(lasso)[-1]),c] <- coef(lasso)[-1]

 }

 if (c <= ncol(EXP_Q) && sum(EXP_Q[,c] != 0) > 0){

  lasso <- ic.glmnet(
   scale(EXP_Q[,-c],FALSE,FALSE),
   scale(EXP_Q[,c],FALSE,FALSE),
   crit="bic"
  )

  BETA_Q[names(coef(lasso)[-1]),c] <- coef(lasso)[-1]

 }

}
```

------------------------------------------------------------------------

# Step 2: Construct Gene Regulatory Networks

``` r
NP_C <- c(0,0,0)
NP_Q <- c(0,0,0)

T_COL <- max(c(ncol(BETA_C),ncol(BETA_Q)))

for (c in 1:T_COL){

 if (ncol(BETA_C)>=c && sum(BETA_C[,c]!=0)>0){

  NP_C <- rbind(
   NP_C,
   cbind(
    rownames(BETA_C)[BETA_C[,c]!=0],
    colnames(BETA_C)[c],
    BETA_C[,c][BETA_C[,c]!=0]
   )
  )

 }

 if (ncol(BETA_Q)>=c && sum(BETA_Q[,c]!=0)>0){

  NP_Q <- rbind(
   NP_Q,
   cbind(
    rownames(BETA_Q)[BETA_Q[,c]!=0],
    colnames(BETA_Q)[c],
    BETA_Q[,c][BETA_Q[,c]!=0]
   )
  )

 }

}
```

------------------------------------------------------------------------

# Step 3: Compute Regulatory Effect

``` r
NP_C <- NP_C[-1,]
colnames(NP_C) <- c("RG","TG","COEF")

NP_Q <- NP_Q[-1,]
colnames(NP_Q) <- c("RG","TG","COEF")

NP_C <- cbind(
 NP_C,
 abs(as.numeric(NP_C[,3]) *
 apply(EXP_C[,NP_C[,1]],2,mean))
)

NP_Q <- cbind(
 NP_Q,
 abs(as.numeric(NP_Q[,3]) *
 apply(EXP_Q[,NP_Q[,1]],2,mean))
)

colnames(NP_C)[4] <- "RE"
colnames(NP_Q)[4] <- "RE"

write.table(NP_C,"MetaCCI_DATA/NP_C.csv",sep=",")
write.table(NP_Q,"MetaCCI_DATA/NP_Q.csv",sep=",")
```

------------------------------------------------------------------------

# Step 4: Identification of Differentially Regulated Genes

Differentially regulated genes are identified using activity scores
derived from:

• Betweenness centrality\
• Regulatory effect\
• Hubness\
• Jaccard distance between GRNs

These scores are combined to compute the **gene activity difference
(Δ)** between query and control networks.

The resulting genes are saved as:

    drGENEs.txt

------------------------------------------------------------------------

# Step 5: Eigen Cell Estimation

Eigen cells summarize gene expression patterns using **principal
component analysis (PCA)**.

``` r
PC_cut <- 0.95

EXP_Q  <- LR_EXP[as.matrix(ALL_OBS[1,]),rownames(drGENEs)]
EXP_T1 <- LR_EXP[as.matrix(ALL_OBS[2,]),rownames(drGENEs)]
EXP_T2 <- LR_EXP[as.matrix(ALL_OBS[3,]),rownames(drGENEs)]
EXP_T3 <- LR_EXP[as.matrix(ALL_OBS[4,]),rownames(drGENEs)]

X_Q <- as.matrix(scale(EXP_Q,TRUE,FALSE))

PCA <- prcomp(X_Q)

No.PC <- c(1:length(PCA$sdev))[cumsum(PCA$sdev)/sum(PCA$sdev)>=PC_cut][1]

EG_CELL_Q <- svd(X_Q)$v[,1:No.PC]
```

------------------------------------------------------------------------

# Step 6: Cell--Cell Communication Inference

CCI significance is evaluated using **correlation networks and
hypergeometric testing**.

``` r
library(Hmisc)

T_CELL <- cbind(EG_CELL_Q)

COR_P <- rcorr(as.matrix(T_CELL))$P

Pvalue <- mean(COR_P < 0.01, na.rm=TRUE)

FDR <- p.adjust(Pvalue,"BH")
```

------------------------------------------------------------------------

# Output

  Output     Description
  ---------- -------------------------------------------
  NP_C.csv   GRN for control group
  NP_Q.csv   GRN for query group
  Pvalue     statistical significance of inferred CCIs
  FDR        multiple-testing corrected q-values

------------------------------------------------------------------------

# Reproducibility

All scripts and example datasets are provided to ensure **full
reproducibility of the MetaCCI pipeline**.

Researchers can apply the workflow to **their own gene expression
datasets** by replacing the input files.

------------------------------------------------------------------------

# Citation

If you use this software, please cite the MetaCCI paper.
