source("4_EmpiricalBayesFunctions.R")

# Compute BF(M_1 : M_0) between a pair of genes using g maximizing posterior
# likelihood of the data
Compute.BF.M1.M0 <- function(gene1, gene2) {
  Y <- Get.Gene.Data.As.Vector(gene1)
  X <- Construct.Design.Matrix(gene1, gene2)
  Xmodel1 <- X[,c(1,5)]
  n <- length(Y)
  model1 <- lm(Y ~ Xmodel1 + 0)
  R2M1 <- summary(model1)$r.squared
  Fstatistic <- R2M1 / ((1-R2M1)/(n-2))
  g <- max(Fstatistic-1, 0)
  cat("R2: ", R2M1, "\t F: ", Fstatistic, "\t g: ", g, "\n")
  BF <- ((1+g)^((n-2)/2)) / ((1+g*(1-R2M1))^((n-1)/2))
  return(BF)
}

# Compute BF(M_2 : M_3) between a pair of genes using g maximizing posterior
# likelihood of the data
Compute.BF.M2.M3 <- function(gene1, gene2) {
  Y <- Get.Gene.Data.As.Vector(gene1)
  X <- Construct.Design.Matrix(gene1, gene2)
  Xmodel2 <- X
  Xmodel3 <- X[,3:5]
  n <- length(Y);  p2 <- 4;  p3 <- 3
  model2 <- lm(Y ~ Xmodel2 + 0)
  model3 <- lm(Y ~ Xmodel3 + 0)
  R2M2 <- summary(model2)$r.squared
  R2M3 <- summary(model3)$r.squared
  FstatM2 <- (R2M2/p2) / ((1-R2M2)/(n-p2-1))
  FstatM3 <- (R2M3/p3) / ((1-R2M3)/(n-p3-1))
  gM2 <- max(FstatM2-1, 0)
  gM3 <- max(FstatM3-1, 0)
  BF.M2.M0 <- (1+gM2)^((n-p2-1)/2) / ((1+gM2*(1-R2M2))^((n-1)/2))
  BF.M3.M0 <- (1+gM3)^((n-p3-1)/2) / ((1+gM3*(1-R2M3))^((n-1)/2))
  BF.M2.M3 <- BF.M2.M0 / BF.M3.M0
  return(BF.M2.M3)
}

# For each pair of genes, compute BF(M_1 : M_0) and BF(M_LL : M_3) and 
# store these in matrices
Compute.BF.Matrices <- function() {
  
}

gene1 <- "IM4";  gene2 <- "IM2"
Plot.Gene.Pair(gene1, gene2)
Compute.BF.M1.M0(gene1, gene2)
Compute.BF.M2.M3(gene1, gene2)
