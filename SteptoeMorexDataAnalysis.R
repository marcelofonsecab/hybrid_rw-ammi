# Load the SteptoeMorexData.RData file into the SxM variable
SxM = readRDS("SteptoeMorexData.RData")

# Load the Methods.R script
source("Methods.R")

# Calculate the error variance of SxM
varerror_SxM = Error_Var(SxM)

# Calculate the weights of SxM using the error variance
Weights_SxM = R.Weights(SxM, varerror_SxM)

# Perform all singular value decompositions (SVDs) using SxM (as not contaminated), SxM (as contaminated), and the calculated weights
allsvds = All_SVDS(SxM, SxM, Weights_SxM, rSVD = TRUE)

# Perform a normality test on the residuals
shapiro.test(allsvds$AMMI$residual)

# Perform a homogeneity of variances test on yield grouped by gen
fligner.test(yield ~ gen, SxM)

# Perform a homogeneity of variances test on yield grouped by env
fligner.test(yield ~ env, SxM)


# Best genotype for each environment
qtd = 1 # Change to 2 or 3 to get the best 2 and 3, respectively
g_n = as.character(unique(SxM$gen))
e_n = as.character(unique(SxM$env))

matrixofbests = cbind(
  # Calculate the best genotypes for each environment using the best_envgen function
  best_envgen(allsvds$AMMI, quant = qtd, g_n, e_n),
  
  best_envgen(allsvds$WAMMI$WAMMI_lmm.Weights.Env, quant = qtd, g_n, e_n),
  best_envgen(allsvds$WAMMI$WAMMI_lmm.Weights.Gen, quant = qtd, g_n, e_n),
  best_envgen(allsvds$WAMMI$WAMMI_lmm.Weights.GendotEnv, quant = qtd, g_n, e_n),
  best_envgen(allsvds$WAMMI$WAMMI_lmm.Weights.GenplusEnv, quant = qtd, g_n, e_n),
  best_envgen(allsvds$WAMMI$WAMMI_Weights.Rlm, quant = qtd, g_n, e_n),
  best_envgen(allsvds$WAMMI$WAMMI_lmm.Weights.RlmdotEnv, quant = qtd, g_n, e_n),
  best_envgen(allsvds$WAMMI$WAMMI_lmm.Weights.RlmdotGen, quant = qtd, g_n, e_n),
  best_envgen(allsvds$WAMMI$WAMMI_lmm.Weights.RlmdotEnv, quant = qtd, g_n, e_n),
  best_envgen(allsvds$WAMMI$WAMMI_lmm.Weights.RlmdotGenplusEnv, quant = qtd, g_n, e_n),
  
  best_envgen(allsvds$RAMMI$RAMMI_rlmm.Weights.Env, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RAMMI$RAMMI_rlmm.Weights.Gen, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RAMMI$RAMMI_rlmm.Weights.GendotEnv, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RAMMI$RAMMI_rlmm.Weights.GenplusEnv, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RAMMI$RAMMI_Weights.Rlm, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RAMMI$RAMMI_rlmm.Weights.RlmdotEnv, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RAMMI$RAMMI_rlmm.Weights.RlmdotGen, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RAMMI$RAMMI_rlmm.Weights.RlmdotEnv, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RAMMI$RAMMI_rlmm.Weights.RlmdotGenplusEnv, quant = qtd, g_n, e_n),

  best_envgen(allsvds$RWAMMI$RWAMMI_rlmm.Weights.Env, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RWAMMI$RWAMMI_rlmm.Weights.Gen, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RWAMMI$RWAMMI_rlmm.Weights.GendotEnv, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RWAMMI$RWAMMI_rlmm.Weights.GenplusEnv, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RWAMMI$RWAMMI_Weights.Rlm, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RWAMMI$RWAMMI_rlmm.Weights.RlmdotEnv, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RWAMMI$RWAMMI_rlmm.Weights.RlmdotGen, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RWAMMI$RWAMMI_rlmm.Weights.RlmdotEnv, quant = qtd, g_n, e_n),
  best_envgen(allsvds$RWAMMI$RWAMMI_rlmm.Weights.RlmdotGenplusEnv, quant = qtd, g_n, e_n))

# Final table of bests genotypes for each environment
t(matrixofbests)


# Biplots from the paper
## AMMI
.envnames = colnames(allsvds$AMMI$estimated.data)

BiplotCreation(allsvds$AMMI$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))

## W-AMMI
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.Env$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$WAMMI$WAMMI_Weights.Rlm$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.RlmdotEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))

## R-AMMI
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.Env$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RAMMI$RAMMI_Weights.Rlm$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.RlmdotEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))

## RW-AMMI
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.Env$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RWAMMI$RWAMMI_Weights.Rlm$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.RlmdotEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))

# Biplots from Supplementary material
## W-AMMI
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.Gen$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.GendotEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.GenplusEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.RlmdotGen$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.RlmdotGendotEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$WAMMI$WAMMI_lmm.Weights.GenplusEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))

## R-AMMI
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.Gen$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.GendotEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.GenplusEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.RlmdotGen$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.RlmdotGendotEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RAMMI$RAMMI_rlmm.Weights.GenplusEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))

## RW-AMMI
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.Gen$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.GendotEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.GenplusEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.RlmdotGen$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.RlmdotGendotEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))
BiplotCreation(allsvds$RWAMMI$RWAMMI_rlmm.Weights.GenplusEnv$SVD,
               env.names = .envnames,
               lims = TRUE,
               lims.U = c(-2.6, 3.7), lims.V = c(-2.94, 3.2))

