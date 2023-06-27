## Required Packages ##
library(dplyr) # For the pipe function '%>%'
library(pbapply) # To run the function in parallel processing
library(lme4) # Linear mixed models (LMM)
library(robustlmm) # Robust linear mixed models (RLMM)
library(MASS) # Robust linear model (RLM)
library(reshape2) # 'melt' function to reshape the data for the created functions
library(rsvddpd) # Robust SVD

## Functions ##

# The sim.amb function simulates a new dataset based on the main effects and multiplicative effects.
sim.amb <- function(seed = NULL, Ngen = 100, Nenv = 8, Ncomp = 2,
                    effectGlobal = c(mean = 15, sd = sqrt(3)),
                    effectGen = c(mean = 5, sd = 1), 
                    effectEnv = c(mean = 8, sd = sqrt(2)), 
                    k = rep(1, Nenv)) {
  
  # Simulate Genotype-Environment data
  
  # Arguments:
  # - seed: Seed value for simulation (default is NULL).
  # - Ngen: Number of genotypes.
  # - Nenv: Number of environments.
  # - Ncomp: Number of principal components to use.
  # - effectGlobal: Mean and standard deviation of the global mean (default is c(mean = 15, sd = sqrt(3))).
  # - effectGen: Mean and standard deviation of genotype effect (default is c(mean = 5, sd = 1)).
  # - effectEnv: Mean and standard deviation of environment effect (default is c(mean = 8, sd = sqrt(2))).
  # - k: Vector of weights for each principal component (default is rep(1, Nenv)).
  
  seed.aux <- ifelse(is.null(seed), sample(1:2^16, 1), seed) # If no seed provided, generate a random seed
  set.seed(seed.aux)
  
  globalMean <- rnorm(1, effectGlobal[1], effectGlobal[2]) # Generate global mean
  alpha <- rnorm(Ngen, mean = effectGen[1], effectGen[2]) # Generate genotype effects
  beta <- rnorm(Nenv, mean = effectEnv[1], effectEnv[2]) # Generate environment effects
  
  rand.mat <- matrix(runif(Ngen * Nenv, min = -0.5, max = 0.5), ncol = Nenv)
  rand.svd <- svd(rand.mat) # SVD of random matrix
  rand.u <- rand.svd$u # Matrix U
  rand.v <- rand.svd$v # Matrix V
  rand.d <- diag(rand.svd$d) # Diagonal vector of matrix D
  
  simulated.amb <- matrix(rep(1, Ngen)) %*% t(matrix(rep(1, Nenv))) * globalMean +
    alpha %*% t(matrix(rep(1, Nenv))) +
    matrix(rep(1, Ngen)) %*% t(beta) # Matrix with only the effects
  # Simulate the environment using N principal components from the random matrix via U(-0.5, 0.5)
  for (j in 1:Ncomp) {
    simulated.amb <- simulated.amb + 
      ((rand.u[, j] * rand.d[j, j]) %*% t(rand.v[, j])) * (k[j]) # Add random effects
  }
  
  # Organize the matrix and transform it into a dataframe
  aux <- matrix(rep(0, Ngen * Nenv * 4), ncol = 4)
  aux <- as.data.frame(aux)
  colnames(aux) <- c("gen", "env", "rep", "yield")
  aux$gen <- rep(paste0("G", sprintf('%0.3d', 1:Ngen)), Nenv)
  aux[order(aux$gen), ]$env <- rep(paste("E", 1:Nenv, sep = ""), Ngen)
  aux$yield <- c(simulated.amb)
  
  # Create a dataframe with all replications and a dataframe of means based on replications
  finaldf <- aux
  finaldf$gen <- as.factor(finaldf$gen)
  finaldf$env <- as.factor(finaldf$env)
  finaldf$rep <- 1
  finaldf$rep <- as.factor(finaldf$rep)
  data <- as.data.frame(finaldf)
  
  return(data)
}

# The Create.Replications function generates replication of genotype-environment data.
Create.Replications = function(data, reps = 2, sig = 1) {
  
  
  # Arguments:
  # - data: Original dataset containing genotype-environment data.
  # - reps: Number of replications to create (default is 2).
  # - sig: Standard deviation of the residuals (default is 1).
  
  Nenv = length(levels(data$env))
  Ngen = length(levels(data$gen))
  N = Nenv * Ngen
  data.tmp = data.aux = NULL
  
  for (i in 1:reps) {
    resids = rnorm(N, mean = 0, sd = sqrt(sig)) # Generate random residuals
    data.aux = data
    data.aux$yield = data.aux$yield + resids # Add residuals to the yield values
    data.aux$rep = i # Assign replication number
    data.tmp = rbind(data.tmp, data.aux) # Append replicated data to the temporary dataframe
  }
  
  return(data.tmp)
}

# The Data_Contamination function contaminates the first replication of the simulated data.
Data_Contamination = function(data, percentage = 5, seed = 1,
                              type = "shift", k = 7, c = 10) {

  # Arguments:
  # - data: Original dataset containing genotype-environment data.
  # - percentage: Percentage of data to contaminate (default is 5%).
  # - seed: Seed value for reproducibility (default is 1).
  # - type: Type of contamination ("shift", "pointmass", or "varinflated", default is "shift").
  # - k: Shift factor for "shift" type contamination (default is 7).
  # - c: Scaling factor for "pointmass" and "varinflated" type contamination (default is 10).
  
  set.seed(seed)
  types_aux = c("shift", "pointmass", "varinflated")
  match.arg(type, types_aux) # Match the type argument to one of the allowed types
  percentage_aux = ifelse(percentage > 1, percentage / 100, percentage) # Convert percentage to decimal if necessary
  df_aux = data[data$rep == 1, ] # Subset the data for the first replication
  matrix_aux = tapply(data$yield,
                      data[, c("gen", "env")],
                      mean) # Calculate means for each genotype-environment combination
  Means_Deviations = cbind(mean = apply(matrix_aux, 2, mean),
                           deviation = apply(matrix_aux, 2, sd)) # Calculate means and deviations
  
  Nenv = levels(df_aux$env) %>% length() # Number of environments
  Ngen = levels(df_aux$gen) %>% length() # Number of genotypes
  Nrow = nrow(df_aux) # Number of rows in the subsetted data
  
  sample_aux = sample(1:Nrow, Nrow * percentage_aux, replace = FALSE) # Randomly sample rows to contaminate
  df_tmp = df_aux[sample_aux,] # Subset the contaminated rows
  
  for (i in 1:Nenv) {
    sample_tmp = which(df_tmp$env == rownames(Means_Deviations)[i]) # Subset rows for each environment
    df_tmp[sample_tmp, "yield"] = rnorm(length(sample_tmp),
                                        mean = Means_Deviations[i, 1] + ifelse(type == "shift", k, 0) * Means_Deviations[i, 2],
                                        sd = Means_Deviations[i, 2] / ifelse(type == "pointmass" | type == "varinflated", sqrt(c), 1)) # Contaminate yield values based on the selected type
  }
  
  df_aux[sample_aux, ] = df_tmp # Replace the contaminated rows in the original data subset
  df_aux = rbind(df_aux, data[data$rep != 1,]) # Combine the contaminated subset with the remaining data
  data.full = list(data.not_contaminated = data, data.contaminated = df_aux) # Create a list containing both the original and contaminated data
  
  return(data.full)
}

# The transform_usable_data function takes a dataframe or matrix as input and applies a specified function to transform the data.
transform_usable_data = function(dataframe, func, type = c("dataframe", "matrix")){
  type_aux = c("dataframe", "matrix")
  match.arg(type, type_aux)  # Matches the 'type' argument with the valid options
  
  if(type == "dataframe"){
    # If the 'type' argument is set to "dataframe"
    # Transform the dataframe by grouping it by 'gen' and 'env'
    # Summarize the 'yield' variable using the specified function
    # Sort the resulting dataframe by 'env'
    df = dataframe %>%
      group_by(gen, env) %>%
      summarise(yield = func(yield)) %>%
      arrange(env)
  } else if(type == "matrix"){
    # If the 'type' argument is set to "matrix"
    # Transform the matrix by applying the function 'func' to the 'yield' variable
    # Use 'gen' and 'env' as indices to define the resulting matrix
    df = tapply(dataframe[, "yield"], dataframe[, c("gen", "env")], func)
  }
  
  return(df)  # Return the transformed dataframe or matrix
}

# The LMM.Error_variance function calculates the error variance for a linear mixed model (LMM) based on the input data.
LMM.Error_variance = function(data){
  quant.gen = length(unique(data$gen))  # Calculate the number of unique genotypes
  quant.env = length(unique(data$env))  # Calculate the number of unique environments
  
  if(quant.env == 1){
    aux = lmer(yield ~ rep + (1|gen), data = data)  # Fit an LMM with random intercept for 'gen'
    tmp = attr(VarCorr(aux), "sc")^2  # Extract and square the error standard deviation component
  }
  else if(quant.gen == 1){
    aux = lmer(yield ~ rep + (1|env), data = data)  # Fit an LMM with random intercept for 'env'
    tmp = attr(VarCorr(aux), "sc")^2  # Extract and square the error standard deviation component
  }
  
  return(tmp)  # Return the error variance
}

# The RLMM.Error_variance function calculates the error variance for a robust linear mixed model (RLMM) based on the input data.
RLMM.Error_variance = function(data){
  quant.gen = length(unique(data$gen))  # Calculate the number of unique genotypes
  quant.env = length(unique(data$env))  # Calculate the number of unique environments
  
  if(quant.env == 1){
    aux = rlmer(yield ~ rep + (1|gen), data = data)  # Fit an RLMM with random intercept for 'gen'
    tmp = attr(VarCorr(aux), "sc")^2  # Extract and square the error standard deviation component
  }
  else if(quant.gen == 1){
    aux = rlmer(yield ~ rep + (1|env), data = data)  # Fit an RLMM with random intercept for 'env'
    tmp = attr(VarCorr(aux), "sc")^2  # Extract and square the error standard deviation component
  }
  
  return(tmp)  # Return the error variance
}

# The Error_Var function calculates the error variance for different groups (genotypes and environments) using either LMM or RLMM, based on the input data.
Error_Var = function(data, cluster = NULL){
  
  if(!is.null(cluster)){
    print("Using parellel loop")  # Print a message indicating parallel processing if 'cluster' is not null
  }
  
  # Calculate error variance for each genotype using LMM
  tmp1 = data %>%
    group_by(gen) %>%
    group_split() %>%
    pblapply(LMM.Error_variance, cl = cluster) %>% 
    unlist() %>%
    tibble(Group = unique(data$gen), Error_Variance = .) 
  
  # Calculate error variance for each environment using LMM
  tmp2 = data %>%
    group_by(env) %>%
    group_split() %>%
    pblapply(LMM.Error_variance, cl = cluster) %>% 
    unlist() %>% 
    tibble(Group = unique(data$env), Error_Variance = .) 
  
  # Calculate error variance for each genotype using RLMM
  tmp3 = data %>%
    group_by(gen) %>%
    group_split() %>%
    pblapply(RLMM.Error_variance, cl = cluster) %>% 
    unlist() %>% 
    tibble(Group = unique(data$gen), Error_Variance = .) 
  
  # Calculate error variance for each environment using RLMM
  tmp4 = data %>%
    group_by(env) %>%
    group_split() %>%
    pblapply(RLMM.Error_variance, cl = cluster) %>% 
    unlist() %>% 
    tibble(Group = unique(data$env), Error_Variance = .) 
  
  # Store the results in a list
  tmp = list(GenotypeLMM = tmp1, EnvironmentLMM = tmp2,
             GenotypeRLMM = tmp3, EnvironmentRLMM = tmp4)
  
  return(tmp)  # Return the list of error variances
}

# The R.Weights function calculates weights for different combinations of LMM and RLMM models based on the input data and error variances.
R.Weights = function(data, errorvariances){
  
  # Calculate the number of genotypes, environments, and repetitions
  Ngen = nlevels(data$gen)
  Nenv = nlevels(data$env)
  Nrep = max(data$rep)
  
  # Calculate weights based on repetitions for LMM and RLMM for environments
  rep.weights = tapply(data[,"rep"],
                       data[,c("gen", "env")],
                       max) / Nrep
  
  # Calculate weights for LMM models based on error variances for environments and genotypes
  lmm.env.weights.tmp = 1/errorvariances$EnvironmentLMM$Error_Variance
  lmm.Weights.env.tmp = matrix(rep(lmm.env.weights.tmp/max(lmm.env.weights.tmp), Ngen),
                               ncol = Nenv, byrow = T)
  lmm.Weights.env = lmm.Weights.env.tmp * rep.weights
  
  lmm.gen.weight.tmp = 1/errorvariances$GenotypeLMM$Error_Variance
  lmm.Weights.gen.tmp = matrix(rep(lmm.gen.weight.tmp/max(lmm.gen.weight.tmp), Nenv),
                               nrow = Ngen, byrow = F)
  lmm.Weights.gen = lmm.Weights.gen.tmp * rep.weights
  
  # Calculate weights for RLM models based on the median-transformed data
  data.rlm = transform_usable_data(data, median, "dataframe")
  model.rlm.tmp = rlm(yield ~ gen + env, data = data.rlm)
  rlm.weights.tmp = model.rlm.tmp$w
  Weights.rlm = matrix(rlm.weights.tmp, ncol = Nenv, byrow = F)
  
  # Combinations LMM
  
  # Calculate weights for different combinations of LMM models
  lmm.Weights.gendotenv = lmm.Weights.gen * lmm.Weights.env
  lmm.Weights.genplusenv = ( (lmm.Weights.gen + lmm.Weights.env) / 2 )
  lmm.Weights.rlmdotenv = Weights.rlm * lmm.Weights.env
  lmm.Weights.rlmdotgen = Weights.rlm * lmm.Weights.gen
  lmm.Weights.rlmdotgendotenv = Weights.rlm * lmm.Weights.gen * lmm.Weights.env
  lmm.Weights.rlmdotgenplusenv = Weights.rlm * ( (lmm.Weights.gen + lmm.Weights.env) / 2 )
  
  # Combinations RLMM
  
  # Calculate weights for different combinations of RLMM models
  rlmm.env.weights.tmp = 1/errorvariances$EnvironmentRLMM$Error_Variance
  rlmm.Weights.env.tmp = matrix(rep(rlmm.env.weights.tmp/max(rlmm.env.weights.tmp), Ngen),
                                ncol = Nenv, byrow = T)
  rlmm.Weights.env = rlmm.Weights.env.tmp * rep.weights
  
  rlmm.gen.weight.tmp = 1/errorvariances$GenotypeRLMM$Error_Variance
  rlmm.Weights.gen.tmp = matrix(rep(rlmm.gen.weight.tmp/max(rlmm.gen.weight.tmp), Nenv),
                                nrow = Ngen, byrow = F)
  rlmm.Weights.gen = rlmm.Weights.gen.tmp * rep.weights
  
  rlmm.Weights.gendotenv = rlmm.Weights.gen * rlmm.Weights.env
  rlmm.Weights.genplusenv = ( (rlmm.Weights.gen + rlmm.Weights.env) / 2 )
  rlmm.Weights.rlmdotenv = Weights.rlm * rlmm.Weights.env
  rlmm.Weights.rlmdotgen = Weights.rlm * rlmm.Weights.gen
  rlmm.Weights.rlmdotgendotenv = Weights.rlm * rlmm.Weights.gen * rlmm.Weights.env
  rlmm.Weights.rlmdotgenplusenv = Weights.rlm * ( (rlmm.Weights.gen + rlmm.Weights.env) / 2 )
  
  # Create a list to store the calculated weights
  weights_list <- list(
    LMM = list(
      lmm.Weights.Env = lmm.Weights.env,
      lmm.Weights.Gen = lmm.Weights.gen,
      lmm.Weights.GendotEnv = lmm.Weights.gendotenv,
      lmm.Weights.GenplusEnv = lmm.Weights.genplusenv,
      Weights.Rlm = Weights.rlm,
      lmm.Weights.RlmdotEnv = lmm.Weights.rlmdotenv,
      lmm.Weights.RlmdotGen = lmm.Weights.rlmdotgen,
      lmm.Weights.RlmdotGendotEnv = lmm.Weights.rlmdotgendotenv,
      lmm.Weights.RlmdotGenplusEnv = lmm.Weights.rlmdotgenplusenv
    ),
    RLMM = list(
      rlmm.Weights.Env = rlmm.Weights.env,
      rlmm.Weights.Gen = rlmm.Weights.gen,
      rlmm.Weights.GendotEnv = rlmm.Weights.gendotenv,
      rlmm.Weights.GenplusEnv = rlmm.Weights.genplusenv,
      Weights.Rlm = Weights.rlm,
      rlmm.Weights.RlmdotEnv = rlmm.Weights.rlmdotenv,
      rlmm.Weights.RlmdotGen = rlmm.Weights.rlmdotgen,
      rlmm.Weights.RlmdotGendotEnv = rlmm.Weights.rlmdotgendotenv,
      rlmm.Weights.RlmdotGenplusEnv = rlmm.Weights.rlmdotgenplusenv
    )
  )
  
  return(weights_list)  # Return the list of calculated weights
}

# The WeightedSvdRes function performs weighted singular value decomposition (SVD) or rSVD on the input data.
WeightedSvdRes <- function(data, weight = NULL, Ncomp = 2, robust = FALSE){
  # Extract the number of genotypes and environments from the data matrix
  Ngen = nrow(data)
  Nenv = ncol(data)
  
  # Convert the weight vector to a numeric vector
  W = c(weight)
  
  # Assign the data matrix to D
  D = data
  
  # Initialize X as an empty matrix with dimensions Ngen x Nenv
  X = matrix(0, Ngen, Nenv)
  
  # Initialize auxiliary variables
  aux = matrix(1, Ngen, Nenv)
  Xold = Inf*aux
  Err = Inf
  eps = 1e-5
  
  # Perform the initial SVD using either svd or rSVDdpd function based on the 'robust' parameter
  Xold = X
  if(!robust){
    wsvd = svd(W*D + (1-W)*X)
  } else {
    wsvd = rSVDdpd(W*D + (1-W)*X, alpha = 0.3, maxiter = 100,
                   initu= svd(W*D + (1-W)*X)$u, initv = svd(W*D + (1-W)*X)$v,
                   eps = 1e-4)
  }
  U = wsvd$u
  S = diag(wsvd$d)
  V = wsvd$v
  
  # Truncate singular values beyond Ncomp
  if(Ncomp + 1 < length(wsvd$d)){
    S[(Ncomp + 1):length(wsvd$d), (Ncomp + 1):length(wsvd$d)] = 0
  }
  
  # Calculate the updated X matrix
  X = U%*%S%*%t(V)
  
  m = 0
  while(m < 100){
    m = m + 1
    Xold = X
    
    # Perform the SVD on the updated weighted matrix using either svd or rSVDdpd function based on the 'robust' parameter
    if(!robust){
      wsvd = svd(W*D + (1-W)*X)
    } else {
      wsvd = rSVDdpd(W*D + (1-W)*X, alpha = 0.3, maxiter = 100,
                     initu= svd(W*D + (1-W)*X)$u, initv = svd(W*D + (1-W)*X)$v,
                     eps = 1e-4)
    }
    U = wsvd$u
    S = diag(wsvd$d)
    V = wsvd$v
    
    # Truncate singular values beyond Ncomp
    if(Ncomp + 1 < length(wsvd$d)){
      S[(Ncomp + 1):length(wsvd$d), (Ncomp + 1):length(wsvd$d)] = 0
    }
    
    # Calculate the updated X matrix
    X = U%*%S%*%t(V)
    
    # Calculate the difference between the updated X and the previous X
    Err = sum(sum((X-Xold)^2))
  }
  
  # Return the final X matrix
  return(X)
}

# The AMMI model
ammi.model <- function(dataframe, Ncomp = 2){
  # Calculate mean values per genotype and environment
  df.mean = transform_usable_data(dataframe, mean, type = "dataframe")
  
  # Extract the number of genotypes and environments from the mean data frame
  Ngen = nlevels(df.mean$gen)
  Nenv = nlevels(df.mean$env)
  
  # Extract the names of environments and genotypes
  env.names = levels(df.mean$env)
  gen.names = levels(df.mean$gen)
  
  # Simple AMMI model
  model.GEI = lm(yield ~ gen + env, data = df.mean)  
  matrix.GEI = matrix(residuals(model.GEI), ncol = Nenv, nrow = Ngen)
  svd.GEI = svd(matrix.GEI)
  svd.u = svd.GEI$u
  svd.d = diag(svd.GEI$d)
  svd.v = svd.GEI$v
  
  aux = matrix(model.GEI$fitted.values, ncol = Nenv)
  
  if(Ncomp >= 1){
    for(i in 1:Ncomp){
      aux = aux + (svd.u[,i]*svd.d[i,i])%*%t(svd.v[,i])
    }
  }
  
  colnames(aux) = env.names
  rownames(aux) = gen.names
  matrix.fitted = matrix(model.GEI$fitted.values, ncol = Nenv)
  
  SVD.red = 0
  
  if(Ncomp == 1){
    SVD.red = 0
  } else {
    SVD.red = svd.u[, 1:Ncomp]%*%(svd.d[1:Ncomp, 1:Ncomp])%*%t(svd.v[, 1:Ncomp])
  }
  
  colnames(SVD.red) = env.names
  rownames(SVD.red) = gen.names
  
  return(list(estimated.data = aux, anovatable = anova(model.GEI), SVD = svd.GEI, 
              residual = matrix.GEI, matrix.fitted = matrix.fitted, reduced.SVD = SVD.red))
}

# The W-AMMI model
wammi.model <- function(data, weight = NULL, Ncomp = 2){
  # Calculate the mean values per genotype and environment
  dataframe = transform_usable_data(data, mean, type = "dataframe")
  
  # Extract the number of environments from the transformed data frame
  Nenv = nlevels(dataframe$env)
  
  # Check if the weight vector is provided
  if(is.null(weight)){
    stop("This function requires a weight vector")
  }
  
  weight = c(weight)
  
  # Perform weighted linear regression
  weighted.lm = lm(yield ~ gen + env, weights = weight, data = dataframe)
  
  fitted.values = matrix(weighted.lm$fitted.values, ncol = Nenv)
  
  residual.matrix = matrix(weighted.lm$residuals, ncol = Nenv)
  
  # Set column and row names for the residual matrix
  colnames(residual.matrix) = levels(dataframe$env)
  rownames(residual.matrix) = levels(dataframe$gen)
  
  # Perform Weighted Singular Value Decomposition (W-SVD)
  SVD.aux = WeightedSvdRes(residual.matrix, weight = weight, Ncomp = Ncomp)
  
  colnames(SVD.aux) = colnames(residual.matrix)
  rownames(SVD.aux) = rownames(residual.matrix)
  
  SVD.red = SVD.aux
  SVD.fin = svd(SVD.red)
  
  # Append the weighted residuals to the transformed data frame
  dataframe = cbind(dataframe, "W.residuals" = melt(residual.matrix)$value)
  
  return(list(dataframe = dataframe,
              residual = residual.matrix,
              matrix.fitted = fitted.values,
              reduced.SVD = SVD.red,
              SVD = SVD.fin))
}

# The R-AMMI model
rammi.model <- function(dataframe, weight = NULL, Ncomp = 2){
  # Calculate the median values per genotype and environment
  data = transform_usable_data(dataframe, median, "dataframe")
  
  # Extract the number of genotypes and environments from the transformed data frame
  Ngen = nlevels(data$gen)
  Nenv = nlevels(data$env)
  
  # Convert weight to a vector
  weight = c(weight)
  
  # Check if the weight vector is provided
  if(is.null(weight)){
    # Perform robust linear regression without weights
    model = rlm(yield ~ gen + env, data = data)
  } else {
    # Perform robust linear regression with weights
    model = rlm(yield ~ gen + env, data = data, w = weight)
  }
  
  
  # Create the residual matrix from the model residuals
  residual.matrix = matrix(model$residuals, ncol = Nenv, nrow = Ngen)
  
  # Create the matrix of fitted values from the model
  fitted.values = matrix(model$fitted.values, ncol = Nenv, nrow = Ngen)
  
  # Perform robust Singular Value Decomposition (rSVD)
  SVD = rSVDdpd(residual.matrix, alpha = 0.3, maxiter = 100,
                initu= svd(residual.matrix)$u, initv = svd(residual.matrix)$v,
                eps = 1e-4)
  
  # Reduce the SVD matrix to the specified number of components
  SVD.red = SVD$u[, 1:Ncomp] %*% diag(SVD$d[1:Ncomp]) %*% t(SVD$v[, 1:Ncomp])
  
  return(list(residual = residual.matrix,
              matrix.fitted = fitted.values,
              reduced.SVD = SVD.red,
              SVD = SVD))
}

# The RW-AMMI model
rwammi.model <- function(dataframe, weight = NULL, Ncomp = 2, rSVD = FALSE){
  # Calculate the median values per genotype and environment
  data = transform_usable_data(dataframe, median, "dataframe")
  
  # Extract the number of genotypes and environments from the transformed data frame
  Ngen = nlevels(data$gen)
  Nenv = nlevels(data$env)
  
  # Convert weight to a vector
  weight = c(weight)
  
  # Perform robust linear regression with weights
  model = rlm(yield ~ gen + env, data = data, w = weight)
  
  # Create the residual matrix from the model residuals
  residual.matrix = matrix(model$residuals, ncol = Nenv, nrow = Ngen)
  
  # Create the matrix of fitted values from the model
  fitted.values = matrix(model$fitted.values, ncol = Nenv, nrow = Ngen)
  
  # Perform weighted SVD on the residual matrix
  SVD.aux = WeightedSvdRes(residual.matrix, weight = weight, Ncomp = Ncomp, robust = rSVD)
  
  # Set the column and row names of the reduced SVD matrix
  colnames(SVD.aux) = colnames(residual.matrix)
  rownames(SVD.aux) = rownames(residual.matrix)
  
  # Assign the reduced SVD matrix
  SVD.red = SVD.aux
  
  # Perform SVD on the reduced SVD matrix
  SVD.fin = svd(SVD.red)
  
  return(list(residual = residual.matrix,
              matrix.fitted = fitted.values,
              reduced.SVD = SVD.red,
              SVD = SVD.fin))
}

# The All_SVDs function stores all decompositions outputs from all models
All_SVDS = function(data_not_contaminated, data_contaminated, weights, Ncomp = 2, rSVD = F){
  # Initialize empty lists to store the models
  WAMMI = RAMMI = RWAMMI = list()
  
  # Fit AMMI model to the real data
  Estimated.real = ammi.model(data_not_contaminated, Ncomp = Ncomp)
  
  # Fit AMMI model to the complete data
  AMMI.model = ammi.model(data_contaminated, Ncomp = Ncomp)
  
  # Iterate over the weights for WAMMI, RAMMI, and RWAMMI models
  for(i in 1:length(weights$LMM)){
    print(i)
    
    # Fit WAMMI model with the current weight
    WAMMI[[paste0("WAMMI_", names(weights$LMM)[i])]] = 
      wammi.model(data_contaminated, weight = weights$LMM[[i]], Ncomp = Ncomp)
    
    # Fit RAMMI model with the current weight (use NULL weight for i = 5)
    if(i == 5){
      RAMMI[[paste0("RAMMI_", names(weights$RLMM)[i])]] = 
        rammi.model(data_contaminated, weight = NULL, Ncomp = Ncomp)
    } else {
      RAMMI[[paste0("RAMMI_", names(weights$RLMM)[i])]] = 
        rammi.model(data_contaminated, weight = weights$RLMM[[i]], Ncomp = Ncomp)
    }
    
    # Fit RWAMMI model with the current weight
    RWAMMI[[paste0("RWAMMI_", names(weights$RLMM)[i])]] = 
      rwammi.model(data_contaminated, weight = weights$RLMM[[i]], Ncomp = Ncomp, rSVD = rSVD)
  }
  
  # Create a list to store all models
  ALLModels = list(RealSVD = Estimated.real, "AMMI" = AMMI.model, "WAMMI" = WAMMI,
                   "RAMMI" = RAMMI, "RWAMMI" = RWAMMI)
  
  return(ALLModels)
}

# The GEI.Metrics calculates all metrics for all models
GEI.Metrics <- function(allSVDs, Ncomp = 2){
  Metrics = list()
  
  # Extract residual and fitted values from the theoric real data 
  realresidual = allSVDs$RealSVD$residual
  realdf = allSVDs$RealSVD$matrix.fitted + realresidual
  
  # Extract residual and fitted values from the AMMI model
  ammiresidual = allSVDs$AMMI$reduced.SVD
  ammidf = allSVDs$AMMI$matrix.fitted + ammiresidual
  
  Ngen = nrow(realresidual)
  Nenv = ncol(realresidual)
  
  # Calculate singular values and metrics for AMMI
  real.singvals = svd(realresidual)$d
  ammi.singvals = svd(ammiresidual)$d
  aux.mse = (ammi.singvals[1:Ncomp] - real.singvals[1:Ncomp])^2
  
  MPEV = (sum(ammi.singvals[1:Ncomp]^2) / sum(real.singvals^2)) * 100
  MSE.Singvals = aux.mse |> round(digits = 5)
  MtSPE = sum(sort((ammidf - realdf)^2)[1:(Ngen * Nenv * 0.9)]) / (Ngen * Nenv * 0.9)
  I = svd(realresidual)$v[, 1:Ncomp]
  P = svd(ammiresidual)$v[, 1:Ncomp]
  small.eigen = min(eigen(t(I)%*%P%*%t(P)%*%I)$values)
  Maxsub = round(acos(sqrt(small.eigen)) / (pi / 2), 5)
  
  Metrics[["AMMI"]] = list("MPEV" = MPEV,
                           "MSE.Singvals1" = MSE.Singvals[1],
                           "MSE.Singvals2" = MSE.Singvals[2],
                           "MtSPE" = MtSPE,
                           "Maxsub" = Maxsub)
  
  # Calculate metrics for other models
  for(i in 3:length(allSVDs)){
    for(k in 1:length(allSVDs[[i]])){
      aux.matrix = allSVDs[[i]][[k]]$reduced.SVD
      aux.svd = svd(aux.matrix)
      aux.singvals = aux.svd$d
      aux.fitted = allSVDs[[i]][[k]]$matrix.fitted
      aux.est = aux.fitted + aux.matrix
      
      MPEV = (sum(aux.singvals[1:Ncomp]^2) / sum(real.singvals^2)) * 100
      
      MSE.Singvals = (aux.singvals[1:Ncomp] - real.singvals[1:Ncomp])^2 |> round(digits = 5)
      
      MtSPE = sum(sort((aux.est - realdf)^2)[1:(Ngen * Nenv * 0.9)]) / (Ngen * Nenv * 0.9)
      MtSPE = MtSPE |> round(digits = 5)

      I = svd(realresidual)$v[, 1:Ncomp]
      P = aux.svd$v[, 1:Ncomp]
      small.eigen = min(eigen(t(I)%*%P%*%t(P)%*%I)$values)
      Maxsub = round(acos(sqrt(small.eigen)) / (pi / 2), 5)
      
      # Store metrics for the current model
      Metrics[[paste0(names(allSVDs[[i]])[k])]] = list("MPEV" = MPEV,
                                                       "MSE.Singvals1" = MSE.Singvals[1],
                                                       "MSE.Singvals2" = MSE.Singvals[2],
                                                       "MtSPE" = MtSPE,
                                                       "Maxsub" = Maxsub)
    }
  }
  
  return(Metrics)
}

# The BiplotCreation function creates a biplot from the scores and loadings obtained from a Singular Value Decomposition (SVD) analysis.
BiplotCreation <- function(SVD, env.names = NULL, gen.names = NULL,
                           nms = c("IPC1", "IPC2"), lims = FALSE,
                           lims.U, lims.V){
  
  Ncomp <- 2
  
  # Extract scores and loadings from the SVD object
  if(is.matrix(SVD$d)){ # Check if diagonal matrix
    diag <- SVD$d
  } else{
    diag <- diag(SVD$d) # Extract diagonal elements
  }
  scores <- SVD$u[, 1:Ncomp] %*% (diag[1:Ncomp, 1:Ncomp] ^ 0.5) # Calculate scores
  loadings <- SVD$v[, 1:Ncomp] %*% (diag[1:Ncomp, 1:Ncomp] ^ 0.5) # Calculate loadings
  
  # Assign default names to environments and genotypes if not provided
  if(is.null(env.names)){
    env.names <- paste0("Env ", 1:nrow(loadings)) # Generate names as "Env 1", "Env 2", ...
  }
  
  if(is.null(gen.names)){
    gen.names <- 1:nrow(scores) # Generate names as 1, 2, ...
  }
  
  expand <- 1.01
  
  # Calculate limits for the plot axes if lims = FALSE
  if(lims == FALSE){
    lims.U <- c(min(loadings[, 1] * expand, scores[, 1] * expand), # Calculate limits for U axis
                max(loadings[, 1] * expand, scores[, 1] * expand))
    lims.V <- c(min(loadings[, 2] * expand, scores[, 2] * expand), # Calculate limits for V axis
                max(loadings[, 2] * expand, scores[, 2] * expand))
  }
  
  par(pty = "s") # Set aspect ratio to 1
  
  # Plot scores
  plot(scores,
       ylim = lims.V, xlim = lims.U, # Set limits for axes
       type = "n", xlab = nms[1], ylab = nms[2]) # Set plot type, labels
  abline(h = 0) # Add horizontal line at y = 0
  abline(v = 0) # Add vertical line at x = 0
  points(scores, cex = 1.2, pch = 21, bg = "blue", col = "black") # Plot scores as points
  
  par(new = TRUE)
  
  # Plot loadings as arrows
  plot(0, 0,
       ylim = lims.V, xlim = lims.U, # Set limits for axes
       type = "n", axes = FALSE, xlab = "", ylab = "") # Set plot type, remove axes, labels
  arrows(0, 0,
         loadings[, 1], loadings[, 2], # Plot loadings as arrows
         length = 0.1, col = "forestgreen", lwd = 2) # Set arrow properties
  text(loadings[, 1] * 1.01, loadings[, 2] * 1.01,
       labels = env.names, col = "darkgreen", cex = 1.4) # Add environment names
  
}

# The best_envgen function creates a matrix which indicates the best, if quant = 1, genotype for every environment
# if quant = 2, it'll pick up the second best genotype
best_envgen = function(x, quant = 1,
                       gen_names0 = NULL, env_names = NULL,
                       colors = FALSE) {
  
  # Best Environment-Genotype function
  
  # Arguments:
  # - x: Object containing SVD results.
  # - quant: Number of top genotypes to select for each environment (default is 1).
  # - gen_names0: Original genotype names (default is NULL).
  # - env_names: Environment names (default is NULL).
  # - colors: Flag indicating whether to use colors for highlighting (default is FALSE).
  
  Ncomp = 2
  
  if (is.matrix(x$SVD$d) == TRUE) {
    diag_D = x$SVD$d
  } else {
    diag_D = diag(x$SVD$d)
  }
  scores = x$SVD$u[, 1:Ncomp] %*% (diag_D[1:Ncomp, 1:Ncomp] ^ 0.5)
  loadings = x$SVD$v[, 1:Ncomp] %*% (diag_D[1:Ncomp, 1:Ncomp] ^ 0.5)
  
  gen_names = paste0(1:nrow(scores))
  
  if (is.null(env_names)) {
    env_names = paste0("ENV", 1:nrow(loadings))
  }
  
  rownames(scores) = gen_names
  rownames(loadings) = env_names
  
  df_aux = matrix(rep("0", length(env_names)),
                  nrow = length(env_names), ncol = 1)
  rownames(df_aux) = env_names
  colnames(df_aux) = paste0(quant)
  
  data_aux = df_aux
  
  for (env in env_names) {
    
    venvironment = env
    quant
    line_projection = c(-loadings[venvironment, 2],
                        loadings[venvironment, 1],
                        loadings[venvironment, 1],
                        loadings[venvironment, 2])
    
    line_projection = matrix(line_projection, nrow = 2)
    gen_pos = NULL
    for (i in 1:length(gen_names)) {
      vec_coord = matrix(c(0, loadings[venvironment, 1] * scores[i, 1] +
                             loadings[venvironment, 2] * scores[i, 2]), ncol = 1)
      solve_matrix = solve(line_projection, vec_coord)
      gen_pos = rbind(gen_pos, t(solve_matrix))  
    }
    
    values = NULL
    
    for (i in 1:length(gen_names)) {
      values = c(values, sum(loadings[venvironment, ] * gen_pos[i, ]) / 
                   sqrt(sum(loadings[venvironment, ]^2)))
    }
    
    aux = data.frame(values, gen_names)
    best_quant = aux[order(aux[, "values"], decreasing = TRUE), 2][quant]
    data_aux[env, ] = best_quant
    
  }
  
  if (is.null(gen_names0)) {
    gen_names0 = paste0(1:nrow(scores))
  }
  gen_names0 = gen_names0[as.numeric(data_aux)]
  
  if (colors) {
    data_cols = changetocolors(data_aux)
    
    data_aux = paste0(gen_names0, cell_spec("", background = t(data_cols)))
  } else {
    data_aux[, 1] = gen_names0
  }
  
  return(data_aux)
  
}

# The appendList function append values to a list, it helps to append metrics of all the simulations
appendList <- function(x, val) {
  stopifnot(is.list(x), is.list(val)) # Ensure x and val are both lists
  xnames <- names(x) # Get the names of x
  for (v in names(val)) {
    x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
      # If the name v exists in x, and both x[[v]] and val[[v]] are lists, recursively call appendList
      appendList(x[[v]], val[[v]])
    else
      # Otherwise, concatenate the values
      c(x[[v]], val[[v]])
  }
  x
}

# Function to apply a metric to a list of data
# Optionally perform hypothesis testing
Metricsapply = function(list, metric, hypothesis.testing = FALSE, h0 = 100, side = "greater") {
  len = length(list)  # Get the length of the list
  tmp = tmp_aux = test_aux = NULL
  for (i in 1:len) {
    tmp_aux = mean(list[[i]][[metric]])  # Calculate the metric for each element in the list
    tmp = c(tmp, tmp_aux)  # Append the calculated metric to a temporary vector
    
    if (hypothesis.testing) {
      test_aux = c(test_aux, t.test(x = list[[i]][[metric]],
                                    alternative = side,
                                    mu = h0)$p.value)  # Perform hypothesis testing and store the p-value
    }
  }
  
  if (hypothesis.testing) {
    Values = list('Metric' = tmp, 'pvalues' = test_aux)  # Create a list with calculated metrics and p-values
  } else {
    Values = list('Metric' = tmp)  # Create a list with calculated metrics only
  }
  
  return(Values)  # Return the list of calculated metrics and p-values (if applicable)
}