############################################################################### 
###                                                                         ### 
###                               SPR Moose Project                         ### 
###                                                                         ### 
###       Excerpt of code written in Program R to perform statistical       ### 
###       population reconstruction of moose in                             ### 
###       the state of Minnesota                                            ### 
###                                                                         ### 
############################################################################### 

############################################################################### 
############################## Set Up Work Space ############################## 
############################################################################### 

# Clear global environment, import required packages, and set working directory 
{ 
  rm(list=ls()) 
  
  require(BB) 
  require(pso) 
  require(readxl) 
  require(ggplot2) 
  require(numDeriv) 
  require(truncnorm) 
  
  #setwd("~/Dropbox/Published Papers/Active - SPR of Otter in IN/Analysis") 
} 

############################################################################### 
############################# Import Data ############################# 
############################################################################### 

# Import three-tier age-at-harvest matrix 
{ 
  ## Import values as a matrix with each year as a row 
  h <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,28,0,0,25,0,0,26,0,0,25,0,0,38,0,0,29,0,0,37,0,0,24), 
              ncol = 3, byrow = T, 
              dimnames = list(c(2013:2023), c("CALVES","COWS","BULLS"))) 
} 

# Import estimates of yearly catch-effort 
{ 
  f <- c(0,0,0,0.56,0.64,0.75,0.75,0.75,0.87,1.05,0.97)
} 

# Import number available u
{
  n_u <- matrix(c(0,0,1,0,0,3,0,0,5,0,0,22,0,0,10,0,0,10,0,0,11,0,0,17,0,0,20,0,0,0,0,0,0), 
                ncol = 3, byrow = T, 
                dimnames = list(c(2013:2023), c("CALVES","COWS","BULLS")))
}

# Import telemetry hunted
{
  u <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,2,0,0,3,0,0,0,0,0,1,0,0,0,0,0,2,0,0,0,0,0,0), 
              ncol = 3, byrow = T, 
              dimnames = list(c(2013:2023), c("CALVES","COWS","BULLS")))
}

# Import number available v
{
  n_v <- matrix(c(0,95,23,0,90,30,0,67,28,0,45,20,0,21, 7,0,20,10,0,35,10,0,37,17,0,40,18,0, 0, 0,0, 0, 0), 
                ncol = 3, byrow = T, 
                dimnames = list(c(2013:2023), c("CALVES","COWS","BULLS")))
}

# Import non-harvest deaths
{
  v <- matrix(c(0,14,2,0,14,4,0,9,6,0,6,3,0,3,1,0,1,1,0,5,4,0,3,3,0,10,3,0,0,0,0,0,0), 
              ncol = 3, byrow = T, 
              dimnames = list(c(2013:2023), c("CALVES","COWS","BULLS")))
}

# Import aerial estimate
{
  a <- matrix(c(714,1623,2013,439,1513,1498,689,1641,1690,588,1634,1487,428,1156,1446,522,1633,2025,502,1394,1254,NA,NA,NA,885,1967,1849,474,1246,1570,623,1222,1625), 
              ncol = 3, byrow = T, 
              dimnames = list(c(2013:2023), c("CALVES","COWS","BULLS")))
}

# Import standard error aerial estimate
{
  s <- matrix(c(149,339,420,84,288,285,102,242,249,82,228,207,78,211,264,89,277,343,93,258,232,NA,NA,NA,191,425,399,91,239,302,29,254,337), 
              ncol = 3, byrow = T, 
              dimnames = list(c(2013:2023), c("CALVES","COWS","BULLS")))
}

#

## Define number of years and age classes 
{ 
  Y <- nrow(h) 
  A <- ncol(h) 
  yearRange <- c(2013:2023) 
} 

# Define function to calculate the log of a binomial coefficient
{
  binom_coeff_log <- function(n, k) {
    log_n <- lgamma(n + 1)
    log_k <- lgamma(k + 1)
    log_nk <- lgamma(n - k + 1)
    return(log_n - (log_k + log_nk))
  }
}

############################################################################### 
########################## Define Objective Function ########################## 
############################################################################### 

# Define objective function using a multinomial likelihood formulation 
objectiveFunction <- function(par) { 
  
  ## Import initial (diagonal) cohort values 
  { 
    N <- matrix(NA, nrow = Y, ncol = A) 
    N[1, 1:A] <- par[1:A] 
    N[2:Y, 1] <- par[(A + 1):(A + Y - 1)] 
  } 
  
  ## Import vulnerability and survival estimates 
   { 
      C <- matrix(par[Y + (A - 1) + c(1)] / scaleFactor,  
                  nrow = Y, ncol = A, byrow = T) 
      S <- matrix(par[Y + (A - 1) + c(2, 3, 4)] / scaleFactor,  
                  nrow = Y, ncol = A, byrow = T) 
    } 

  
  ## Define probability of harvest 
  { 
    P <- (1 - exp(-C * f)) 
  } 
  
  ## Define expected population sizes 
  { 
    for (i in 1:(Y - 1)) { 
      for (j in 1:(A - 2)) { 
        N[i + 1, j + 1] <- N[i, j] * (1 - P[i, j]) * S[i, j] 
      }  
      
      for (j in (A - 1)) { 
        N[i + 1, j + 1] <- N[i, j]     * (1 - P[i, j])     * S[i, j] + 
          N[i, j + 1] * (1 - P[i, j + 1]) * S[i, j + 1] 
      } 
    } 
  } 
  
  # Return population size estimate (if requested) 
  { 
    if (returnPopulationAbundance) return (N) 
  } 
  
  # Define expected harvest values 
  { 
    E <- N * P 
  } 
  
  ## Calculate contributions of each likelihood component 
  { 
    
    ### Contribution of the age-at-harvest component 
    { 
      logL_AAH <- binom_coeff_log(N,h) + h*log(P) + (N-h)*log(1-P)
      logL_AAH[h==0] = 0
    } 
    
    ### Radio telemetry
    {
     logL_R1 <- binom_coeff_log(n_u,u) + u*log(P) + (n_u-u)*log(1-P)
     logL_R1[u==0] = 0
     logL_R1[n_u==0] = 0
     
     logL_R2 <- binom_coeff_log(n_v,v) + v*log(1-S) + (n_v-v)*log(S)
     logL_R2[v==0] = 0
     logL_R2 [n_v==0] = 0
    }
    
    ### Aerial survey
    {
      logL_AS <- -0.5*log(2*pi) - log(s) - 0.5*((N-a)/s)^2
      logL_AS[is.na(a)] = 0
      logL_AS[is.na(s)] = 0
    }
    
  } 
  
  ## Return value of objective function if valid 
  { 
    logLikelihood <- -c(logL_AAH, logL_R1, logL_R2, logL_AS) 
    
    if (any(is.na(logLikelihood)))     return (9000003) else 
      if (any(logLikelihood == -Inf))  return (9000002) else 
        if (any(logLikelihood == Inf)) return (9000001) else 
          return(sum(logLikelihood)) 
  } 
} 

############################################################################### 
############### Perform NUMERICAL OPTIMIZATION of Observed Data ############### 
############################################################################### 

# Initialize data frames to store results 
{ 
  pointEstimatesBFG <- cbind(data.frame(N11 = rep(NA, 4), N12 = NA, N13 = NA), 
                             data.frame(matrix(NA, nrow = 4, ncol = (Y - 1))), 
                             data.frame(C1 = rep(NA,  4), 
                                        S1 = NA, S2 = NA, S3 = NA, VAL = NA)) 
  populationSizeBFG <- data.frame(matrix(NA, nrow = 4, ncol = Y)) 
} 

# Initialize starting values for parameterization 
{ 
  scaleFactor <- 1e5 
  initialValues <- c(rep(c(mean(h) * 10), Y + (A - 1)), 
                     (c(0.10, 0.70, 0.70, 0.70)) * scaleFactor) 
} 

# Optimize parameter space using a numerical optimization method
{ 
      
  ### Optimize objective function using the BFGS method 
  { 
    returnPopulationAbundance <- FALSE 
    returnChiSquareCorrection <- FALSE 
    optimized <- optim(par = initialValues, fn = objectiveFunction,  
                       method = "L-BFGS-B", lower = 0.001, upper = scaleFactor,
                       control = list(maxit = 1e7, trace = 1)) 
    
    returnPopulationAbundance <- TRUE 
    
    estimatedN <- objectiveFunction(optimized$par) 
    
    pointEstimatesBFG[1, ] <- c(optimized$par, 
                                         optimized$value) 
    populationSizeBFG[1, ] <- rowSums(estimatedN) 
    returnPopulationAbundance <- FALSE 
  } 
 
}

# Extract uncertainty estimates for best-fit model 
{ ## Construct data frame to compile subsequent results 
  { 
    bestFitModel <- data.frame(Year = yearRange, 
                               Estimate = NA, 
                               Lo95CI = NA, Hi95CI = NA, 
                               Measure = rep(c("Abundance",  
                                               "Recruitment"), 
                                             each = length(yearRange))) 
  } 
  
  ## Calculate standard errors for best-fit model parameters 
  { 
    bestParameters <- as.numeric(pointEstimatesBFG[4, 1:(A + Y - 1 + A + A)]) 
    vulnerability <- 3 
    survivability <- 3 
    
    returnPopulationAbundance <- FALSE 
    hessianMatrix <- abs(hessian(x = bestParameters,  
                                 func = objectiveFunction)) 
    
    standardErrors <- sqrt(abs(diag(solve(hessianMatrix)))) 
    
    returnChiSquareCorrection <- TRUE 
    chiSquare <- objectiveFunction(bestParameters) 
    returnChiSquareCorrection <- FALSE 
    
    standardErrors <- standardErrors * sqrt(chiSquare / length(bestParameters)) 
  } 
  
  ## Calculate stochastic abundance estimates using standard errors 
  { 
    abundanceProjection <- matrix(NA, nrow = 1000, ncol = Y) 
    recruiterProjection <- matrix(NA, nrow = 1000, ncol = Y) 
    tempParameters <- bestParameters 
    
    returnPopulationAbundance <- TRUE 
    
    for (i in 1:nrow(abundanceProjection)) { 
      for (j in 1:(A + Y - 1)) { 
        tempParameters[j] <- rtruncnorm(1, a = 0, 
                                        mean = bestParameters[j], 
                                        sd = standardErrors[j]) 
      } 
      for (j in (A + Y - 1 + 1:(2 * A))) { 
        tempParameters[j] <- rtruncnorm(1, a = 0, b = scaleFactor, 
                                        mean = bestParameters[j], 
                                        sd = standardErrors[j]) 
      } 
      
      returnPopulationAbundance <- TRUE 
      abundanceProjection[i,] <- rowSums(objectiveFunction(tempParameters)) 
      recruiterProjection[i,] <- (objectiveFunction(tempParameters))[,1] 
    } 
  } 
  
  ## Compute confidence intervals 
  { 
    bestFitModel$Estimate[bestFitModel$Measure == "Abundance"] <- 
      as.numeric(populationSizePSO[4,]) 
    
    bestFitModel$Lo95CI[bestFitModel$Measure == "Abundance"] <- 
      apply(abundanceProjection, 2, quantile, prob = 0.025) 
    
    bestFitModel$Hi95CI[bestFitModel$Measure == "Abundance"] <- 
      apply(abundanceProjection, 2, quantile, prob = 0.975) 
    
    bestFitModel$Estimate[bestFitModel$Measure == "Recruitment"] <- 
      as.numeric(pointEstimatesPSO[4,c(1, (A - 1) + c(2:Y))]) 
    
    bestFitModel$Lo95CI[bestFitModel$Measure == "Recruitment"] <- 
      apply(recruiterProjection, 2, quantile, prob = 0.025) 
    
    bestFitModel$Hi95CI[bestFitModel$Measure == "Recruitment"] <- 
      apply(recruiterProjection, 2, quantile, prob = 0.975) 
  } 
} 

############################################################################### 
############### Construct VISUALIZATION Best-Fit Reconstruction ############### 
############################################################################### 

# Generate figure of abundance and recruitment estimates 
{ 
  Estimate <- populationSizeBFG[4,]
  ggplot(aes(y = Estimate, x = Year, fill = Measure),  
         data = bestFitModel) +  
    geom_ribbon(aes(ymin = Lo95CI, ymax = Hi95CI,  
                    fill = Measure), alpha = 0.25) + 
    geom_point(aes(colour = Measure)) +  
    geom_path(aes(colour = Measure, linetype = Measure), linewidth = 1) + 
    scale_x_continuous(breaks = seq(min(yearRange), max(yearRange), 1)) + 
    scale_fill_manual(values = c("#F8766D", "#00BA38")) +  scale_color_manual(values = c("#F8766D", "#00BA38")) 
  
} 

############################################################################### 
############### Synthesize ANALYSIS of Reconstruction Estimates ############### 
############################################################################### 

# Compile summary statistics under the likelihood objective function 
{ 
  lambda <- ((bestFitModel$Estimate[bestFitModel$Measure == "Abundance"])[A] /  
               (bestFitModel$Estimate[bestFitModel$Measure == "Abundance"])[1]) ^ 
    (1 / (Y - 1)) 
  
  print(paste("      BEST-FIT MODEL USING PARTICLE SWARM OPTIMIZATION      ")) 
  print(paste("------------------------------------------------------------")) 
  
  print(paste("The model detected", vulnerability,  
              "separate vulnerability coefficient(s):", sep = " ")) 
  print(paste(round(bestParameters[(A + Y - 1 + 1:A)] / scaleFactor, 3), 
              " (SD = ", 
              round(standardErrors[(A + Y - 1 + 1:A)] / scaleFactor, 3), 
              ")", 
              sep = "")) 
  
  print(paste("These corresponded to harvest rates between ", 
              round((1 - exp(-mean(bestParameters[(A + Y - 1 + 1:A)] /  
                                     scaleFactor) * min(f))) * 100, 1), 
              "% and ", 
              round((1 - exp(-mean(bestParameters[(A + Y - 1 + 1:A)] /  
                                     scaleFactor) * max(f))) * 100, 1), 
              "%", sep = "")) 
  
  print(paste("The model detected", survivability,  
              "separate survival coefficient(s):", sep = " ")) 
  print(paste(round(bestParameters[(A + Y - 1 + A + 1:A)] / scaleFactor, 3), 
              " (SD = ", 
              round(standardErrors[(A + Y - 1 + A + 1:A)] / scaleFactor, 3), 
              ")", 
              sep = "")) 
  
  if (lambda > 1) { 
    print(paste("The model detected a positive annual growth rate of",  
                round(lambda, 3), sep = " ")) 
  } else { 
    print(paste("The model detected a negative annual growth rate of",  
                round(lambda, 3), sep = " ")) 
  } 
  
  print(paste("------------------------------------------------------------")) 
}