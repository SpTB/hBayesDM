choiceRT_lba2 = function (data = "choose", niter = 3000, nwarmup = 1000, nchain = 2, 
          ncore = 2, nthin = 1, inits = "random", indPars = "mean", 
          saveDir = NULL, modelRegressor = FALSE, vb = FALSE, inc_postpred = FALSE, 
          adapt_delta = 0.95, stepsize = 1, max_treedepth = 10) 
{
  if (modelRegressor) {
    stop("** Model-based regressors are not available for this model **\n")
  }
  else {
    modelPath <- system.file("stan", "choiceRT_lba2.stan", 
                             package = "hBayesDM")
  }
  startTime <- Sys.time()
  if (data == "example") {
    data <- system.file("extdata", "choiceRT_exampleData.txt", 
                        package = "hBayesDM")
  }
  else if (data == "choose") {
    data <- file.choose()
  }
  if (file.exists(data)) {
    rawdata <- read.table(data, header = T)
  }
  else {
    stop("** The data file does not exist. Please check it again. **\n  e.g., data = '/MyFolder/SubFolder/dataFile.txt', ... **\n")
  }
  subjList <- unique(rawdata[, "subjID"])
  numSubjs <- length(subjList)
  numPars <- 4
  POI <- c("mu_d", "mu_A", "mu_v", "mu_tau", "sigma_d", "sigma_A", 
           "sigma_v", "sigma_tau", "d", "A", "v", "tau", "log_lik")
  if (inc_postpred) {
    POI <- c(POI, "y_pred")
  }
  modelName <- "choiceRT_lba2"
  cat("\nModel name = ", modelName, "\n")
  cat("Data file  = ", data, "\n")
  cat("\nDetails:\n")
  if (vb) {
    cat(" # Using variational inference # \n")
  }
  else {
    cat(" # of chains                   = ", nchain, "\n")
    cat(" # of cores used               = ", ncore, "\n")
    cat(" # of MCMC samples (per chain) = ", niter, "\n")
    cat(" # of burn-in samples          = ", nwarmup, "\n")
  }
  cat(" # of subjects                 = ", numSubjs, "\n")
  Tsubj <- as.vector(rep(0, numSubjs))
  for (i in 1:numSubjs) {
    curSubj <- subjList[i]
    Tsubj[i] <- sum(rawdata$subjID == curSubj)
  }
  maxTrials <- max(Tsubj)
  cat(" # of (max) trials across subjects = ", maxTrials, "\n\n")
  num_choices <- length(unique(rawdata$choice))
  num_cond <- length(unique(rawdata$condition))
  n_tr_cond <- array(NA, dim = c(numSubjs, num_cond))
  for (j in 1:num_cond) {
    n_tr_cond[, j] <- with(rawdata, aggregate(condition == 
                                                j, by = list(subjID = subjID), FUN = sum)[["x"]])
  }
  max_tr <- max(n_tr_cond)
  RT <- array(-1, dim = c(numSubjs, num_cond, 2, max_tr))
  for (subj in 1:numSubjs) {
    for (cond in 1:num_cond) {
      for (choice in 1:num_choices) {
        tmp <- subset(rawdata, rawdata$subjID == subjList[subj] & 
                        rawdata$condition == cond & rawdata$choice == 
                        choice)
        tmp_trials <- n_tr_cond[subj, cond]
        RT[subj, cond, 1, 1:tmp_trials] <- tmp$RT
        RT[subj, cond, 2, 1:tmp_trials] <- tmp$choice
      }
    }
  }
  dataList <- list(N = numSubjs, N_tr_cond = n_tr_cond, N_choices = num_choices, 
                   N_cond = num_cond, RT = RT, Max_tr = max_tr)
  if (inits[1] != "random") {
    if (inits[1] == "fixed") {
      inits_fixed <- c(0.25, 0.75, 2, 0.2)
    }
    else {
      if (length(inits) == numPars) {
        inits_fixed <- inits
      }
      else {
        stop("Check your inital values!")
      }
    }
    genInitList <- function() {
      list(mu_d = inits_fixed[1], mu_A = inits_fixed[2], 
           mu_v = inits_fixed[3], mu_tau = inits_fixed[4], 
           sigma_d = 1, sigma_A = 1, sigma_v = 1, sigma_tau = 1, 
           d = rep(inits_fixed[1], numSubjs), A = rep(inits_fixed[2], 
                                                      numSubjs), v = rep(inits_fixed[3], numSubjs), 
           tau = rep(inits_fixed[4], numSubjs))
    }
  }
  else {
    genInitList <- "random"
  }
  rstan::rstan_options(auto_write = TRUE)
  if (ncore > 1) {
    numCores <- parallel::detectCores()
    if (numCores < ncore) {
      options(mc.cores = numCores)
      warning("Number of cores specified for parallel computing greater than number of locally available cores. Using all locally available cores.")
    }
    else {
      options(mc.cores = ncore)
    }
  }
  else {
    options(mc.cores = 1)
  }
  cat("************************************\n")
  cat("** Building a model. Please wait. **\n")
  cat("************************************\n")
  m = rstan::stan_model(modelPath)
  if (vb) {
    fit = rstan::vb(m, data = dataList, pars = POI, init = genInitList)
  }
  else {
    fit = rstan::sampling(m, data = dataList, pars = POI, 
                          warmup = nwarmup, init = genInitList, iter = niter, 
                          chains = nchain, thin = nthin, control = list(adapt_delta = adapt_delta, 
                                                                        max_treedepth = max_treedepth, stepsize = stepsize))
  }
  parVals <- rstan::extract(fit, permuted = T)
  if (inc_postpred) {
    parVals$y_pred[parVals$y_pred == -1] <- NA
  }
  d <- parVals$d
  A <- parVals$A
  v <- parVals$v
  tau <- parVals$tau
  allIndPars <- array(NA, c(numSubjs, numPars + (num_cond * 
                                                   num_choices) - 1))
  allIndPars <- as.data.frame(allIndPars)
  for (i in 1:numSubjs) {
    if (indPars == "mean") {
      allIndPars[i, ] <- c(mean(d[, i]), mean(A[, i]), 
                           as.vector(apply(v[, i, , ], c(2, 3), mean)), 
                           mean(tau[, i]))
    }
    else if (indPars == "median") {
      allIndPars[i, ] <- c(median(d[, i]), median(A[, i]), 
                           as.vector(apply(v[, i, , ], c(2, 3), median)), 
                           median(tau[, i]))
    }
    else if (indPars == "mode") {
      allIndPars[i, ] <- c(estimate_mode(d[, i]), estimate_mode(A[, 
                                                                  i]), estimate_mode(v[, i]), estimate_mode(tau[, 
                                                                                                                i]))
    }
  }
  allIndPars <- cbind(allIndPars, subjList)
  colnames(allIndPars) <- c("d", "A", apply(expand.grid(paste0("v_cd", 
                                                               1:num_cond), paste0("_ch", 1:num_choices)), 1, paste, 
                                            collapse = ""), "tau", "subjID")
  modelData <- list(modelName, allIndPars, parVals, fit, rawdata)
  names(modelData) <- c("model", "allIndPars", "parVals", "fit", 
                        "rawdata")
  class(modelData) <- "hBayesDM"
  endTime <- Sys.time()
  timeTook <- endTime - startTime
  if (!is.null(saveDir)) {
    currTime <- Sys.time()
    currDate <- Sys.Date()
    currHr <- substr(currTime, 12, 13)
    currMin <- substr(currTime, 15, 16)
    timeStamp <- paste0(currDate, "_", currHr, "_", currMin)
    dataFileName = sub(pattern = "(.*)\\..*$", replacement = "\\1", 
                       basename(data))
    save(modelData, file = file.path(saveDir, paste0(modelName, 
                                                     "_", dataFileName, "_", timeStamp, ".RData")))
  }
  cat("\n************************************\n")
  cat("**** Model fitting is complete! ****\n")
  cat("************************************\n")
  return(modelData)
}