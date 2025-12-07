parHlme = function(
    data, 
    fixedi, # Fixed effects
    randomi, # Random effects
    mixturei = NULL, # Mixture effects
    subject = "ID", 
    ng, # ng = K classes
    ncores = parallel::detectCores() - 1, 
    modelCache = NULL,
    seeds = 1:4 # Fits model 4 times. 
) {
  
  if (is.null(modelCache)) modelCache <- list()
  
  if (ng == 1) { 
    m1 =
      hlme(
        fixed   = fixedi,
        random  = randomi,
        subject = subject,
        ng      = 1,
        data    = data
      )
    
    modelCache[[1]] = m1
    
    return(list(
      models = list(m1),
      `Initial model` = NULL,
      cache = modelCache
    )
    )
  }
  
  
  message("Fitting model with ng = ", ng - 1, " to obtain valid initialization...")
  
  recurRun = parHlme(
    data      = data,
    fixedi    = fixedi,
    randomi   = randomi,
    mixturei  = mixturei,
    subject   = subject,
    ng        = ng - 1,
    ncores    = ncores,
    seeds     = seeds,
    modelCache = modelCache
  )
  
  modelCache <- recurRun$cache
  init_model <- recurRun$cache[[ng - 1]]

  cl = parallel::makeCluster(ncores)
  parallel::clusterEvalQ(cl, library(lcmm))
  
  parallel::clusterExport(
    cl,
    varlist = c("data","fixedi","randomi","mixturei","subject", "ng", "init_model"),
    envir = environment()
  )
  
  models = parallel::parLapply(cl, seeds, function(s) {
    set.seed(s)
    try(
      hlme(
        fixed   = fixedi,
        random  = randomi,
        mixture = mixturei,
        subject = subject,
        ng      = ng,
        data    = data,
        B       = init_model
      ),
      silent = TRUE
    )
  })
  
  stopCluster(cl)
  
  goodModels <- Filter(function(x) {
    inherits(x, "hlme") && !is.null(x$loglik) && is.finite(x$loglik)
  }, models)
  
  if (length(goodModels) == 0)
    stop("All seeds failed for ng = ", ng)
  
  modelBIC = sapply(goodModels, function(m) {m$BIC})
  modelAIC = sapply(goodModels, function(m) {m$AIC})
  modelLoglik = sapply(goodModels, function(m) {m$loglik})
  
  bestBIC = goodModels[[which.min(modelBIC)]]
  bestAIC = goodModels[[which.min(modelAIC)]]
  bestLoglik = goodModels[[which.max(modelLoglik)]]
  
  modelCache[[ng]] <- bestBIC
  
  return(parHlmeResults = list(
    bestModels = list(
      byAIC = list(
        model = bestAIC,
        score = min(modelAIC)
        ),
      byBIC = list(
        model = bestBIC,
        score = min(modelBIC)
      ),
      byLoglik = list(
        model = bestLoglik,
        score = max(modelLoglik)
      )
    ),
    `Initial model` = init_model,
    cache = modelCache
  )
  )
}
