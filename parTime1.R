parHlme = function(
    data, 
    fixedi, 
    randomi, 
    mixturei = NULL,
    subject = "ID",
    ng,
    ncores = parallel::detectCores() - 1,
    seeds = 1:12
) {
  start.time = Sys.time()
  if (ng == 1) {
    return(
      hlme(
        fixed   = fixedi,
        random  = randomi,
        subject = subject,
        ng      = 1,
        data    = data
      )
    )
    message("Runtime: ", round(Sys.time() - start.time, 3), " seconds")
  }
  
  
  message("Fitting model with ng = ", ng - 1, " to obtain valid initialization...")
  
  if (ng - 1 == 1) {
    init_model = hlme(
      data = data,
      fixed = fixedi,
      random = randomi,
      subject = subject,
      ng = 1
    )
  } else {
  init_model = parHlme(
    data    = data,
    fixedi  = fixedi,
    randomi = randomi,
    mixturei = mixturei,
    subject = subject,
    ng      = ng - 1,
    ncores  = ncores
  )
  }
  
  message("The valid initial model is set!")
  message("Runtime: ", round(Sys.time() - start.time, 3), " seconds")
  
  cl = parallel::makeCluster(ncores)
  parallel::clusterEvalQ(cl, library(lcmm))
  
  parallel::clusterExport(
    cl,
    varlist = c("data","fixedi","randomi","mixturei","subject", "ng", "init_model"),
    envir = environment()
  )
  
  message("Fitting model with the valid initial model...")
  
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
  
  good = models[sapply(models, function(x) class(x)[1] == "hlme")]
  
  if (length(good) == 0)
    stop("All seeds failed for ng = ", ng)
  
  message("Runtime: ", round(Sys.time() - start.time, 3), " seconds")
  return(parHlme_results = list(
    models = good,
    `Initial model` = init_model,
    runtime = round(Sys.time() - start.time, 3)
  )
  )
}
