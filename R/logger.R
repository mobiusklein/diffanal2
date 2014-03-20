require(log4r)
.Logger <- setRefClass("Logger", fields=list(name = 'character', .verbose='logical', vars = "list"), methods = list(
  verbose = function(...){
    if(.self$.verbose){
      cat(paste(.self$name, "> ", sep=''), paste(...), '\n')
      if("logger" %in% names(.self$vars)){
        info(.self$vars$logger, paste("verbose:", ...))
      }
    }
  },
  log = function(...){
    if(!"logger" %in% names(.self$vars)) stop(paste("No logfile is configured for", .self$name))
    info(.self$vars$logger, ...)
  },
  log.err = function(...){
    if(!"logger" %in% names(.self$vars)) stop(paste("No logfile is configured for", .self$name))
    error(.self$vars$logger, ...)    
  },
  log.debug = function(...){
    if(!"logger" %in% names(.self$vars)) stop(paste("No logfile is configured for", .self$name))
    debug(.self$vars$logger, ...)
  })
)

Logger <- function(name, verbose = F, logfile = F, vars = list()){
  logger <- NULL;
  if(is.character(logfile)){
    # Create directory path to logfile if it does 
    # not yet exist
    dir.create(dirname(logfile), showWarnings=F)
    logger <- create.logger(logfile)
    # Log everything, from FATAL to DEBUG
    level(logger) <- log4r:::DEBUG
    # `logger` is not a formal class so it must get
    # attached to the vars list rather than be its own 
    # data member.
    vars$logger <- logger
  }
  logger <- .Logger(name = name, .verbose = verbose, vars = vars)
  return(logger)
}