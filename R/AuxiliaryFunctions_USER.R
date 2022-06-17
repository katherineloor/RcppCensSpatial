
dist2Dmatrix = function(coords){

  if (!is.null(coords)){
    coords = as.matrix(coords)
    if (!all(c(is.finite(coords)))) stop ("coords must contain only finite values")
    if (ncol(coords) != 2) stop("coords must contain 2 columns")
    if (nrow(coords) <= 1) stop("coords must be a matrix")
  } else { stop("coords must be specified") }

  distancesM = crossdist(coords)
  out.ST = (distancesM + t(distancesM))/2
  return(out.ST)
}


CovMat = function(phi, tau2, sigma2, dist, type="exponential", kappa=NULL){

  if (length(c(phi))>1 | !is.numeric(phi)) stop("phi must be specified")
  if (phi <= 0) stop("The spatial parameter (phi) must be non-negative")

  if (length(c(tau2))>1 | !is.numeric(tau2)) stop("tau2 must be specified")
  if (tau2 <= 0) stop("The nugget effect (tau2) must be non-negative")

  if (length(c(sigma2))>1 | !is.numeric(sigma2)) stop("sigma2 must be specified")
  if (sigma2 <= 0) stop("The partial sill (sigma2) must be non-negative")

  dist = as.matrix(dist)
  if (!all(c(is.finite(dist)))) stop ("dist must contain only finite values")
  if (ncol(dist) != nrow(dist)) stop("Distance matrix must be specified")
  if (!isSymmetric(dist)) stop("Distance matrix must be symmetric")

  if (is.null(type)) stop("type must be specified")
  if (type!="matern" & type !="gaussian" & type != "pow.exp" & type != "exponential"){
    stop("type should be one of matern, gaussian, pow.exp, exponential")}

  if (type!="exponential" & type!="gaussian"){
    if (is.null(kappa)) stop("kappa mut be provided for pow.exp and matern")
    if (length(c(kappa))>1 | !is.numeric(kappa)) stop("kappa must be specified")
    if (type=="pow.exp" & (kappa > 2| kappa<=0)) stop("kappa must be a real in (0,2]")
    if (type=="matern" & kappa <= 0) stop("kappa must be a real number in (0,Inf)")
    if (type=="matern" & is.infinite(kappa)) stop("kappa must be a real number in (0,Inf)")
  } else { kappa = 0 }

  Covariance = varianceMat(phi, tau2, sigma2, kappa, dist, type)$Sigma
  out.ST = (Covariance + t(Covariance))/2
  return(out.ST)
}


#' @export
predict.sclm = function(object, locPre, xPre, ...){

  if (is.null(object)) stop("object must be specified")

  if (!is.null(locPre) & !is.null(xPre)){
    locPre = as.matrix(locPre)
    xPre = as.matrix(xPre)
    if (!all(c(is.finite(locPre)))) stop("locPre must contain only finite values")
    if (!all(c(is.finite(xPre)))) stop("xPre must contain only finite values")
    if(ncol(xPre)!=length(c(object$beta))) stop("Non-conformable dimensions between xPred and beta")
    if (nrow(locPre)!=nrow(xPre) | ncol(locPre)!=2) stop("Non-conformable dimensions between locPre and xPre")
  } else { stop("locPre and xPre must be specified") }

  ypred = predict.new(object,xPre,locPre)
  out.ST = ypred

  return(out.ST)
}


#' @export
summary.sclm = function(object, ...){
  cat('----------------------------------------------------------------\n')
  cat('           Censored Linear Spatial Regression Model   \n')
  cat('----------------------------------------------------------------\n')
  cat("Call:\n")
  print(object$call)
  cat('\nEstimated parameters:\n')
  # Estimates
  if (all(object$X[,1]==1)) {
    namesx = paste0('\u03b2',0)
    if (ncol(object$X) > 1){ for (i in 2:ncol(object$X)){ namesx = cbind(namesx, paste0('\u03b2',i-1)) } }
  } else {
    namesx = paste0('\u03b2',1)
    if (ncol(object$X) > 1){ for (i in 2:ncol(object$X)){ namesx = cbind(namesx, paste0('\u03b2',i)) } }
  }
  greeks = c(sigma='\u03c3\u00B2', phi='\u03D5', tau='\u03C4\u00B2')
  lab = c(namesx,greeks)

  if (object$show.SE) {
    tab = round(rbind(c(object$theta),c(object$SE)),4)
    colnames(tab) = lab
    rownames(tab) = c("","s.e.")
  } else {
    tab = round(rbind(c(object$theta)),4)
    colnames(tab) = lab
    rownames(tab) = c("")
  }
  print(tab)
  cat('\n')
  cat(paste('The effective range is',round(object$range,4),'spatial units.\n'))
  cat('\nModel selection criteria:\n')
  critFin <- c(object$loglik, object$AIC, object$BIC)
  critFin <- round(t(as.matrix(critFin)),2)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC"))
  print(critFin)
  cat('\nDetails:\n')
  cat('Number of censored/missing values:', sum(object$ci), '\n')
  cat('Convergence reached?:', (object$Iterations < object$MaxIter), '\n')
  cat('Iterations:', object$Iterations,'/',object$MaxIter, '\n')
  cat('Processing time:', round(object$time,4), units(object$time), '\n')
}


#' @export
print.sclm = function(x, ...){
  cat('----------------------------------------------------------------\n')
  cat('           Spatial Censored Linear Regression Model   \n')
  cat('----------------------------------------------------------------\n')
  cat("Call:\n")
  print(x$call)
  cat('\nEstimated parameters:\n')
  # Estimates
  if (all(x$X[,1]==1)) {
    namesx = paste0('\u03b2',0)
    if (ncol(x$X) > 1){ for (i in 2:ncol(x$X)){ namesx = cbind(namesx, paste0('\u03b2',i-1)) } }
  } else {
    namesx = paste0('\u03b2',1)
    if (ncol(x$X) > 1){ for (i in 2:ncol(x$X)){ namesx = cbind(namesx, paste0('\u03b2',i)) } }
  }
  greeks = c(sigma='\u03c3\u00B2', phi='\u03D5', tau='\u03C4\u00B2')
  lab = c(namesx,greeks)

  if (x$show.SE) {
    tab = round(rbind(c(x$theta),c(x$SE)),4)
    colnames(tab) = lab
    rownames(tab) = c("","s.e.")
  } else {
    tab = round(rbind(c(x$theta)),4)
    colnames(tab) = lab
    rownames(tab) = c("")
  }
  print(tab)
  cat('\n')
  cat(paste('The effective range is',round(x$range,4),'spatial units.\n'))
  cat('\nModel selection criteria:\n')
  critFin <- c(x$loglik, x$AIC, x$BIC)
  critFin <- round(t(as.matrix(critFin)),2)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC"))
  print(critFin)
  cat('\nDetails:\n')
  cat('Number of censored/missing values:',sum(x$ci),'\n')
  cat('Convergence reached?:',(x$Iterations < x$MaxIter),'\n')
  cat('Iterations:',x$Iterations,'/',x$MaxIter,'\n')
  cat('Processing time:',round(x$time,4),units(x$time),'\n')
}


#' @export
plot.sclm = function(x, ...){
  plot.convergence(x)
}
