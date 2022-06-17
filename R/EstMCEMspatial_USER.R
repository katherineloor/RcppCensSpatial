
MCEM.sclm = function(y, x, ci, lcl=NULL, ucl=NULL, coords, init.phi, init.nugget, type="exponential", kappa=NULL,
                     lower=c(0.01,0.01), upper=c(30,30), MaxIter=500, nMin=20, nMax=5000, error=1e-5, show.SE=TRUE){
  n = length(c(y))
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.numeric(x)) stop("x must be a numeric matrix")
  if (!is.matrix(x)) x = as.matrix(x)
  if (!is.matrix(y)) y = as.matrix(y)
  if (det(t(x)%*%x)==0) stop("the columns of x must be linearly independent")

  #No data
  if ((length(x) == 0) | (length(y) == 0) | (length(ci) == 0)) stop("All parameters must be provided")

  #Validating if exists NA's
  if (sum(ci%in%c(0,1)) < length(ci)) stop("The elements of the vector ci must be 0 or 1")
  if (sum(is.na(x)) > 0) stop("There are some NA values in x")
  if (!all(c(is.finite(x)))) stop("x must contain only finite values.")
  if (sum(is.na(ci)) > 0) stop("There are some NA values in ci")
  miss = which(is.na(y))
  if (sum(ci[miss]) != length(miss)) stop ("NA values in y must be specified through arguments ci, lcl, and ucl")
  if (!all(c(is.finite(coords)))) stop("coords must contain only finite values.")

  #Validating dims data set
  if (ncol(as.matrix(y)) > 1) stop("y must have just one column")
  if (ncol(as.matrix(ci)) > 1) stop("ci must have just one column")
  if (nrow(as.matrix(x)) != n) stop("x does not have the same number of lines than y")
  if (length(ci) != n) stop("ci does not have the same length than y")
  if (nrow(coords)!=n | ncol(coords)!=2) stop("Non-conformable dimensions between coords and y.")

  if (sum(ci) > 0){
    if (is.null(lcl) | is.null(ucl)) stop("lcl and ucl must be provided for censored data")
    if (!is.numeric(lcl) | !is.numeric(ucl)) stop("lcl and ucl must be numeric vectors")
    if (length(miss)>0){
      censor = (ci==1 & !is.na(y))
      if (any(is.infinite(lcl[censor])) & any(is.infinite(ucl[censor]))) stop("lcl or ucl must be finite for censored data")
    } else {
      if (any(is.infinite(lcl[ci==1])) & any(is.infinite(ucl[ci==1]))) stop("lcl or ucl must be finite for censored data")
    }
    if (length(lcl) != n) stop("lcl does not have the same length than y")
    if (length(ucl) != n) stop("ucl does not have the same length than y")
    if (ncol(as.matrix(lcl)) > 1) stop("lcl must have just one column")
    if (ncol(as.matrix(ucl)) > 1) stop("ucl must have just one column")
    if (sum(is.na(lcl))>0 | sum(is.na(ucl))>0) stop("There are some NA values in lcl or ucl")
    if (!all(lcl[ci==1]<ucl[ci==1])) stop ("lcl must be smaller than ucl")
  }

  #Validating supports
  if (length(c(init.phi))>1 | !is.numeric(init.phi)) stop("Initial value for phi must be provided")
  if (init.phi <= 0) stop("init.phi must be non-negative")
  if (length(c(init.nugget))>1 | !is.numeric(init.nugget)) stop("Initial value for nugget effect must be provided")
  if (init.nugget <= 0) stop("init.nugget must be non-negative")
  if (is.null(type)) stop("type must be specified")
  if (type!="matern" & type!="gaussian" & type!="pow.exp" & type!="exponential"){
    stop("type should be one of matern, gaussian, pow.exp or exponential")}
  if (type!="exponential" & type!="gaussian"){
    if (is.null(kappa)) stop("kappa must be provided")
    if (length(c(kappa))>1 | !is.numeric(kappa)) stop("kappa must be specified")
    if (type=="pow.exp" & (kappa>2| kappa<=0)) stop("kappa must be a real in (0,2]")
    if (type=="matern" & kappa<=0) stop("kappa must be a real number in (0,Inf)")
    if (type=="matern" & is.infinite(kappa)) stop("kappa must be a real number in (0,Inf)")
  } else { kappa=0 }
  if (length(c(lower))!=2 | length(c(upper))!=2) stop("lower and upper must be vectors of length 2")
  if (any(is.na(lower)) | !is.numeric(lower)) stop("lower is not specified or contains NA")
  if (any(is.na(upper)) | !is.numeric(upper)) stop("upper is not specified or contains NA")
  if (any(lower>=upper)) stop("lower must be smaller than upper")
  if (any(lower<=0)) stop("Values in lower must be non-negative")
  if (!all(is.finite(upper))) stop("upper must contain only finite values")
  if (length(c(MaxIter))>1 | !is.numeric(MaxIter)) stop("MaxIter must be a positive integer value")
  if (MaxIter<=0 | MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
  if (length(c(nMin))>1 | !is.numeric(nMin)) stop("nMin must be a positive integer")
  if (nMin<=1 | nMin%%1!=0) stop("nMin must be a positive integer greater than 1")
  if (length(c(nMax))>1 | !is.numeric(nMax)) stop("nMax must be a positive integer")
  if (nMax<=0 | nMax%%1!=0) stop("nMax must be a positive integer")
  if (nMin > nMax) stop("nMax must be greater than or equal to nMin")
  if (length(c(error))>1 | !is.numeric(error)) stop("error must be specified")
  if (error <= 0) stop("error must be a positive value (suggested to be small)")
  if (!is.logical(show.SE)) stop("show.SE must be logical (TRUE/FALSE).")

  ci = as.matrix(ci)
  lcl = as.matrix(lcl)
  ucl = as.matrix(ucl)
  coords = as.matrix(coords)

  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  out.Sp = MCEM_Spatial(y, x, ci, lcl, ucl, coords, init.phi, init.nugget, type,
                        kappa, lower, upper, MaxIter, nMin, nMax, error, show.SE)
  out.Sp$call = match.call()
  out.Sp$range  = Effective.range(0.05, out.Sp$phi, kappa, type)
  out.Sp$show.SE = show.SE
  out.Sp$ci = ci
  out.Sp$MaxIter = MaxIter

  class(out.Sp) <- "sclm"
  return(out.Sp)
}
