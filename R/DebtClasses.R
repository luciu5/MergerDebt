#'@title  S4 classes to model the effect of debt on mergers
#'@name Debt-classes
#'@import antitrust
#'@import methods
#'@aliases LogitDebt-class
#'@description Extend classes from the \pkg{antitrust} package to accommodate horizontal mergers with debt.
#'
#'@description The \dQuote{LogitDebt} class has the information for a Nash Bertand Pricing Game with Logit demand and firm-level debt. Note that code supports firms operating in multiple markets.
#'@description Let k denote the number of products produced by all firms across all markets. Let f denote the number of firms producing products across all markets. Let m denote the number of markets.
#'
#'@section Objects from the Class:
#'For LogitDebt, objects can be created by using the constructor function \code{\link{logit.debt}}.
#'
#'@slot marketID A length k character vector identifying which market a product is sold in.
# #'@slot marketWgt A length m numeric vector identifying the relative market size of a market, compared to the focal market.
#'@slot productID A length k character vector identifying a product.
#'@slot debtPre A length f vector of pre-merger debt payments.
#'@slot debtPost A length f vector of post-merger debt payments. Default is \sQuote{debtPre}.
#'@slot density A length 1 character vector equal to either "beta' (default) or "uniform". If "beta" calibration routine assumes that firm outside shares are sampled from a \code{\link[base]{Beta}} distribution and uses margin information to calibrate \sQuote{shape1} and \sQuote{shape2}. If "uniform", then \sQuote{shape1}=\sQuote{shape2}=1.
#'@slot focal An integer (default 1) indicating which market is the "focal" market.
#'@slot parmsStart A length $m + 2$ vector whose first $m$ elements are starting values for the price coefficients for each of the m markets. The final 2 elements are starting values for \sQuote{shape1} and \sQuote{shape2} of the \code{\link[base]{Beta}} distribution.
#'@slot shareOutParm TBD
#'@slot hfocalIntegral TBD
#'@slot gfocalIntegral  TBD
#'@slot tOtherIntegral TBD
#'@slot rOtherIntegral TBD
#'@slot control.cub TBD
#'
#'@section Extends:
#'LogitDebt: Class \code{\linkS4class{Logit}}, directly.
#'Class \code{\linkS4class{Bertran}}, by class \code{\linkS4class{Logit}}, distance 2.
#'

#'@author Charles Taragin \email{ctaragin+mergerdebtr@gmail.com}
#'@examples
#'showClass("LogitDebt")           # get a detailed description of the class

#'@rdname Debt-classes
#'@export
setClass(
  Class   = "LogitDebt",
  contains="Logit",
  representation=representation(
    marketID = "character",
    #marketWgt= "numeric",
    productID =  "character",
    debtPre  =   "numeric",
    debtPost =   "numeric",
    density=    "character",
    focal =      "numeric",
    parmsStart = "numeric",
    shareOutParm = "numeric",
    hFocalIntegral = "function",
    gFocalIntegral = "function",
    tOtherIntegral = "function",
    rOtherIntegral = "function",
    control.cub ="list"
  ),
  prototype=prototype(
    normIndex         =  1,
    shareInside       = 1,
    focal             =  1,
    shareOutParm      = c(NA_real_,NA_real_),
    density=      "uniform",
    control.slopes = list(
      factr = 1e7
    ),
    control.cub=list(method="pcubature",
                     nVec=1024L,
                     relTol=1e-5,
                     absTol=2e-5,
                     maxEval=10^6)
  ),

  validity=function(object){

    nProducts =length(object@productID)
    nFirms=length(unique(object@ownerPre))
    nMarkets=length(unique(object@marketID))

    if(length(object@debtPre) != nProducts ||
       length(object@debtPre) != length(object@debtPost)){
      stop("'debtPre' and 'debtPost' must equal the number of products")}

    if(length(object@mktSize) != nMarkets){stop("'mktSize' must be of length ",nMarkets)}
    #if(any(is.na(object@mktSize) | object@mktSize < 0)  ){stop("elements of mktSize must not be NA or negative")}

    #if(length(object@marketWgt) != nMarkets){stop("'marketWgt' must be of length ",nMarkets)}
    #if(any(is.na(object@marketWgt) | object@marketWgt < 0)  ){stop("elements of marketWgt must not be NA or negative")}
    #if(!isTRUE(all.equal(object@marketWgt[object@focal],1,check.attributes = FALSE))){stop("'marketWgt' values must be expressed relative to focal market")}

    nMargins  <- length(object@margins[!is.na(object@margins)])

    if(nMargins< nMarkets+2){stop("At least", nMarkets+2, "elements of 'margins' must not be NA in order to calibrate demand parameters")}

    if(!any(is.na(object@shareOutParm)) && (max(object@shareOutParm)<=0 || length(object@shareOutParm) !=2)){stop("'shareOutParm' must be a length-2 numeric vector with positive values")}
    if(!isTRUE(all.equal(as.numeric(tapply(object@shares, object@marketID,sum,na.rm=TRUE)),rep(1,nMarkets),check.attributes = FALSE,tolerance=1e-3))){
      stop("sum of 'shares' must equal 1")
    }
    if(any(object@margins<0|object@margins>1,na.rm=TRUE)){stop("'margins should be expressed as a proportion between 0 and 1'")}
    if(length(object@parmsStart)!=nMarkets+2){
      stop("'parmsStart' must be a vector of length ",nMarkets+2)
    }
    return(TRUE)
  }
)

