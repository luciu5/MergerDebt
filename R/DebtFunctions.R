#' @title Merger Simulation with Debt
#' @name Debt-functions
#' @aliases debt logit.debt
#' @description Calibrates consumer demand using Logit
#' and then simulates the price effect of a merger between two firms
#' under the assumption that all firms in the market are playing a
#' differentiated products Bertrand pricing game with firm-level debt.
#' @description  Let k denote the number of products produced by all firms across all markets. Let f denote the number of firms producing products across all markets. Let m denote the number of markets.
#' @param shares A length k vector of product (quantity) shares. Values must be
#'   between 0 and 1 .
#' @param prices  A length k vector of  product prices.
#' @param margins  A length k vector of product margins in proportions (bounded between 0 and 1), some of which may
#'   equal NA.
#' @param ownerPre  A length k vector whose values
#'   indicate which firm produced a product pre-merger.
#' @param ownerPost  A length k vector whose values
#'   indicate which firm produced a product post-merger.
#' @param debtPre a length k vector of pre-merger FIRM-level debt payments. Note that sum of product-level debt payments should equal firm level debt.
#' @param debtPost a length k vector of post-merger debt payments. Default equals \sQuote{debtPre}.
#' @param marketID A length k vector whose values
#'   indicate which market a product is sold in.
#' @param productID A length k vector whose values
#'   indicate a product's identity.
#' @param mcDelta A length k vector where each element equals the
#'   proportional change in a downstream firm's product-level marginal costs due to
#'     the merger. Default is 0, which assumes that the merger does not
#'     affect any products' marginal cost.
#' @param mktSize A length m vector equal to market size.
#'   If shares sum to one, this also equals the size of the market.
#@param subset  A length k vector where each element equals TRUE if
#the product indexed by that element should be included in the
#post-merger simulation and FALSE if it should be excluded.Default is a
#matrix where each element is TRUE.
#' @param priceOutside A length 1 vector indicating the price of the
#' outside good. Default is 0.
#' @param priceStart  A length k vector of  starting values used to solve for
#' downstream equilibrium prices. Default is the \sQuote{prices} matrix.
#' @param control.slopes A list of  \code{\link{optim}}
#' control parameters passed to the calibration routine
#' optimizer (typically the \code{calcSlopes} method).
#' @param control.equ A list of  \code{\link[BB]{BBsolve}} control parameters
#' passed to the non-linear equation solver (typically the \code{calcPrices} method).
#' @param labels A k-length vector of labels. Default is "Prod#", where
#' \sQuote{#} is a number between 1 and the length of \sQuote{prices}.
#' @param ... Additional options to feed to the \code{\link[BB]{BBsolve}}
#' optimizer used to solve for equilibrium prices.
#'
#' @details Using product prices, quantity shares and all of the
#' product margins from at least one firm, \code{logit} is able to
#' recover the price coefficient and product mean valuations in a
#' Logit demand model. \code{logit} then uses these
#' calibrated parameters to simulate a merger between two firms.
#'
#'
#' @return \code{logit.debt} returns an instance of class
#' \code{\linkS4class{LogitDebt}}.
#' @author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#' @references
#' @examples
#' ## 3 single-product firms, each selling in 2 markets.
#' ## Firms 1 and 2 merge
#' nprods <- 3
#' nmarkets <- 2
#' marketId <- rep(1:nmarkets,nprods)
#' shares <- c(0.156, 0.418, 0.616, 0.424, 0.229, 0.158)
#' prices <- c(97.303, 119.083, 94.103, 147.905, 93.113, 100.935)
#' margins <-c(17.163, 84.026, 26.451, 103.05, 17.316, 64.105)/prices
#' ownerPre <- ownerPost <- as.character(rep(c(1, 2, 3),each=nmarkets))
#' ownerPost[ownerPost=='2'] <- '1'
#' debt <- c(22.69,NA, 51.002,NA, 2.922,NA)
#' mktSize <- 1
#'
#'

#' simres_debt <- logit.debt(prices,shares,margins,
#'                              ownerPre,ownerPost,debtPre=debt,marketID=marketId,
#'                              mktSize=mktSize)
#'
#'
#' print(simres_debt)
#' summary(simres_debt)
#'
#'@include Debt-methods.R
NULL

#' @rdname Debt-functions
#' @export

logit.debt <- function(prices,shares,margins,
                       ownerPre,ownerPost,
                       debtPre = rep(0, length(prices)),
                       debtPost = debtPre,
                       marketID=as.character(rep(1,length(prices))),
                       productID,
                       #marketWgt=rep(1,unique(marketID)),
                       mktSize=rep(1,unique(marketID)),
                       density=c("beta","uniform"),
                       normIndex=1,
                       focal=1,
                       mcDelta=rep(0, length(prices)),
                       priceOutside=0,
                       priceStart = prices,
                       parmsStart,
                       shareOutParm,
                       isMax=FALSE,
                       control.slopes,
                       control.equ,
                       control.cub,
                       labels,
                       ...
){

  density=match.arg(density)
  if(missing(shareOutParm)) shareOutParm <- rep(NA_real_,2)
  if(density=="uniform"){shareOutParm <- rep(1,2)}

  if(missing(productID)){
    nFirmProds <- tapply(prices,ownerPre,length)
    productID <- unlist(sapply(nFirmProds,function(x){return(1:x)}),use.names = FALSE)
    productID <- paste0(ownerPre,productID)
  }

  mktData <- data.frame(marketID,productID,prices,shares,margins,ownerPre,ownerPost,debtPre,debtPost,mcDelta#,subset
                        )

  if(missing(parmsStart)){
    ## set starting value of alpha to average of single-product equilibrium margin
    parmsStart <- with(mktData,tapply(-1/(margins*prices*(1-shares)),marketID,mean,na.rm=TRUE))
    ## set starting values of shareOutParm to Beta(4,4)
    parmsStart <- c(as.numeric(parmsStart),4,4)
  }
  ## fill in all missing market/owner/product combinations
  ## necessary for reshaping
  # allPerms <- expand.grid(marketID=unique(marketID),
  #                         ownerPre=unique(ownerPre),
  #                         productID=unique(productID)
  # )
  #
  #
  # mktData <- merge(allPerms,mktData)

  mktData <- mktData[with(mktData,order(ownerPre,productID,marketID)),]

  if(missing(labels)){labels <- with(mktData,paste(ownerPre,productID,marketID,sep=":"))}
  nMarkets= length(unique(marketID))

  hFocalIntegral <- function(u,s1=shareOutParm[1],s2=shareOutParm[2]){
    beta(s1+1,s2+1)/beta(s1,s2)*
      pbeta(u,shape1=s1+1,shape2=s2+1,lower.tail = TRUE)
  }

  gFocalIntegral <- function(u,s1=shareOutParm[1],s2=shareOutParm[2]){
    beta(s1,s2+1)/beta(s1,s2)*
      pbeta(u,shape1=s1,shape2=s2+1,lower.tail = TRUE)
  }
  tOtherIntegral <- function(u){NULL}
  rOtherIntegral <- function(u){NULL}



  ## Create LogitDebt  container to store relevant data
  result <- new("LogitDebt",prices=mktData$prices, shares=mktData$shares,
                margins=mktData$margins,
                debtPre=mktData$debtPre,
                debtPost=mktData$debtPost,
                #insideSize=insideSize,
                #marketWgt=marketWgt,
                mktSize=mktSize,
                normIndex=normIndex,
                focal=focal,
                density=density,
                shareOutParm=shareOutParm,
                ownerPre=as.character(mktData$ownerPre),
                ownerPost=as.character(mktData$ownerPost),
                mcDelta=mktData$mcDelta,
                subset=rep(TRUE, length(mktData$prices)),
                marketID=as.character(mktData$marketID),
                productID=as.character(mktData$productID),
                #priceOutside=priceOutside,
                priceStart=priceStart,
                parmsStart=parmsStart,
                hFocalIntegral=hFocalIntegral,
                gFocalIntegral=gFocalIntegral,
                tOtherIntegral=tOtherIntegral,
                rOtherIntegral=rOtherIntegral,
                labels=labels)

  if(!missing(control.slopes)){
    result@control.slopes <- control.slopes
  }
  if(!missing(control.equ)){
    result@control.equ <- control.equ
  }

  if(!missing(control.cub)){
    result@control.cub <- control.cub  }

  ## Convert ownership vectors to ownership matrices
  #result@ownerPre  <- ownerToMatrix(result,TRUE)
  #result@ownerPost <- ownerToMatrix(result,FALSE)
  ## Calculate Demand Slope Coefficients
  result <- MergerDebt::calcSlopes(result)



  ## Calculate marginal cost
  result@mcPre <-  MergerDebt::calcMC(result,TRUE)
  result@mcPost <- MergerDebt::calcMC(result,FALSE)

  ## Solve Non-Linear System for Price Changes
  result@pricePre  <- MergerDebt::calcPrices(result,preMerger=TRUE,isMax=isMax,...)
  result@pricePost <- MergerDebt::calcPrices(result,preMerger=FALSE,isMax=isMax,...)


  return(result)

}
