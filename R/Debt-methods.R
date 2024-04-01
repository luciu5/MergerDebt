#'@title Additional methods for LogitDebt Class
#'@description \code{calcShares}, \code{calcSlopes}, \code{calcPrices}, \code{calcMargins}, \code{calcMC},
#'\code{summary} methods for the \code{LogitDebt} class
#' @name LogitDebt-methods
#' @include DebtClasses.R
#' @aliases calcShares, LogitDebt-methods
#' calcSlopes, LogitDebt-methods
#' calcPrices, LogitDebt-methods
#' calcMargins, LogitDebt-methods
#' calcMC, LogitDebt-methods
#' summary, LogitDebt-methods
#' @param object an instance of class \code{LogitDebt}
#' @param preMerger when TRUE, computes equilibrium prices under the pre-merger regime. When FALSE, calculates
#' post-merger equilibrium prices.Default is TRUE.
#' @param level when TRUE, computes margins in dollars. When FALSE, calculates
#' margins as a proportion of prices. Default is FALSE.
#' @param revenue when TRUE, computes revenue shares. When FALSE, calculates
#' quantity shares. Default is FALSE.
#' @param ... harmlessly pass the arguments used in other calcPrices methods to
#' methods for \code{LogitDebt}.
#' @return \code{calcSlopes} return a  \code{LogitDebt}  object containing estimated slopes. \code{CalcQuantities} returns
#' a matrix of equilibrium quantities under `preMerger'.
NULL

setGeneric (
  name= "plotBeta",
  def=function(object,...){standardGeneric("plotBeta")}
)

#'@rdname LogitDebt-methods
#'@export
setMethod(
  f= "plotBeta",
  signature= "LogitDebt",
  definition=function(object){

    shareOutParm <- object@shareOutParm

    #define range
    p = seq(0,1, length=100)

    #create custom plot of Beta distribution
    plot(p, dbeta(p, shareOutParm[1], shareOutParm[2]), ylab='density',
         type ='l', col='purple', main=paste0('Beta(',
                                              round(shareOutParm[1],2),","
                                              ,round(shareOutParm[2],2),")")
         )

  })
#'@rdname LogitDebt-methods
#'@export
setMethod(
  f= "calcShares",
  signature= "LogitDebt",
  definition=function(object,preMerger=TRUE,revenue=FALSE,insideOnly = FALSE,outSideOnly=FALSE){

    if( preMerger) {
      # 1:nFirms
      debt    <- object@debtPre
      prices <- object@pricePre
      owner <- factor(object@ownerPre,levels=unique(object@ownerPre))
      mc <- object@mcPre
    }
    else{
      #firmVec <- unique(object@ownerPost) #2:nFirms
      debt    <- object@debtPost
      prices <- object@pricePost
      owner <- factor(object@ownerPost,levels=unique(object@ownerPost))
      mc <- object@mcPost
    }


    firmVec <-levels(owner)
    nFirms <- nlevels(owner)
    nMarkets <- length(unique(object@marketID))
    slopes <- object@slopes$alpha
    mval <- object@slopes$meanval
    M <- object@mktSize
    shareOutParm <- object@shareOutParm
    outPrice=object@priceOutside
    focal=object@focal
    control.cub <- object@control.cub
    density=object@density
    if(all(!is.na(shareOutParm)) &&
      shareOutParm[1]==1 &&
      shareOutParm[1]==shareOutParm[2]){
        density="uniform"
      }



    #z_crit=mkt$z_crit

    debt <- tapply(debt,owner,sum,na.rm=TRUE)
    debt <- setNames(as.numeric(debt),names(debt))


    if(nMarkets >1){
      prices <- matrix(prices,nrow=nMarkets)
      mc <- matrix(mc,nrow=nMarkets)
      mval <- matrix(mval,nrow=nMarkets)
      owner <- matrix(as.character(owner),nrow=nMarkets)[1,]
    }


    hFocalIntegral <- object@hFocalIntegral
    gFocalIntegral <- object@gFocalIntegral
    tOtherIntegral <- object@tOtherIntegral
    rOtherIntegral <- object@rOtherIntegral


    lowerB <- rep(0,nMarkets-1)
    upperB <- rep(1,nMarkets-1)



    bDens <- function(u,...) {dbeta(u, shape1 = shareOutParm[1],shape2 = shareOutParm[2],...)}
    bProb <- function(u,...) {pbeta(u, shape1 = shareOutParm[1],lower.tail=TRUE,shape2 = shareOutParm[2],...)}

    integrateDerivOtherFun <- function(u,f,pre=preMerger){s0=s0Focal(u,f,pre); u*(1-u)*bProb(s0)}
    integrateInsideOtherFun <- function(u,f,pre=preMerger){s0=s0Focal(u,f,pre); (1-u)*bProb(s0)}

    s0Focal <- function(sOut=rep(NA,nMarkets),i,pre=preMerger,n=nMarkets,l=focal){
      ## allow firm index f to be a vector of length>1.
      ## choose firm 2's product in market 1 to be the reference good
      if(n==1){ res <- 1- debt[i]/profitsCond[,i]; return(res)}

      #res <- (sum(profitsCond[,i,drop=FALSE]) - sum(debt[i]))/profitsCond[l,i[1],drop=TRUE] -
      #  sum(as.vector(sOut) * profitsCond[-l,i,drop=FALSE])/profitsCond[l,i[1],drop=TRUE]

      res <- (sum(profitsCond[,i,drop=FALSE]) - debt[i])/profitsCond[l,i,drop=TRUE] -
        sum(as.vector(sOut) * profitsCond[-l,i,drop=FALSE])/profitsCond[l,i,drop=TRUE]



      return(res)
    }



    if(nMarkets>1) theseStates <-  expand.grid(lapply(1:(nMarkets-1),function(x){1:0}))




    margin <- prices - mc

    sCond <- exp(mval + slopes * prices)
    if(nMarkets>1){sCond <- sCond/rowSums(sCond)}
    else{
         sCond <- sCond/sum(sCond)
         sCond <- t(as.matrix(sCond))
    }

    profitsCond <- M*margin*sCond


    ## aggregate profits to the firm level
    profitsCond <- t(rowsum(t(profitsCond), group = owner, na.rm = TRUE))
    #colnames(profitsCond) <- owner



    if(nMarkets ==1){

    hFun <- sapply(firmVec,function(q){s0_focal=s0Focal(NA,q);hFocalIntegral(s0_focal,s1=shareOutParm[1],s2=shareOutParm[2])})
    gFun <- sapply(firmVec,function(q){s0_focal=s0Focal(NA,q);gFocalIntegral(s0_focal,s1=shareOutParm[1],s2=shareOutParm[2])})

    }
    else if(nMarkets>=2 && density=="uniform"){
      hFun <- sapply(firmVec,function(q,pre=preMerger) {
        s0_focal <- apply(theseStates,1,s0Focal,i=q)

        #betaStates <- hFocalIntegral(s0_focal,s1=shareOutParm[1],s2=shareOutParm[2])

        betaStates <- 6/factorial(nMarkets+1)*s0_focal^(nMarkets+1)
        betaStates <- beta(2,2)*(betaStates - 12/factorial(nMarkets+2)*s0_focal^(nMarkets+2))

        st=ncol(theseStates)
        while(st>=1){
          lastState <- theseStates[,st]
          betaStates <- betaStates[lastState==1]-betaStates[lastState==0]
          theseStates <- unique(theseStates[,-st,drop=FALSE])
          st <- st-1
        }
        #betaStates <- rev((-1)^(1:length(betaStates)))*betaStates
        #betaStates <- sum(betaStates)
        thisProfits <- rowSums(profitsCond[,q,drop=FALSE])

        res <- (-1)^(nMarkets-1)*(prod(thisProfits[-1])/(thisProfits[1]^(nMarkets-1)))*betaStates

        return(res)
      })


      gFun <- sapply(firmVec,function(q,pre=preMerger) {


        s0_focal <- apply(theseStates,1,s0Focal,i=q)
        #betaStates <- gFocalIntegral(s0_focal,s1=shareOutParm[1],s2=shareOutParm[2])
        betaStates <- 2/factorial(nMarkets)*s0_focal^(nMarkets)
        betaStates <-  beta(1,2)*(betaStates - 2/factorial(nMarkets+1)*s0_focal^(nMarkets+1))

        st=ncol(theseStates)
        while(st>=1){
          lastState <- theseStates[,st]
          betaStates <- betaStates[lastState==1]-betaStates[lastState==0]
          theseStates <- unique(theseStates[,-st,drop=FALSE])
          st <- st-1
        }

        #if(!pre && q==2){q <- 2:1}
        thisProfits <- rowSums(profitsCond[,q,drop=FALSE])

        res <- (-1)^(nMarkets+1)*(prod(thisProfits[-1])/(thisProfits[1]^(nMarkets-1)))*betaStates
        return(res)
      })



      if(nMarkets>1){



        otherFun <- function(q,wgt=.5,pre=preMerger) {
          theseStates <-theseStates[theseStates$Var1==1,,drop=FALSE]

          s0_focal <- apply(theseStates,1,s0Focal,i=q)
          if(!pre && q==2){q <- 2:1}
          thisProfits <- rowSums(profitsCond[,q,drop=FALSE])
          betaStates <- s0_focal^(nMarkets - 1)/factorial(nMarkets - 1)

          st=ncol(theseStates)
          while(st>=2){
            lastState <- theseStates[,st]
            betaStates <- betaStates[lastState==1]-betaStates[lastState==0]
            theseStates <- unique(theseStates[,-st,drop=FALSE])
            st <- st-1
          }

          if(nMarkets==2){res <- betaStates}
          else{res <- (-1)^(nMarkets-1)*(prod(thisProfits[-(1:2)])/(thisProfits[1]^(nMarkets-2)))*betaStates}

          res <- res + wgt*(thisProfits[2]/thisProfits[1])
          return(res)
        }

        tFun <- beta(2,2)*sapply(firmVec,otherFun,wgt= 1/2)
        rFun <- beta(1,2)*sapply(firmVec,otherFun,wgt= 2/3)



      }


    }

    else{
      hFun <-
        sapply(firmVec,function(g) {cubature::cubintegrate(function(x,f=g){
          x <- as.matrix(x)
          matrix(apply(x, 2, function(s){
            s0_focal=s0Focal(s,f)
            res <- hFocalIntegral(s0_focal,s1=shareOutParm[1],s2=shareOutParm[2])*prod(bDens(s))
            return(res)
          }),ncol=col(x)
          )}
          ,lower=lowerB,upper=upperB,nVec=control.cub$nVec,
          method= control.cub$method,
          relTol = control.cub$relTol, absTol = control.cub$absTol
          ,maxEval = control.cub$maxEval)$integral
        })

      gFun <-
        sapply(firmVec,function(g) {cubature::cubintegrate(function(x,f=g){
          x <- as.matrix(x)
          matrix(apply(x, 2, function(s){
            s0_focal=s0Focal(s,f)

            gFocalIntegral(s0_focal,s1=shareOutParm[1],s2=shareOutParm[2])*prod(bDens(s))}),ncol=col(x)
          )}
          ,lower=lowerB,upper=upperB,nVec=control.cub$nVec,
          method= control.cub$method,
          relTol = control.cub$relTol, absTol = control.cub$absTol
          ,maxEval = control.cub$maxEval)$integral
        })

      if(nMarkets>1){
        {
          tFun <-
            sapply(firmVec,function(g) {cubature::cubintegrate(function(x,f=g){
              x <- as.matrix(x)
              matrix(apply(x, 2, function(s){

                integrateDerivOtherFun(s,f)*prod(bDens(s))}),ncol=col(x)
              )}
              ,fDim=1,
              lower=lowerB,upper=upperB,nVec=control.cub$nVec,
              method= control.cub$method,
              relTol = control.cub$relTol, absTol = control.cub$absTol
              ,maxEval = control.cub$maxEval)$integral
            })

          rFun <-
            sapply(firmVec,function(g) {cubature::cubintegrate(function(x,f=g){
              x <- as.matrix(x)
              matrix(apply(x, 2, function(s){

                integrateInsideOtherFun(s,f)*prod(bDens(s))}),ncol=col(x)
              )}
              ,fDim=1,
              lower=lowerB,upper=upperB,nVec=control.cub$nVec,
              method= control.cub$method,
              relTol = control.cub$relTol, absTol = control.cub$absTol
              ,maxEval = control.cub$maxEval)$integral
            })

         # tFun <- matrix(rep(tFun,nMarkets-1),byrow = TRUE,nrow = nMarkets-1)
         #  rFun <- matrix(rep(rFun,nMarkets-1),byrow = TRUE,nrow = nMarkets-1)

        }
      }

    }

    if(nMarkets >1){

      #shareOutOther <- exp(log(abs(tFun))-log(abs(rFun)))
      shareOutOther <- tFun/rFun
      shareOutOther[shareOutOther>1] <- 1
      shareOutOther[shareOutOther<0] <- 0
    }



    #replace hFun/gFun with log difference because of stability concerns
    #shareOutFocal <- matrix(exp(log(abs(hFun))-log(abs(gFun))),nrow=1)

    shareOutFocal <- hFun/gFun
    shareOutFocal[shareOutFocal > 1] <- 1
    shareOutFocal[shareOutFocal < 0] <- 0
    #if(nMarkets==2){shareOutOther <- matrix(shareOutOther,nrow=nMarkets - 1)}



    # if(!preMerger){
    #   #shareOutFocal <- c(shareOutFocal[1],shareOutFocal)
    #   shareOutFocal <- shareOutFocal[owner]
    #   if(nMarkets>1) {
    #     #shareOutOther <- cbind(shareOutOther[,1],shareOutOther)
    #     shareOutOther <- shareOutOther[,owner]
    #     }
    # }

    SharesOutInt <- matrix(NA,nrow=nMarkets,ncol=nFirms,dimnames=list(1:nMarkets,firmVec))
    SharesOutInt[focal,]=shareOutFocal
    if(nMarkets>1) SharesOutInt[-focal,]=shareOutOther


    sharesCand <- sCond
    if(outSideOnly){return(SharesOutInt)}
    if(!insideOnly) sharesCand <- sharesCand*(1-SharesOutInt[,owner])

    if(revenue){sharesCand <- prices*sharesCand/sum(prices*sharesCand,object@priceOutside*(1-sum(sharesCand,na.rm=TRUE)),na.rm=TRUE)}

    sharesCand <- as.numeric(sharesCand)
    names(sharesCand) <- object@labels

    return(sharesCand)


  }
)




#'@rdname LogitDebt-methods
#'@export
setMethod(
  f= "calcSlopes",
  signature= "LogitDebt",
  definition=function(object){

    ## Uncover Demand Coefficents

    nMarkets <- length(unique(object@marketID))

    ownerPre     <-  object@ownerPre
    ownerPreMat <- matrix(ownerPre,nrow=nMarkets)
    ownerPreMat <- ownerPreMat[1,]
    nFirms <- length(unique(ownerPre))
    shares       <-  matrix(object@shares,nrow = nMarkets)
    margins      <-  matrix(object@margins,nrow = nMarkets)
    prices       <-  matrix(object@prices,nrow = nMarkets)
    idx          <-   object@normIndex
    focal <- object@focal
    M <- object@mktSize

    debt <- tapply(object@debtPre,ownerPre,sum,na.rm=TRUE)
    debt <- setNames(as.numeric(debt),names(debt))


    shareOutParm <- object@shareOutParm
    parmsStart <- object@parmsStart

    margins <- margins*prices

    nprods <- ncol(object@shares)

    profitsCond <- M*margins*shares


    ## aggregate profits to the firm level
    profitsCond <- t(rowsum(t(profitsCond), group = ownerPreMat, na.rm = TRUE))


    s0Focal <- function(i,sOut){
      if(nMarkets==1){ res <- 1- debt[i]/profitsCond[,i]; return(res)}

      res <- (sum(profitsCond[,i,drop=FALSE]) - debt[i])/profitsCond[focal,i,drop=TRUE] -
        sum(as.vector(sOut) * profitsCond[-focal,i,drop=FALSE])/profitsCond[focal,i,drop=TRUE]

      return(res)
    }


     thisobj <- object

    thisobj@pricePre <- as.numeric(prices)
    thisobj@mcPre <- -1*as.numeric(margins-prices)


    minD <- function(theta){

      alpha <- theta[1:nMarkets]
      if(all(is.na(shareOutParm))){outParm  <- theta[-(1:nMarkets)] }
      else{outParm <- shareOutParm}

      thisobj@slopes <- list(alpha=alpha,
         meanval = as.numeric(log(shares)- log(shares[,idx]) -alpha*(prices - prices[,idx]))
      )

      thisobj@shareOutParm <- outParm

      marginsCand <- MergerDebt::calcMargins(thisobj,preMerger=TRUE,level=TRUE)
      #sharesOutCand <- MergerDebt::calcShares(thisobj,preMerger=TRUE,outSideOnly=TRUE)


      sOutMean <- sapply(ownerPreMat,s0Focal,sOut=outParm[1]/sum(outParm[1]))
      sOutMean[sOutMean<0]=0
      sOutMean[sOutMean>1]=1




      m2 <- m3 <- NA
      m1 <- (as.numeric(margins) - marginsCand)/prices
      if(all(is.na(shareOutParm))){
        m2 <- mean(sOutMean,na.rm=TRUE) - outParm[1]/sum(outParm) #mean of Beta
        m3 <- var(sOutMean,na.rm=TRUE) -  prod(outParm)/(sum(outParm)^2 *(sum(outParm)+1)) #variance of Beta
      }
        #m2 <- mean(sharesOutCand,na.rm=TRUE) - outParm[1]/sum(outParm) #mean of Beta
      #m3 <- var(sharesOutCand,na.rm=TRUE) -  prod(outParm)/(sum(outParm)^2 *(sum(outParm)+1)) #variance of Beta
      #m2 <-as.numeric(shares) - sharesInCand

      measure <- sum((c(m1
                        ,m2
                        ,m3
                        )*100)^2,na.rm=TRUE)

      return(measure)
    }

    ## Constrain optimizer to look  alpha <0,  1e5 >=Beta Parameters >=  1e-5
    lowerB <- c(rep(-1e3,nMarkets),1e-5,1e-5)
    upperB <- c(rep(-1e-10,nMarkets),1e5,1e5)



    if(all(!is.na(shareOutParm))){
      lowerB <- lowerB[1:nMarkets]
      upperB <- upperB[1:nMarkets]
       parmsStart <- parmsStart[1:nMarkets]
      }
    minTheta <- optim(parmsStart,minD,
                      method="L-BFGS-B",
                      lower= lowerB,upper=upperB,
                      control=object@control.slopes)

    minAlpha <- minTheta$par[1:nMarkets]
    minOutParm  <- minTheta$par[-(1:nMarkets)]
    meanval <- log(shares)- log(shares[,idx]) -minAlpha*(prices - prices[,idx])
    meanval <- as.numeric(meanval)
    names(meanval)   <- object@labels

    object@slopes      <- list(alpha=minAlpha,meanval=meanval)
    if(all(is.na(shareOutParm))){
      object@shareOutParm <- minOutParm

      if(any(minOutParm == upperB[-(1:nMarkets)] | minOutParm == lowerB[-(1:nMarkets)])){
        warning("'shareOutParm' estimates are at either upper or lower bound")
      }
      if(any(minOutParm == parmsStart[-(1:nMarkets)])){
        warning("'shareOutParm' estimates are equal to starting values")
      }
      }



    # object@pricePre <- object@prices
    # sharesPre <- calcShares(object,preMerger=TRUE)
    # sharesCondPre <- calcShares(object,preMerger=TRUE,insideOnly=TRUE)
    # sOut <- 1 - sharesPre/sharesCondPre


    return(object)

  }

)



#'@rdname LogitDebt-methods
#'@export
setMethod(
  f= "calcPrices",
  signature= "LogitDebt",
  definition=function(object,preMerger=TRUE,isMax=FALSE,...){


    priceStart <- object@priceStart

    if(preMerger){
      owner <- object@ownerPre
      mc    <- object@mcPre
    }
    else{
      owner <- object@ownerPost
      mc    <- object@mcPost
    }




    #priceEst <- rep(NA,nprods)

    thisObj <- object

    ##Define system of FOC as a function of prices
    FOC <- function(priceCand){

      if(preMerger){ thisObj@pricePre <- priceCand}
      else{          thisObj@pricePost <- priceCand}

      thisMargin <- MergerDebt::calcMargins(thisObj,preMerger=preMerger,level=TRUE)

      thisFOC <- priceCand -mc - thisMargin

      return(thisFOC)
    }

    ## Find price changes that set FOCs equal to 0
    minResult <- BB::BBsolve(priceStart,FOC,quiet=TRUE,control=object@control.equ,...)

    if(minResult$convergence != 0){warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}


    if(isMax){

      hess <- genD(FOC,minResult$par) #compute the numerical approximation of the FOC hessian at optimium
      hess <- hess$D[,1:hess$p]
      hess <- hess * (owner>0)   #0 terms not under the control of a common owner

      state <- ifelse(preMerger,"Pre-merger","Post-merger")

      if(any(eigen(hess)$values>0)){warning("Hessian of first-order conditions is not positive definite. ",state," price vector may not maximize profits. Consider rerunning 'calcPrices' using different starting values")}
    }


    priceEst     <- minResult$par
    names(priceEst) <- object@labels

    return(priceEst)

  }
)
## compute cosntant marginal costs
#'@rdname LogitDebt-methods
#'@export
setMethod(
  f= "calcMC",
  signature= "LogitDebt",
  definition=function(object,preMerger=TRUE){
    if(preMerger){
    object@pricePre <- object@prices

    startMC <- -1*(object@margins - object@prices)
    startMC[is.na(startMC)] <- max(startMC,na.rm=TRUE)

    minD <- function(c){
      object@mcPre <- c
      thisMargin <- MergerDebt::calcMargins(object,preMerger=TRUE,level=TRUE)

      measure <- object@pricePre - c - thisMargin
      measure <- (measure/object@pricePre)^2
      return(sum(measure))
    }

    lowerB <- rep(0,length(startMC))
    upperB <- object@prices

    minMC <- optim(startMC,minD,method = "L-BFGS-B",lower=lowerB,upper = upperB)

    mc <- minMC$par

    }



    if(!preMerger){
      mc <- object@mcPre
      mc <- mc*(1+object@mcDelta)
    }



    names(mc) <- object@labels

    return(mc)
  })

## compute margins
#'@rdname LogitDebt-methods
#'@export
setMethod(
  f= "calcMargins",
  signature= "LogitDebt",
  definition=function(object,preMerger=TRUE,level=FALSE){

    nMarkets <- length(unique(object@marketID))
    alpha <- object@slopes$alpha

    if(preMerger){
      prices <- object@pricePre
      owner <- object@ownerPre
    }

    else{
      prices <- object@pricePost
      owner <- object@ownerPost
      }


    shares <- MergerDebt::calcShares(object,preMerger=preMerger,revenue=FALSE)


    if(nMarkets>1){

      shares <- matrix(shares,nMarkets)
      prices <- matrix(prices,nMarkets)
      owner <- matrix(owner,nMarkets)[1,]


      shares <- t(rowsum(t(shares), group = owner, na.rm = TRUE))
      shares <- shares[,owner]
    }
    else{shares <- tapply(shares, owner,sum,na.rm=TRUE)
    shares <- shares[owner]
    }



    margins <- -1/(alpha*(1-shares))


    if(!level) {
      margins <- margins / prices
    }

    margins <- as.numeric(margins)

    names(margins) <- object@labels

    return(margins)

  }

)

#'@rdname LogitDebt-methods
#'@export
setMethod(
  f= "summary",
  signature= "LogitDebt",
  definition=function(object,revenue=FALSE,shares=TRUE,levels=FALSE,parameters=FALSE,market=FALSE,insideOnly = TRUE,digits=2,...){

    curWidth <-  getOption("width")


    pricePre   <-  object@pricePre
    pricePost  <-  object@pricePost

    priceDelta <- calcPriceDelta(object,levels=levels)

    if(!levels) priceDelta <- priceDelta *100

    if(!shares && !all(is.na(object@prices))){
      outPre  <-  calcQuantities(object,preMerger=TRUE)
      outPost <-  calcQuantities(object,preMerger=FALSE)

      if(revenue){
        outPre <- pricePre*outPre
        outPost <- pricePost*outPost
      }

      sumlabels=paste("quantity",c("Pre","Post"),sep="")
    }

    else{
      if(!shares){warning("'shares' equals FALSE but 'calcQuantities' not defined. Reporting shares instead of quantities")}

      outPre  <-  MergerDebt::calcShares(object,preMerger=TRUE,revenue=revenue,insideOnly=insideOnly) * 100
      outPost <-  MergerDebt::calcShares(object,preMerger=FALSE,revenue=revenue,insideOnly=insideOnly) * 100


      sumlabels=paste("shares",c("Pre","Post"),sep="")
    }

    mcDelta <- object@mcDelta



    if(levels){outDelta <- outPost - outPre}
    else{outDelta <- (outPost/outPre - 1) * 100}


    isParty <- which(object@ownerPost != object@ownerPre)
    partyID <- unique(c(object@ownerPre[isParty],
                        object@ownerPost[isParty]))
    isParty <- object@ownerPre %in% partyID

    isParty <- factor(as.numeric(isParty),levels=0:1,labels=c(" ","*"))

    results <- data.frame(isParty,
                          Market=object@marketID,
                          Product=object@productID,
                          ownerPre=object@ownerPre,
                          ownerPost=object@ownerPost,
                          pricePre=pricePre,pricePost=pricePost,
                          priceDelta=priceDelta,outputPre=outPre,
                          outputPost=outPost,outputDelta=outDelta)

results <- results[order(-1*as.numeric(factor(results$Market)),results$outputPre,decreasing = TRUE),]

    if(sum(abs(mcDelta))>0) results <- cbind(results,mcDelta=mcDelta)

    #rownames(results) <- paste(isParty,object@labels)

    sharesPost <- calcShares(object,FALSE,revenue)

    if(market){

      results$marketID <- object@marketID
      results$ownerPre <- object@ownerPre
      results$ownerPost <- object@ownerPost


      firmIDPre <- firmIDPost <- sort(unique(object@ownerPre))
      firmIDPost[firmIDPost %in% partyID[-1]] <- partyID[1]

      debts <- dplyr::bind_cols(ownerPre= object@ownerPre,
                     ownerPost= object@ownerPost,
                     debtPre=object@debtPre,
                     debtPost=object@debtPost
                     )
      debts <- dplyr::group_by(debts,ownerPre,ownerPost) %>% dplyr::summarise_all(.funs=sum,na.rm=TRUE)

      if(all(debts$debtPre==debts$debtPost)){debts$debtPost <- NULL}

      results <- dplyr::inner_join(results,debts,by=c("ownerPre","ownerPost"))



      profitable <- results[isParty=="*",]%>%
        dplyr::group_by(marketID) %>%
        dplyr::summarise(profitsPost= sum(outputPost * (pricePost-object@mcPost)),
                         profitsPre= sum(outputPre * (pricePre-object@mcPre))
        ) %>%
        dplyr::mutate(leverage=profitsPre,
               isProfitable=(profitsPost-profitsPre)>0
               ) %>%
        dplyr::select(!starts_with("profit"),leverage) %>% ungroup()

      hhisum <- dplyr::ungroup(results)
      hhisum$isParty <- isParty
      hhisum <- dplyr::group_by(hhisum,ownerPre,marketID,isParty)  %>%
        dplyr::summarize(
        #AvgPrePrice=sum(sharesIn*pricePre),
        paascheNum=sum(outputPost * pricePost),
        paascheDen=sum(outputPost * pricePre),
        `Pre-Merger HHI` =sum(outputPre),
        `Pre-Merger Debt HHI` =sum(debtPre)) %>% ungroup()


      hhisum <-   dplyr::group_by(hhisum,marketID)  %>%
        dplyr::summarize(
        `HHI Change`=2*prod(`Pre-Merger HHI`[isParty=="*"]/sum(`Pre-Merger HHI`)*100),
        `HHI Debt Change`= 2*prod((`Pre-Merger Debt HHI`/sum(`Pre-Merger Debt HHI`))[isParty=="*"]*100),
        `Pre-Merger Debt HHI` =sum((`Pre-Merger Debt HHI`/sum(`Pre-Merger Debt HHI`)*100)^2),
        `Pre-Merger HHI` =sum((`Pre-Merger HHI`/sum(`Pre-Merger HHI`)*100)^2),
         paasche=(sum(paascheNum)/sum(paascheDen) - 1)*100,

         ) %>%

        dplyr::mutate( `Post-Merger HHI`= `Pre-Merger HHI`+`HHI Change`,
                `Post-Merger Debt HHI`= `Pre-Merger Debt HHI`+`HHI Debt Change`
        ) %>% dplyr::ungroup() %>%
        dplyr::relocate(paasche,`Pre-Merger HHI`,`Post-Merger HHI`,.after=marketID )



      outAll <- dplyr::inner_join(hhisum,profitable,by = "marketID") %>% dplyr::ungroup() #%>% mutate(leverage=leverage/debt)
      results <- outAll
#
#
#       #thiscmcr <- thiscv <- NA_real_
#
#       #try(thiscmcr <- cmcr(object,levels=levels), silent=TRUE)
#
#       #try(thiscv <- CV(object),silent = TRUE)
#
#       thispsdelta  <- NA_real_
#       try(thispsdelta  <- sum(calcProducerSurplus(object,preMerger=FALSE) - calcProducerSurplus(object,preMerger=TRUE),na.rm=TRUE),silent=TRUE)
#
#       isparty <- isParty == "*"
#
#
#       results <- with(results,
#                       data.frame(
#                         'HHI Change' = as.integer(HHI(outputPre/sum(outputPre),owner=object@ownerPost) - HHI(outputPre/sum(outputPre),owner=object@ownerPre)),
#                         'Industry Price Change (%)' = sum(priceDelta * outputPost/sum(outputPost, na.rm = TRUE),na.rm=TRUE),
#                         'Merging Party Price Change (%)'= sum(priceDelta[isparty] * outputPost[isparty], na.rm=TRUE) / sum(outputPost[isparty], na.rm=TRUE),
#                         # 'Compensating Marginal Cost Reduction (%)' = sum(thiscmcr * outputPost[isparty]) / sum(outputPost[isparty], na.rm=TRUE),
#                         # 'Consumer Harm ($)' = thiscv,
#                         'Producer Benefit ($)' = thispsdelta,
#                         # 'Difference ($)'= thiscv - thispsdelta,
#                         check.names=FALSE
#                       ))
#
#       if(levels){colnames(results) <- gsub("%","$/unit",colnames(results))}
#
#
   }
#
     colnames(results)[colnames(results) %in% c("outputPre","outputPost")] <- sumlabels
#

    #options("width"=ifelse(market,25,100)) # this width ensures that everything gets printed on the same line
    print(results,digits=digits, row.names= FALSE#ifelse(market, FALSE, TRUE)
            )
    #options("width"=curWidth) #restore to current width


    if(!market){

      #results <- cbind(isParty, results)
      #rownames(results) <- object@labels

    cat("\nMerger simulation results under '",class(object),"' demand:\n\n",sep="")



      cat("\n\tNotes: '*' indicates merging parties' products.\n ")
      if(levels){cat("\tDeltas are level changes.\n")}
      else{cat("\tDeltas are percent changes.\n")}
      if(revenue){cat("\tOutput is based on revenues.\n")}
      else{cat("\tOutput is based on units sold.\n")}

    }



    cat("\n\n")


    if(parameters){

      cat("\nDemand Parameter Estimates:\n\n")
      if(is.list(object@slopes)){
        print(lapply(object@slopes,round,digits=digits))
      }
      else{
        print(round(object@slopes,digits))
      }
      cat("\n\n")

      if(.hasSlot(object,"ShareOutParm")){

        cat("\nBeta Shape Parameters:\n\n")
        print(round(object@shareOutParm,digits))
        cat("\n\n")

      }
    }



    return(invisible(results))

  })
