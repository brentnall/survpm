##################################################
## LIBRARIES needed 
##################################################
library("survival")

##################################################
## Global functions
## General formatting functions
##################################################

## return a number with given number of decimal places, rounded to ndigits
fn.format<-function(ind, ndigit=2)
{
    format(round(ind,ndigit), nsmall=ndigit)
}

## format a p-value for reporting depending on number of decimal places and precision wanted
fn.formatp<-function(ind, ndig=3)
{
    if(length(ind)==1){
        if(ind>=0.5){
            roundstop<-1}
        else {
                roundstop<-min(ndig,ndig-sum(sapply(1:ndig, function(idx) 10^idx*ind>=1))+2)
            }
        
        myout<-fn.format(ind, roundstop)
    } else {
        roundstop<-sapply(ind, function(ind2) min(ndig,ndig-sum(sapply(1:ndig, function(idx) 10^idx*ind2>=1))+2))
        roundstop[ind>=0.5]<-1
        myout<-sapply(1:length(ind), function(idx) fn.format(ind[idx], roundstop[idx]))
    }
    
    
    t.inc<-ind < (1/(10^ndig))
    if(sum(t.inc)>0){
        myout[t.inc]<-sapply(myout[t.inc], function(ind) paste("<0",substr(ind,2,ndig+1),"1", sep=""))
    }
    
    myout
}

## format a point estimate and confidence interval of form estimate (95%CI X to Y).
fn.formatCI<-function(ind, ndig=2, ci="95%CI"){

    mynum<-fn.format(ind, ndig)

    paste(mynum[1], " (", ci, " ", mynum[2], " to ", mynum[3], ")", sep="")

}

########################################
## simulate data for spm model
## haz1, haz2 = assumed piecewise constant hazards (Used to simulate data)
## nsim = sample size
fn.spmsim<-function(haz1, haz2, nsim){

    ## simulate piecewise constant survival time
    ## returns maxT + 1 if no event yet
    fn.piecewise<-function(haz, maxT)
    {
        for(idx in 1:maxT){
            thisT<-rexp(1, haz[idx])
            if(thisT<1){
                return((idx-1) + thisT)
            }
        }
        return(maxT)
    }
    

    maxT<-length(haz1)
    
    ##Time 1
    myT1<-sapply(1:nsim, function(ind) fn.piecewise(haz1, maxT))

    ##Time 2
    myT2<-sapply(1:nsim, function(ind) fn.piecewise(haz2, maxT))

    myT<-pmin(myT1, myT2)

    myD <- apply(cbind(maxT, myT1, myT2), 1, which.min)-1

    myT<-pmin(myT, maxT-0.00001)

    myhazmod<-t(matrix(rep(c(cumsum(haz1), cumsum(haz2)), nsim), nrow=length(haz1) + length(haz2)) )

    mysumdta<-data.frame(myid=seq(1, nsim), t = myT, mycause = myD, mymod=myhazmod)

    colnames(mysumdta)<-c("id", "t", "d", paste("H1-", 1:20, sep=""), paste("H2-", 1:20, sep=""))

    mysumdta
}



########################################
## SurvPM class definitions
########################################

####
## * crData format needed:
##     - id, time, cause, cumulative hazard (number discrete time units) cause 1, 2, ..., then covariates used in formula (if supplied)
##Time in same units as hazard. e.g. 1 y hazard, time in years
## inX - covariates, including grp for risk subgroups

## Check validity of class
checkSurvPM <- function(object) {

    errors <- character()
    
    if (length(errors) == 0) TRUE else errors
}



setClass(Class="SurvPM",
         representation=representation(
             crData="data.frame",
             crhaz="matrix",
             crH="matrix",
             m="integer",
             n="integer",
             maxT = "integer",
             pmH = "matrix",
             resOEtot="data.frame",
             resOEtime="list",
             crHt="matrix",
             inX="data.frame",
             nX="integer",
             adjform="character",
             ctime="logical"
             ),
         validity = checkSurvPM
         )



## initialise class object
setMethod(
    f="initialize",
    signature="SurvPM",
    definition=function(.Object, crData, m, maxT, ctime, adjform){

        ##check passed the correct arguments

        ##1. crData shoud have id, survival time, cause, then cumulative hzarad for each competing risk cause
        if(ncol(crData) < 3 + maxT*m){

            myerror <- paste("Input error: crData should have at least", 3+maxT*m, "columns, crData supplied only has", ncol(crData), "columns. Need to follow format with columns for id, survival time, cause, then cumulative hazards for each cause (m=)", m, "supplied),  up to maxT (", maxT, "supplied)." )

            stop(myerror)

        }

        if(sum( colnames(crData)[1:3] == c("id", "t", "d") )< 3){

            mywarning<- paste("Need data frame for crData with the first three columns identifier (use colname id), time (t) and event type (d). You passed", colnames(crData)[1], ",", colnames(crData)[2], ",", colnames(crData)[3],".")

            stop(mywarning)

        }

       
        ## id, survival time and cause
        .Object@crData <- crData[,1:3] 

        ## expected cumulative hazard for each competing risk model
        .Object@crH <- as.matrix(crData[,4:(3+maxT*m)])
       
        ## number of competing risks (e.g. m=1 is survival analysis)
        .Object@m<-as.integer(m)

        ##### error check
        errors <- character()

        errorcount <- 0

        for(idx in 1 : m){
            
            thischeck <- apply( .Object@crH[,seq( (idx - 1)*maxT + 1, (idx - 1)*maxT + maxT)], 1, diff)<0
            
            if(sum(thischeck)>0){
                
            errorcount <- errorcount + 1
                
                errors[errorcount] <- paste("Cumulative hazard should INCREASE! Please check cause m=", idx, "\n")
            }
        }

        if (length(errors) > 0) stop(errors)

######

        
        ## number of individuals
        .Object@n = as.integer(nrow(crData))

        ## maximum number of discrete time periods over which expected cumulative hazard given
        .Object@maxT<-as.integer(maxT)

        ## whether analysis through time sought
        .Object@ctime <- as.logical(ctime)
        
        ## adjustment formula
        .Object@adjform <- as.character(adjform)
        
        ## covariates for group check etc, if supplied
        if(ncol(crData) > 3+maxT*m) {

            .Object@inX <- crData[ , (3+maxT*m + 1) : ncol(crData)]

            .Object@nX <- as.integer( ncol(.Object@inX) )

        }
        else {
            
            .Object@nX <- as.integer(0)

        }

        ## Method to calculated overall expected number events, each cause based on cumulative hazard approach
        ##print("[debug] calcH")
        allhaz<-calcH(.Object)

        ## Overal expected number of events assigned
        .Object@pmH<-allhaz[[1]]

        ## Overall expected number of events per time period assigned
        .Object@crhaz<-allhaz[[2]]

        ## Calculate overall calibration table based on exact Poisson etc
        ##print("[debug] calcOE")
        .Object@resOEtot<-calcOE(.Object)

        ## Calculate time dependent expected cumulative hazard for each cause (i.e. through followup time, at each event or when model hazard may change (discrete time interval).
        ##TODO: this takes a long time when n large (even 1000)
        print("[debug] tdepcal")
        .Object@crHt <-tdepcal(.Object)
        
        if(ctime==TRUE){
            print("[debug] calcOEt")
            .Object@resOEtime <- calcOEt(.Object)
        }
        ## Initialise
        return(.Object)
        
        }
)

## Friendly interface to S4 object
## - df: data frame with specific strucutre
## - m: number of competing risk causes
## - maxT: maximum number of discrete time periods in model expected risks
## - ctime: whether to do the performance analysis through follow-up time or not
## - adjform: additional adjustment terms for Poisson regression (when ctime TRUE)
spm<- survpm <- function(df, m, maxT, ctime=FALSE, adjform = ""){
    thisspm<- new(
        Class="SurvPM",
        crData=df,
        m=as.integer(m),
        maxT=as.integer(maxT),
        ctime=as.logical(ctime),
        adjform=as.character(adjform)
    )

    return(thisspm)
    }

## Class method define. Not intended for users (internal)
## Will define a function that can calculate cumulative hazards for cause-specific model projections
setGeneric(
    name="calcH",
    def=function(object){standardGeneric("calcH")}
    )

## Class method implement
## Calculate cumulative Hazards overall and by time interval
setMethod(
    f="calcH",
    signature="SurvPM",
    definition=function(object){
        
        lastunit <- floor(object@crData$t)
        
        extra <- object@crData$t - lastunit
        
        ##initialise
        nid<-nrow(object@crData)
        
        mycumH<-matrix(0, ncol=object@m, nrow= nid)

        myhaz<-matrix(0, ncol=object@m*object@maxT, nrow=nid) 
        
        for(idx in 1:(object@m)){ ## competing risks
            
            for(idy in 1:nid){ ## individuals
            
                if(lastunit[idy]>0){
                    
                    mycumH[idy, idx] <- object@crH[idy, ((1+(idx-1)*object@maxT) + lastunit[idy]-1)]
                    
                    mycumH[idy, idx] <- mycumH[idy, idx] + extra[idy] * (object@crH[idy, (1+(idx-1)*object@maxT) + lastunit[idy]] - mycumH[idy, idx])
                    
                }
                
            else { ## in first time period
                    
                    mycumH[idy, idx] <- extra[idy] * (object@crH[idy, (1+(idx-1)*object@maxT) + lastunit[idy]])
                    
                }

                ##hazard
                myhaz[idy,  (1 + ( (idx-1) * object@maxT) ) : (1 + ( (idx-1) * object@maxT) + object@maxT-1 )] <- diff( c(0, unlist( object@crH[idy, (1+(idx-1)*object@maxT) : (1+(idx-1)*object@maxT + object@maxT-1)] ) ))
                
           
            }
            
        }

        return(list(mycumH, myhaz))
 }       
)

## Class method define. Not intended for users (internal)
## Overall calibration (O/E) based on cumulative hazards for cause-specific model projections
setGeneric(
    name="calcOE",
    def=function(object){standardGeneric("calcOE")}
)
## Class method implement
## Calculate calibration statistics, O/E overall
setMethod(
    f="calcOE",
    signature="SurvPM",
    definition=function(object){
        myo <- hist(object@crData$d, breaks= seq(0, object@m+1) -0.5, plot=FALSE)$count
        
        mye<-colSums(object@pmH)
        
        pois.exact<-lapply(1:object@m, function(idx) poisson.test(x=myo[idx+1], r=mye[idx]))
        
        myci.fmt<-sapply(pois.exact, function(ind) fn.formatCI(c(ind$estimate, ind$conf.int[1], ind$conf.int[2])/ind$null.value))
        
        myp.fmt<-sapply(pois.exact, function(ind) fn.formatp(ind$p.value,4))
        
        myres <- data.frame(O=myo[2:(object@m+1)], E=mye, OvE=myci.fmt, P=myp.fmt)
        
        colnames(myres)<-c("Observed (O)", "Expected (E)", "O/E", "P")

        return(myres)

    }
)

## Class method define. Not intended for users (internal)
## Time dependent calibration (O/E) based on cumulative hazards for cause-specific model projections
setGeneric(
    name="calcOEt",
    def=function(object){standardGeneric("calcOEt")}
)
## Class method implement
## Calculate time dependent calibration statistics, O/E 
setMethod(
    f="calcOEt",
    signature="SurvPM",
    definition=function(object){

        ## unique times when  hazard changes in model
        myz<-0:object@maxT 
        
        nz<-object@maxT

        lastrisk<-floor(object@crData$t)

        ntot<-sum(lastrisk+1)
        
        myzrem<-object@crData$t - lastrisk

        mydtafit <- data.frame( matrix(NA, ncol = 3 + object@m + object@nX, nrow=ntot ) ) #id, time, hazard cr1, hazard cr2, .. (etc if m>1), outcome, fixed X

        colnames(mydtafit) <- c("id", "t", paste("h", 1 : object@m, sep=""), "J", colnames( object@inX) )

        ## id (as integer)
        mydtafit[, 1] <- rep(1:object@n, lastrisk + 1) 

        ## time period (integer)
        mydtafit[, 2] <- unlist(sapply(lastrisk + 1, function(ind) 1 : ind)) 

        ## covariates
        mydtafit[, (3 + object@m) : (3 + object@m + object@nX) ] <- object@inX[ rep( row.names( object@inX ), lastrisk+1),]
        
        ## add hazards
        for(idx in 1:object@m){ #loops thru competing risks

            mycounter<-0 ## initialise - mycounter is current position on time-dept survival model data frame
            
            for(idy in 1:object@n){ ##loop through individuals

                if( lastrisk[idy]>0 ){ #if one complete time period at start
                    for(idz in 1 : (lastrisk[idy]) ){ #loop through time periods at risk 
        
                        mycounter <- mycounter + 1

                        mydtafit[mycounter, 2 + idx] <- object@crhaz[idy, idz + (idx - 1) * object@maxT]##hazard
                    }
                }

               

                ##last period which might be not whole
                mycounter <- mycounter + 1

                mydtafit[mycounter, 2 + idx] <- myzrem[idy] * object@crhaz[idy, lastrisk[idy] + 1 + (idx - 1) * object@maxT]##haazard while at risk

            }
            
        }
        
        ##event for each person position in the data
        myposlast<-cumsum(1+lastrisk)

        ## competing risk cause, init
        mydtafit$J <- 0

        ## competing risk cause last time period
        mydtafit$J[myposlast] <- object@crData$d

        ## data for fitting calibration slope (only)
        mydtafit2 <- data.frame(J=object@crData$d, t=object@crData$t, H= object@pmH )
        
        ##formulas for calibration check models
        myformula<-vector("list", 4)

        ##1. calibration in the large
        myformula[[1]]<- as.formula("(mydtafit$J==idx) ~ offset(log (mydtafit[,2+idx]) )")

        ##2. calibration slope
##        myformula[[2]] <- as.formula("(mydtafit$J==idx) ~ log (mydtafit[,2+idx])")
        myformula[[2]] <- as.formula("(mydtafit2$J==idx) ~ log (mydtafit2[,2+idx])")

        ##3. overall calibration plus time
        myformula[[3]]<- as.formula( "(mydtafit$J==idx) ~  mydtafit$t + offset(log(mydtafit[,2+idx]))" )

        ##4. subgroups / user adjustments plus time, if user supplied formula
        if(object@adjform!=""){
            
            myformula[[4]] <- as.formula( paste( c("(mydtafit$J==idx) ~ -1 + t +", object@adjform, "+ offset(log(mydtafit[,2+idx]))") ))
            
        }
        

        #return coefficients and confidence intervals only
        fn.formatglm<-function(inglm){

            ## if problems with convergence, return NA (i.e. v small or large calib coefs)
            if( (sum( coef(inglm) < -6) > 0 )  | sum( coef(inglm) > 5 ) >0 ){

                myout<-matrix(NA, ncol=3, nrow=length( coef(inglm) ) )

            } else {

                myout<-cbind(coef(inglm), confint.default(inglm) )

            }
            return(myout)
            
        }
            
        ##initilise
        myglm <- vector("list", object@m)
        
        for(idx in 1 : object@m){

            ## initialise
            myglm[[idx]]<- vector("list", 5)
            
            for(idy in c(1,3)){

                thisglm<-  glm( myformula[[idy]] , mydtafit, family=poisson)

                myglm[[idx]][[idy]] <- fn.formatglm( thisglm )

            }

            idy<-2
            
            thisglm<-  glm( myformula[[idy]] , mydtafit2, family=poisson)

            myglm[[idx]][[idy]] <- fn.formatglm( thisglm )
                        
            ##IF adjmodel provided
            if(object@adjform != ""){

                myglm[[idx]][[4]] <- fn.formatglm( glm( myformula[[5]], mydtafit, family=poisson) )

            }

        }
            
        return(myglm)
        
    }
)


## Class method definition. Not intended for users (internal)
## Defines function that can calc cumulative hazards for cause-specific model projections over follow-up time (each event or when hazard function changes)
## needed??????????
setGeneric(
    name="calcHt",
    def=function(object){standardGeneric("calcHt")}
    )
## details
setMethod(
    f="calcHt",
    signature="SurvPM",
    definition=function(object){
        
    ## unique times when event or hazard changes in model
    myz<-sort(unique(c(object@crData$t, 1:object@maxT))) ##add when hazard changes each time unit

    myzidx<-    floor(myz)

    myzrem<-    myz-floor(myz)

    ##add t=0
    myz<-c(0,myz); myzidx<-c(0,myzidx); myzrem<-c(0,myzrem)

    nz<-length(myz)

    lastrisk<-sapply(object@crData$t, function(idx) which(myz==idx))

    ## hazard from time point idz to idz+1
    fn.intervalH<-function(myH, myh, idz){
        
        myHstart<- myH[ myzidx[idz]+1 ] + myzrem[idz] * myh[myzidx[idz]+1]

        myHend<- myH[ myzidx[idz+1]+1 ] + myzrem[idz+1] * myh[myzidx[idz+1]+1]

        myHinterval<-myHend-myHstart

        myHinterval
    
    }

    myzH<-matrix(0, nrow=nz, ncol=object@m)

    for(idx in 1:object@m){

        for(idy in 1:object@n){

            for(idz in 1:(lastrisk[idy]-1)){

                myzH[idz ,idx] <- myzH[idz ,idx] + fn.intervalH(
                    object@crH[idy, (1 + (idx-1)*object@maxT) : ((idx-1)*object@maxT + object@maxT)],
                    object@crhaz[idy, (1 + (idx-1)*object@maxT) : ((idx-1)*object@maxT + object@maxT)],
                    idz
                    )
            }
        }   
    }

    ## cumulative at risk hazard (eqn 22)
    myzHcuml <- apply(myzH,2, function(ind) cumsum(ind/atrisk))
    
    return(cbind(myz, myzHcuml))

      ##number at risk
        atrisk <- rep(object@n, nz)

        ## need to consider order of last risks
        ##ties
        mytimes <- sort(unique(lastrisk))

        ntimes <- length(mytimes)
        
        changetimes<-cbind(mytimes, hist(lastrisk, c(0,mytimes), plot=FALSE)$count)
        
        for(idz in 1:( ntimes -1 )){

            atrisk[ (mytimes[idz]+1) : (mytimes[idz+1]) ] <- atrisk[ mytimes[idz] ] - changetimes[idz, 2]
        }

        ## last time point
        atrisk[( mytimes[ntimes] + 1):nz]  <- atrisk[ mytimes[ntimes] ] - changetimes[ntimes, 2] 
        
        atrisk[atrisk==0]<-1 # avoid NaN
        
        ## cumulative at risk hazard (eqn 22)
        myzHcuml <- apply(myzH,2, function(ind) cumsum(ind/atrisk))
        
        return(cbind(myz, myzHcuml))
  


        
    }
)

## Class method definition. Not intended for users (internal)
## Will define a function to calculations that enable a time-dependent calibration analysis
setGeneric(
    name="tdepcal",
    def=function(object){standardGeneric("tdepcal")}
    )
## Class method implement.
setMethod(
    f="tdepcal",
    signature="SurvPM",
    definition=function(object){

        ## unique times when event or hazard changes in model
        myz<-sort(unique(c(object@crData$t, 1:object@maxT))) ##add when hazard changes each time unit
        
        myzidx<-    floor(myz)

        myzrem<-    myz-floor(myz)
        
        ##add t=0
        myz<-c(0,myz); myzidx<-c(0,myzidx); myzrem<-c(0,myzrem)
        
        nz<-length(myz)
        
        lastrisk<-sapply(object@crData$t, function(idx) which(myz==idx))
        
        ## hazard from time point idz to idz+1
        fn.intervalH<-function(myH, myh, idz){
                    
            myHstart<- myH[ myzidx[idz]+1 ] + myzrem[idz] * myh[myzidx[idz]+1]
            
            myHend<- myH[ myzidx[idz+1]+1 ] + myzrem[idz+1] * myh[myzidx[idz+1]+1]
            
            myHinterval<-myHend-myHstart
            
            myHinterval
            
        }

        print("[debug] tdepcal 2")
        
        

        ## cumulative hazard should start at zero
        thisH<-list()

        for(idx in 1:object@m){

            thisH[[idx]] <-  cbind( rep(0, object@n), object@crH[,(1 + (idx-1)*object@maxT) : ( (idx-1)*object@maxT + object@maxT)] )  

        }        

        ##initialise
        myzH<-matrix(0, nrow=nz, ncol=object@m)

        ## this is rate limiting step - want to speed up
        ## todo: some people may have exactly same predicted hazard - could use this to speed up (do unique thisH and count number at risk each period)
        for(idx in 1:object@m){

            for(idy in 1:object@n){
                    
                for(idz in 1:(lastrisk[idy]-1)){
                    
                    myzH[idz+1 ,idx] <- myzH[idz+1 ,idx] + fn.intervalH(
                      thisH[[idx]][idy, ],
                       object@crhaz[idy, (1 + (idx-1)*object@maxT) : ((idx-1)*object@maxT + object@maxT)],
                       idz
                    )
            }
            }   
        }

        
        ##number at risk
        atrisk <- rep(object@n, nz)

        ## need to consider order of last risks
        ##ties
        mytimes <- sort(unique(lastrisk))

        ntimes <- length(mytimes)
        
        changetimes<-cbind(mytimes, hist(lastrisk, c(0,mytimes), plot=FALSE)$count)
        
        print("[debug] tdepcal 3")
        
        for(idz in 1:( ntimes -1 )){
##            print(idz)
            atrisk[ (mytimes[idz]+1) : (mytimes[idz+1]) ] <- atrisk[ mytimes[idz] ] - changetimes[idz, 2]
        }

        ## last time point
        atrisk[( mytimes[ntimes] + 1):nz]  <- atrisk[ mytimes[ntimes] ] - changetimes[ntimes, 2] 
        
        atrisk[atrisk==0]<-1 # avoid NaN
        
        ## cumulative at risk hazard (eqn 22)
        myzHcuml <- apply(myzH,2, function(ind) cumsum(ind/atrisk))
        
        return(cbind(myz, myzHcuml))
    
        }
)



## Summary method implementation (for users).
setMethod(
    f="summary",
    signature="SurvPM",
    definition=function(object){

        ##function to format matrix with calib coefs
        fn.formatCIm<-function(ind, ndig=2, notxt=TRUE)
            {
                ind.fmt <- fn.format(exp(ind), ndig)

                if(nrow(ind.fmt)==1){

                    thisout<- paste(ind.fmt[1], " (", ind.fmt[2], " to ", ind.fmt[3], ")", sep="") 
                    
                } else {
                
                    thisout <- paste(ind.fmt[,1], " (", ind.fmt[,2], " to ", ind.fmt[,3], ")", sep="") 
                }

                return(thisout)
                
            }

        
        
        print("**Calibration, overall, exact Poisson test**")

        print(object@resOEtot)

        ## analysis through time done
        if(object@ctime){

            for(idx in 1 : object@m){
                print(paste(c("---->> Cause", idx) ))
                
                print("**Calibration in the large, O/E (95%CI), Poisson regression model**")
                thisout<-fn.formatCIm( object@resOEtime[[idx]][[1]] )
                print(thisout)

                ##Note this is not sensible if constant hazard..
                print("**Calibration slope , O/E (95%CI)**")
                thisout<-fn.formatCIm( object@resOEtime[[idx]][[2]] )
                print(thisout[2])

                print("**Calibration through time, trend test O/E (95%CI)**")
                thisout<-fn.formatCIm( object@resOEtime[[idx]][[3]] )
                print(thisout[2])

                if(object@adjform != ""){
                    print("**Calibration, bespoke, Poisson regression model**")
                    thisout<-fn.formatCIm( object@resOEtime[[idx]][[4]] )
                    thisrownames<-rownames(object@resOEtime[[idx]][[4]])
                    thisoutm<-data.frame(c("Overall", "Follow-up", thisrownames[3:length(thisrownames)]), thisout)
                    colnames(thisoutm)<-c("Parameter", "Coefficient (95CI)")
                    print(thisoutm)
                }
            }

        }

    }

)


## Plot method implementation (for users)
setMethod(
    f="plot",
    signature="SurvPM",
    definition=function(x,y=1,m=1,...){

        if(missing(y)) y<-1

        if(y!=1 & y!=2) {
            
            print("Invalid option for y (type of chart). Must be 1 or 2.")
            
            return("Error - not plotted")
        }
        
        ##NA fit
        plotdta<-data.frame(T=x@crData$t, CC= as.integer(x@crData$d==m))
        
        myNA2<-survfit(Surv(T,CC)~1, plotdta)
        
        myrt<-myNA2$n.event / myNA2$n.risk
        
        mysig<-sqrt(cumsum(myNA2$n.event / (myNA2$n.risk^2)))
        
        myH2<-cumsum(myrt)

        if(y==1){

            plot.ylim<-c(0, max(c( max((myH2+1.96*mysig)), max(x@crHt[, 1+m]))))
            
            plot(x@crHt[,1], x@crHt[,1+m], xlim=c(0,x@maxT), ylim=plot.ylim,  type="l", xlab="Time (y)", ylab="Cumulative Hazard", col=2, lty=2, main="", lwd=3)
            
            grid()

            lines(myNA2$time, myH2, col=1, lwd=2)

            lines(myNA2$time, (myH2+1.96*mysig), col=1, lty=3)

            lines(myNA2$time, (myH2-1.96*mysig), col=1, lty=3)

            legend("topleft", c("Observed", "Expected"), col=c(1,2), lty=c(1,2), bty="n")

        }
        if(y==2){
            ## martingale plot

            eventtime<-myNA2$time

            myfillin<-sapply(eventtime, function(ind) which(x@crHt[,1]==ind) )

            nz<-nrow(x@crHt) 
            
            myobs.NA<-rep(0, nz)

            myobs.sig.NA<-rep(0, nz)
            
            for(idy in 1:(length(eventtime)-1)){
                
                myobs.NA[myfillin[idy]:(myfillin[idy+1]-1)] <- myH2[idy]               

                myobs.sig.NA[myfillin[idy] : (myfillin[idy+1]-1)] <- mysig[idy]               
            }
            
            myobs.NA[max(myfillin):length(myobs.NA)]<-max(myH2)

            myobs.sig.NA[max(myfillin):length(myobs.sig.NA)]<-max(mysig)

            ##martingale
            mytimei<-floor(myNA2$time)
            
            mytimer<-myNA2$time-mytimei
            
            mymodH2<-c(0, x@crHt[,1+m])

            plotdta.mart<-cbind(x@crHt[,1], (myobs.NA/x@crHt[,1+m]))[2:length(myobs.NA),]

            plotylim<-c(-0.5,2)
            
            plot(plotdta.mart, type="l", ylab="O/E (Cumulative Hazard)", xlab="Time (y)" ,   main="", ylim=plotylim)
            
            abline(h=1, lty=2, col=2,lwd=3)

            abline(h=0, lty=3, col="gray", lwd=3)

            grid()
            

            lines(x@crHt[,1], (myobs.NA+1.96*myobs.sig.NA)/x@crHt[,1+m],  lty=3, col=1)

            lines(x@crHt[,1], (myobs.NA-1.96*myobs.sig.NA)/x@crHt[,1+m],  lty=3, col=1)


            
        }
    }
)



##Simulate data 
#maxT <- 20
#haz1 <-  (1:20)/100
#haz2 <-  rep(2,20)/100
#myn<-200

#mysumdta <- fn.spmsim(haz1, haz2, myn)
#myspm<-survpm(mysumdta, m=2, maxT, ctime=TRUE)

#par(mfrow=c(2,2))
#plot(myspm); plot(myspm,2)
#plot(myspm,1,2); plot(myspm,2,2)
#summary(myspm)
