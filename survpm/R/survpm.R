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
## SurvPM class definitions
########################################

####
## * crData format needed:
##     - id, time, cause, hazard (number discrete time units) cause 1, 2, ..., then covariates used in formula (if supplied)
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

        ## id, survival time and cause
        .Object@crData <- crData[,1:3] 

        ## expected cumulative hazard for each competing risk model
        .Object@crH <- as.matrix(crData[,4:(3+maxT*m)])

        ## number of competing risks (e.g. m=1 is survival analysis)
        .Object@m<-as.integer(m)

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
        allhaz<-calcH(.Object)

        ## Overal expected number of events assigned
        .Object@pmH<-allhaz[[1]]

        ## Overall expected number of events per time period assigned
        .Object@crhaz<-allhaz[[2]]

        ## Calculate overall calibration table based on exact Poisson etc
        .Object@resOEtot<-calcOE(.Object)

        ## Calculate time dependent expected cumulative hazard for each cause (i.e. through followup time, at each event or when model hazard may change (discrete time interval).
        .Object@crHt <-tdepcal(.Object)
        
        if(ctime==TRUE){
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

        ##formulas for calibration check models
        myformula<-vector("list", 5)

        ##1. calibration in the large
        myformula[[1]]<- as.formula("(mydtafit$J==idx) ~ offset(log (mydtafit[,2+idx]) )")

        ##2. calibration slope
        myformula[[2]] <- as.formula("(mydtafit$J==idx) ~ log (mydtafit[,2+idx])")

        ##3. overall calibration plus time
        myformula[[3]]<- as.formula( "(mydtafit$J==idx) ~ mydtafit$t + offset(log(mydtafit[,2+idx]))" )

        ##4. calibration slope plus time
        myformula[[4]]<- as.formula( "(mydtafit$J==idx) ~ mydtafit$t + log (mydtafit[,2+idx])" )

        ##5. subgroups / user adjustments plus time, if user supplied formula
        if(object@adjform!=""){
            
            myformula[[5]] <- as.formula( paste( c("(mydtafit$J==idx) ~ t +", object@adjform, "+ offset(log(mydtafit[,2+idx]))") ))
            
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
            
            for(idy in 1 : 4){

                thisglm<-  glm( myformula[[idy]] , mydtafit, family=poisson)

                myglm[[idx]][[idy]] <- fn.formatglm( thisglm )

            }

            ##IF adjmodel provided
            if(object@adjform != ""){

                myglm[[idx]][[5]] <- fn.formatglm( glm( myformula[[5]], mydtafit, family=poisson) )

            }

        }
            
        return(myglm)
        
    }
)


## Class method definition. Not intended for users (internal)
## Defines function that can calc cumulative hazards for cause-specific model projections over follow-up time (each event or when hazard function changes)
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

    ##number at risk
    atrisk <- numeric(nz)
    
    atrisk[1:lastrisk[1]] <- object@n
    
    for(idz in 1:(length(lastrisk)-1) ){
        atrisk[ (lastrisk[idz]+1) : (lastrisk[idz+1]) ] <- atrisk[lastrisk[idz]] - 1
    }

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
        
        ##number at risk
        atrisk <- numeric(nz)
        
        atrisk[1:lastrisk[1]] <- object@n
        
        for(idz in 1:(length(lastrisk)-1) ){
            atrisk[ (lastrisk[idz]+1) : (lastrisk[idz+1]) ] <- atrisk[lastrisk[idz]] - 1
        }
        
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
                
                print("**Calibration, overall, Poisson regression model**")
                thisout<-fn.formatCIm( object@resOEtime[[idx]][[1]] )
                print(thisout)

                print("**Calibration, slope, Poisson regression model**")
                thisout<-fn.formatCIm( object@resOEtime[[idx]][[2]] )
                thisoutm<-data.frame(c("Overall", "Slope"), thisout)
                colnames(thisoutm)<-c("Parameter", "Coefficient (95CI)")
                print(thisoutm)

                print("**Calibration, overall + time, Poisson regression model**")
                thisout<-fn.formatCIm( object@resOEtime[[idx]][[3]] )
                thisoutm<-data.frame(c("Overall", "Follow-up"), thisout)
                colnames(thisoutm)<-c("Parameter", "Coefficient (95CI)")
                print(thisoutm)

                print("**Calibration, overall + time, Poisson regression model**")
                thisout<-fn.formatCIm( object@resOEtime[[idx]][[4]] )
                thisoutm<-data.frame(c("Overall", "Follow-up", "Slope"), thisout)
                colnames(thisoutm)<-c("Parameter", "Coefficient (95CI)")
                print(thisoutm)

                if(object@adjform != ""){
                    print("**Calibration, bespoke, Poisson regression model**")
                    thisout<-fn.formatCIm( object@resOEtime[[idx]][[5]] )
                    thisrownames<-rownames(object@resOEtime[[idx]][[5]])
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
    definition=function(x,y,m=1,idx=1,...){

        ##NA fit
        plotdta<-data.frame(T=x@crData$t, CC= as.integer(x@crData$d==m))
        
        myNA2<-survfit(Surv(T,CC)~1, plotdta)
        
        myrt<-myNA2$n.event / myNA2$n.risk
        
        mysig<-sqrt(cumsum(myNA2$n.event / (myNA2$n.risk^2)))
        
        myH2<-cumsum(myrt)

        if(idx==1){
            ##Cum Haz plot
            plot(x@crHt[,1], x@crHt[,1+m]*100, xlim=c(0,x@maxT), ylim=c(0, max(c( max(100*(myH2+1.96*mysig)), max(x@crHt[,2]*100)))),  type="l", xlab="Time (y)", ylab="Cumulative Hazard (%)", col=2, lty=2, main="", lwd=3)
                                                                           
            grid()

            lines(myNA2$time, myH2*100, col=1, lwd=2)

            lines(myNA2$time, 100*(myH2+1.96*mysig), col=1, lty=3)

            lines(myNA2$time, 100*(myH2-1.96*mysig), col=1, lty=3)

            legend("topleft", c("Observed", "Expected"), col=c(1,2), lty=c(1,2), bty="n")

        }
        if(idx==2){
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

    
            plot(x@crHt[,1], (myobs.NA/x@crHt[,1+m]), type="l", ylab="O/E (Cumulative Hazard)", xlab="Time (y)" ,   main="", ylim=c(min((myobs.NA-1.96*myobs.sig.NA)/x@crHt[,1+m]), max((myobs.NA+1.96*myobs.sig.NA)/x@crHt[,1+m])))
            
            abline(h=1, lty=2)

            abline(h=0, lty=3, col=2)

            grid()
            

            lines(x@crHt[,1], (myobs.NA+1.96*myobs.sig.NA)/x@crHt[,1+m],  lty=3, col=1)

            lines(x@crHt[,1], (myobs.NA-1.96*myobs.sig.NA)/x@crHt[,1+m],  lty=3, col=1)


            
        }
    }
)


