
#' @title Sample size in terms of number of subjects in the experimental group
#' @description Calculates the total number of subjects for the experimental group using the total number of events(d2max:the output from design functions) and the estimated failure probability based on the person level historical control data and proportional hazard assumption.
#' @param time event time vector from person level historical control data.
#' @param event numeric vector indicating the status of event from person level historical control data.
#' @param d2max maximum number of events in the experimental group calculated from the design function.
#' @param opt the method of fitting survival curve-"log_normal" or "KM" (log-normal or Kaplan Meier). Default is "KM".
#' @param event_ind numeric value indicating the occurrence of event.
#' @param delta hazard ratio.
#' @param ta enrollment time.
#' @param tf follow-up time.
#' @author Tushar Patni, Yimei Li, Jianrong Wu, and Arzu Onar-Thomas.
#' @return  Returns the value of sample size.
#' @examples
#' time<-c(20,65,12,50,58,65,45,44)
#' event<-c(1,0,0,0,1,1,1,1)
#' d2max=57
#' gg<-SM(time,event,d2max,opt="log_normal",ta=4,tf=3,delta=0.57,event_ind=1)
#' @import stats
#' @importFrom survival Surv survfit
#' @importFrom flexsurv flexsurvreg
#' @references
#' \insertRef{doi:10.1002/pst.1756}{HCTDesign}
#' @references
#' \insertRef{doi:10.1080/10543406.2019.1684305}{HCTDesign}
#' @importFrom Rdpack reprompt
#' @importFrom diversitree set.defaults
#' @export


SM<-function(time,event,d2max,opt,event_ind,ta,tf,delta) {
  if (opt=="log_normal") {
    s1<-Surv(time, event == event_ind)
    model<-flexsurvreg(s1 ~ 1, dist="lognormal" ) # fit the log-normal model #
    summary(model)
    meanlog=model$res[1,1]
    sdlog=model$res[2,1]
    lognormal=function(t){1-pnorm((log(t)-meanlog)/sdlog)}
    S2=function(t){lognormal(t)^delta}
    p1new=1-1/ta*integrate(S2,tf,tf+ta)$value
    p1new
    samplesize=d2max/p1new
    samplesize=ceiling(samplesize)
  }
  else if (opt=="KM") {

    SurvObj <- Surv(time, event==event_ind)
    model<- survfit(SurvObj ~ 1, conf.type = "log-log")
    p0<-c(1, summary(model)$surv)   # KM survival probability #
    t0<-c(0, summary(model)$time)   # ordered failure times #
    outKM<-data.frame(t0=t0,p0=p0)
    KM<-function(t){
      t0=outKM$t0; p0=outKM$p0; k=length(t0)
      if ( t<0) {ans<-0}
      if (t>=t0[k] ) {S0=p0[k]}
      for (i in 1:(k-1)){
        if (t>=t0[i] & t<t0[i+1]) {S0=p0[i]}}
      return(S0)}
    S0=function(t){ans=KM(t); return(ans)}
    S1=function(t){ans=KM(t)^delta; return(ans)}
    b1=(ta/6)*(S1(tf)+4*S1(0.5*ta+tf)+S1(ta+tf))
    p1=1-b1/ta
    samplesize=ceiling(d2max/p1)
  }

  return(data.frame(Samplesize=samplesize))

}

SM<-set.defaults(SM,opt="KM")
############################################


