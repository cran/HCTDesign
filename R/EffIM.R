
#' @title Monitoring the trial at interim looks for a trial with efficacy monitoring only
#' @description Calculates one-sided efficacy boundary values at the observed number of events.
#' @details The number of events have to be entered sequentially. See example.
#' @param d2 vector of number of events at which you want to monitor the trial.
#' @param dmax maximum number of events in the experimental group calculated from design function.
#' @param alpha type I error.
#' @param beta type II error.
#' @param last.look logical which indicates whether the current look is the last look or not. Default is FALSE. If true, the post hoc power is calculated.
#' @param d1 total number of events in the historical control group.
#' @param opt type of spending function: "OBF", "Gamma", "Rho" or "Pocock". Default is "OBF".
#' @param param Parameter for "gamma family" or rho family. Default value is 4.
#' @param etam value of the drift parameter obtained from design function.
#' @author  Tushar Patni, Yimei Li, Jianrong Wu, and Arzu Onar-Thomas.
#' @return A list containing efficacy boundary values along with the p-values and transformed information time for the current look. Post-hoc power is also calculated in case of early stopping of the trial.
#' @examples
#' #Interim look for the trial when the number of events is 13(first look).
#' gg<-EffIM(c(13),dmax=57,alpha=0.05,beta=0.1,etam=3.0726,d1=65,opt="OBF",last.look=FALSE)
#' #Interim look for the trial when the number of events is 35(second look).
#' gg<-EffIM(c(13,35),dmax=57,alpha=0.05,beta=0.1,etam=3.0726,d1=65,opt="OBF",last.look=FALSE)
#' @import stats
#' @references
#' \insertRef{doi:10.1002/pst.1756}{HCTDesign}
#' @references
#' \insertRef{doi:10.1080/10543406.2019.1684305}{HCTDesign}
#' @importFrom diversitree set.defaults
#' @export


EffIM<-function(d2,dmax,last.look,d1,etam,alpha,beta,opt,param){
  ti=d2/dmax
  if (last.look==FALSE) {
    k<-length(ti)
    if (k==1){
      R<-dmax/d1
      ts<-(1+R)*ti/(1+R*ti)

      if (opt=="OBF"){
        alpha1<-2-2*pnorm(qnorm(1-alpha/2)/sqrt(ts))
        beta1<-2-2*pnorm(qnorm(1-beta/2)/sqrt(ts))
      }
      if (opt=="Gamma") {
        gamma=-(param)
        alpha1=alpha*(1-exp(gamma*(ts)))/(1-exp(gamma))
        beta1=beta*(1-exp(gamma*(ts)))/(1-exp(gamma))
      }
      if(opt=="Rho"){
        rho=param
        alpha1=alpha*(ts)^rho
        beta1=beta*(ts)^rho
      }
      if(opt=="Pocock"){
        alpha1=alpha*log(1+(exp(1)-1)*ts)
        beta1=beta*log(1+(exp(1)-1)*ts)
      }

      ub<-qnorm(1-alpha1)
      lb<-NA
      p1<-NULL
      p2<-NULL
      p1[1]<-pnorm(ub[length(ti)],lower.tail = F)
      p2[1]<-NA
      k1<-NULL
      k2<-NULL
      k1[1]<-pnorm(ub[1],lower.tail = F)
      k2[1]<-NA
      posthoc<-NA
    }

    else {
      R<-dmax/d1
      alpha1=beta1=NULL
      alpha1[1]=beta1[1]=0
      ts=NULL

      fx1<-function(x,ub,covm,tprob){#fx1 for solving lower bound# 310
        kn=length(ub)
        lb=rep(-Inf,kn)
        pmv=mvtnorm::pmvnorm(lower=c(lb,x),upper=c(ub,Inf),sigma=covm)[1]
        tprob-pmv}

      fx2<-function(x,lb,etam,ts,covm,tprob){#fx2 for solving upper bound# 315
        kn=length(lb)
        ub=rep(Inf,kn)
        lmean=etam*sqrt(ts[1:(kn+1)])
        pmv=mvtnorm::pmvnorm(lower=c(lb,-Inf),upper=c(ub,x),mean=lmean,sigma=covm)[1]
        tprob-pmv}
      covmatrix=matrix(rep(0,(length(ti)+1)*(length(ti)+1)),ncol=length(ti)+1,byrow=T)

      for (i in 1:length(ti)) {
        ts[i]=(1+R)*ti[i]/(1+R*ti[i])

        if(opt=="OBF") {
          alpha1[i+1]=2-2*pnorm(qnorm(1-alpha/2)/sqrt(ts[i]))
          beta1[i+1]=2-2*pnorm(qnorm(1-beta/2)/sqrt(ts[i]))
        }
        if (opt=="Gamma"){
          gamma=-(param)
          alpha1[i+1]=alpha*(1-exp(gamma*(ts[i])))/(1-exp(gamma))
          beta1[i+1]=beta*(1-exp(gamma*(ts[i])))/(1-exp(gamma))
        }
        if(opt=="Rho"){
          rho=param
          alpha1[i+1]=alpha*(ts[i])^rho
          beta1[i+1]=beta*(ts[i])^rho}
        if(opt=="Pocock"){
          alpha1[i+1]=alpha*log(1+(exp(1)-1)*ts[i])
          beta1[i+1]=beta*log(1+(exp(1)-1)*ts[i])}
      }

      tij=c(0, ts)
      for (i in 1:(length(ti)+1)) {
        for (j in 1:(length(ti)+1)) {
          covmatrix[i,j]=min(tij[i],tij[j])/sqrt(tij[i]*tij[j])}}
      ub=lb=rep(0,length(ti))
      ub[1]=qnorm(1-alpha1[2])
      ubi=NULL

      for (i in c(2:length(ti))){
        ubi=c(ubi,ub[i-1])
        ub[i]=uniroot(fx1,lower=-10,upper=10,ub=ubi,
                      covm=covmatrix[2:(i+1),2:(i+1)],tprob=alpha1[i+1]-alpha1[i])$root}

      lb[1]=NA
      #lbi=NULL
      #for (i in 2:length(ti)) {
       # lbi=c(lbi,lb[i-1])
      #  lb[i]=uniroot(fx2,lower=-10,upper=10,lb=lbi,etam=etam,ts=ts,
                  #    covm=covmatrix[2:(i+1),2:(i+1)],tprob=beta1[i+1]-beta1[i])$root}

      p1<-NULL
      p2<-NA
      p1<-pnorm(ub[length(ti)],lower.tail = F)
      #p2<-pnorm(lb[length(ti)],lower.tail = F)

      k1<-NULL
      k2<-NA
      k1[1]<-pnorm(ub[1],lower.tail = F)
      #k2[1]<-pnorm(lb[1],mean = etam*sqrt(ts[1]))
      for (i in 2:length(ub)){
        k1[i]<-mvtnorm::pmvnorm(lower=c(rep(-Inf,length(lb[1:(i-1)])),ub[i]),upper=c(ub[1:(i-1)],Inf),sigma=covmatrix[2:(i+1),2:(i+1)])[1]
        #k2[i]<-mvtnorm::pmvnorm(lower=c(lb[1:i-1],-Inf),upper=c(ub[1:i-1],lb[i]),sigma=covmatrix[2:(i+1),2:(i+1)],mean = etam*sqrt(ts[1:i]))[1]
      }
      posthoc<-NA
    }
  }

  else{

    R<-dmax/d1
    alpha1=beta1=NULL
    alpha1[1]=beta1[1]=0
    ts=NULL
    covmatrix=matrix(rep(0,(length(ti)+1)*(length(ti)+1)),ncol=length(ti)+1,byrow=T)

    fx1<-function(x,ub,covm,tprob){#fx1 for solving lower bound# 310
      kn=length(ub)
      lb=rep(-Inf,kn)
      pmv=mvtnorm::pmvnorm(lower=c(lb,x),upper=c(ub,Inf),sigma=covm)[1]
      tprob-pmv}

    fx2<-function(x,lb,etam,ts,covm,tprob){#fx2 for solving upper bound# 315
      kn=length(lb)
      ub=rep(Inf,kn)
      lmean=etam*sqrt(ts[1:(kn+1)])
      pmv=mvtnorm::pmvnorm(lower=c(lb,-Inf),upper=c(ub,x),mean=lmean,sigma=covm)[1]
      tprob-pmv}

    for (i in 1:length(ti)) {
      ts[i]=(1+R)*ti[i]/(1+R*ti[i])
      if(opt=="OBF") {
        alpha1[i+1]=2-2*pnorm(qnorm(1-alpha/2)/sqrt(ts[i]))
        beta1[i+1]=2-2*pnorm(qnorm(1-beta/2)/sqrt(ts[i]))
      }
      if (opt=="Gamma") {
        gamma=-(param)
        alpha1[i+1]=alpha*(1-exp(gamma*(ts[i])))/(1-exp(gamma))
        beta1[i+1]=beta*(1-exp(gamma*(ts[i])))/(1-exp(gamma))
      }
      if(opt=="Rho"){
        rho=param
        alpha1[i+1]=alpha*(ts[i])^rho
        beta1[i+1]=beta*(ts[i])^rho}
      if(opt=="Pocock"){
        alpha1[i+1]=alpha*log(1+(exp(1)-1)*ts[i])
        beta1[i+1]=beta*log(1+(exp(1)-1)*ts[i])}

    }
    tij=c(0, ts)
    alpha1[length(ti)+1]<-alpha
    if(length(ti)==1){
      ub<-qnorm(1-alpha1[2])
      lb<-ub
      p1<-NULL
      p2<-NULL
      p1[1]<-pnorm(ub,lower.tail = F)
      p2[1]<-pnorm(lb,lower.tail = F)
      k1<-NULL
      k2<-NULL
      k1[1]<-pnorm(ub,lower.tail = F)
      k2[1]<-pnorm(lb,mean = etam*sqrt(ts[1]))
      posthoc<-1-k2}

    else{
    for (i in 1:(length(ti)+1)) {
      for (j in 1:(length(ti)+1)) {
        covmatrix[i,j]=min(tij[i],tij[j])/sqrt(tij[i]*tij[j])}}
    ub=lb=rep(0,length(ti))
    ub[1]=qnorm(1-alpha1[2])
    ubi=NULL
    for (i in c(2:length(ti))){
      ubi=c(ubi,ub[i-1])
      ub[i]=uniroot(fx1,lower=-10,upper=10,ub=ubi,
                    covm=covmatrix[2:(i+1),2:(i+1)],tprob=alpha1[i+1]-alpha1[i])$root}

    #lb[1]=qnorm(beta1[2])+etam*sqrt(ts[1])
    #lbi=NULL
    #for (i in 2:length(ti)) {
     # lbi=c(lbi,lb[i-1])
      #lb[i]=uniroot(fx2,lower=-10,upper=10,lb=lbi,etam=etam,ts=ts,
                 #   covm=covmatrix[2:(i+1),2:(i+1)],tprob=beta1[i+1]-beta1[i])$root}
    lb[1]<-ub[length(ti)]

    p1<-NULL
    p2<-NULL
    p1[1]<-pnorm(ub[length(ti)],lower.tail = F)
    p2[1]<-pnorm(lb[1],lower.tail = F)


    k1<-NULL
    k2<-NULL
    k1[1]<-pnorm(ub[1],lower.tail = F)
    ii<-length(ts)
    k2[1]<-pnorm(lb[1],mean = etam*sqrt(ts[length(ts)]))
    for (i in 2:length(ub)){
      k1[i]<-mvtnorm::pmvnorm(lower=c(rep(-Inf,length(lb[1:(i-1)])),ub[i]),upper=c(ub[1:(i-1)],Inf),sigma=covmatrix[2:(i+1),2:(i+1)])[1]
      #k2[i]<-mvtnorm::pmvnorm(lower=c(lb[1:i-1],-Inf),upper=c(ub[1:i-1],lb[i]),sigma=covmatrix[2:(i+1),2:(i+1)],mean = etam*sqrt(ts[1:i]))[1]
    }
    if(ts[length(ts)]!=1){
    posthoc<-1-mvtnorm::pmvnorm(lower = c(rep(-Inf,length(ub))),upper = c(ub),mean = etam*sqrt(ts[1:length(ti)]),sigma = covmatrix[2:(length(ti)+1),2:(length(ti)+1)])[1]
    }
    else{posthoc<-NA}}

  }


  ans=list(Efficacy=data.frame("Efficacy boundary in z-score scale"=ub[length(ti)],"Efficacy boundary in p-value scale"=p1,check.names = F),Futility=data.frame("Futility boundary in z-score scale"=lb[1],"Futility boundary in p-value scale"=p2,check.names = F),"Information time"=ts[length(ti)],post_hoc_power=posthoc)
  return(ans)


}

EffIM<-set.defaults(EffIM,opt="OBF",last.look=FALSE,param=4)
###########################################################



