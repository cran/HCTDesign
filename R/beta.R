 # Helper function
fx2DF<-function(x,lb,etam,ts,covm,tprob){#fx2 for solving upper bound#
  kn=length(lb)
  ub=rep(Inf,kn)
  lmean=etam*sqrt(ts[1:(kn+1)])
  pmv=mvtnorm::pmvnorm(lower=c(lb,-Inf),upper=c(ub,x),mean=lmean,sigma=covm)[1]
  tprob-pmv}

# Helper function
findboundDF<-function(d2,d1,alpha,beta,k,option,param){
  ti<-k
  eta0=qnorm(1-alpha)+qnorm(1-beta)
  eta1=sqrt(2)*eta0
  etam=(eta0+eta1)/2
  R=d2/d1
  alpha1=beta1=NULL
  alpha1[1]=beta1[1]=0
  ts=NULL
  if(option=="OBF"){
    for (i in 1:length(ti)) {
      ts[i]=(1+R)*ti[i]/(1+R*ti[i])
      #alpha1[i+1]=2-2*pnorm(qnorm(1-alpha/2)/sqrt(ts[i]))
      beta1[i+1]=2-2*pnorm(qnorm(1-beta/2)/sqrt(ts[i]))}}
  if(option=="Gamma"){
    gamma=-(param)
    for (i in 1:length(ti)) {
      ts[i]=(1+R)*ti[i]/(1+R*ti[i])
      #alpha1[i+1]=alpha*(1-exp(gamma*(ts[i])))/(1-exp(gamma))
      beta1[i+1]=beta*(1-exp(gamma*(ts[i])))/(1-exp(gamma))}}
  if(option=="Rho"){
    rho=param
    for (i in 1:length(ti)) {
      ts[i]=(1+R)*ti[i]/(1+R*ti[i])
      #alpha1[i+1]=alpha*(ts[i])^rho
      beta1[i+1]=beta*(ts[i])^rho}}
  if(option=="Pocock"){
    for (i in 1:length(ti)) {
      ts[i]=(1+R)*ti[i]/(1+R*ti[i])
      #alpha1[i+1]=alpha*log(1+(exp(1)-1)*ts[i])
      beta1[i+1]=beta*log(1+(exp(1)-1)*ts[i])}}

  covmatrix=matrix(rep(0,(length(ti)+1)*(length(ti)+1)),ncol=length(ti)+1,byrow=T)
  tij=c(0, ts)
  for (i in 1:(length(ti)+1)) {
    for (j in 1:(length(ti)+1)) {
      covmatrix[i,j]=min(tij[i],tij[j])/sqrt(tij[i]*tij[j])}}
  ub=lb=rep(0,length(ti))
  ub[1]=qnorm(1-alpha)

  ctn=0
  flag2=0 #flag for convergence
  repeat{
    ctn=ctn+1
    flag=0 #flag for lb>ub
    etam=(eta0+eta1)/2
    lb[1]=qnorm(beta1[2])+etam*sqrt(ts[1])
    lbi=NULL
    for (i in 2:length(ti)) {
        lbi=c(lbi,lb[i-1])
        lb[i]=uniroot(fx2DF,lower=-10,upper=10,lb=lbi,etam=etam,ts=ts,
                      covm=covmatrix[2:(i+1),2:(i+1)],tprob=beta1[i+1]-beta1[i])$root
      }

       lb[length(ti)]=ub[1]
        pv=NULL
        pv[1]=pnorm(lb[1],mean=etam*sqrt(ts[1]))
        for (i in 2:length(ti)){
          lmean=etam*sqrt(ts[1:i])
          covm=covmatrix[2:(i+1),2:(i+1)]
          pv[i]=mvtnorm::pmvnorm(lower=c(lb[1:(i-1)],-Inf),upper=c(rep(Inf,length(lb[1:(i-1)])),lb[i]),
                                 mean=lmean,sigma=covm)[1]}
        betaK=sum(pv)
        if (betaK<beta) {eta1=etam} else {eta0=etam}
        if (abs(beta-betaK)<1e-05) {flag2=1}
        if(flag2==1) {break}

      }



  ans=list(lower=lb,upper=ub[1],etam=round(etam,4), ts=round(ts,4),flag,
           flag2,alpha=alpha,beta=beta1[2:length(beta1)],sig=covmatrix)
  return(ans)}


#' @title HCT design with interim monitoring for futiity only
#' @description The group sequential design for historical controlled survival outcome trials with futility boundaries only.
#' @param k vector of time fraction for all planned looks: k=c(1/3,2/3,1) if the three planned looks will be carried out at 1/3, 2/3 and all of the total events in the experiment arm.
#' @param alpha type I error.
#' @param beta type II error.
#' @param delta hazard ratio: hazard of experiment group over hazard of control group.
#' @param delta0 Non-inferiority margin.
#' @param d1 total number of events in the historical control group.
#' @param option type of spending function: "OBF", "Gamma", "Rho" or "Pocock". Default is "OBF.
#' @param param Parameter for Gamma family or Rho family. Default value is 4.
#' @param trial Type of trial: "Superiority" or "Non-inferiority". Default is "Superiority".
#' @author  Tushar Patni, Yimei Li and Jianrong Wu.
#' @return List of dataframes and vectors containing the details about the following: design of the trial which includes the number of looks and events;
#' details about futility and efficacy boundaries which include transformed information time at each look, cumulative beta and alpha respectively, p-values and crossing probabilities;
#' etam(drift parameter); d2max(maximum number of events in the experimental group); delta_used(hazard ratio used in the design).
#' @examples
#' #Sequential superiority trial for three equally spaced looks for OBF spending function.
#' gg<-FutDesign(k=c(0.3,0.6,1),alpha=0.05,beta=0.1,delta=0.57,d1=65,option="OBF",trial="Superiority")
#' @import stats
#' @import crayon
#' @references
#' \insertRef{doi:10.1002/pst.1756}{HCTDesign}
#' @references
#' \insertRef{doi:10.1080/10543406.2019.1684305}{HCTDesign}
#' @importFrom Rdpack reprompt
#' @importFrom diversitree set.defaults
#' @export


FutDesign<-function(k, alpha, beta, delta, d1, option,param,trial,delta0){
  intdelt<-delta
  j<-0
  while(d1>0) {
    if (trial=="Superiority"){
      temp1=round(exp(sqrt(1/d1*(qnorm(1-alpha)+qnorm(1-beta))^2)),2)
      temp2=round(exp(-sqrt(1/d1*(qnorm(1-alpha)+qnorm(1-beta))^2)),2)
    }

    else if (trial=="Non-inferiority") {
      temp1=round(exp(sqrt(1/d1*(qnorm(1-alpha)+qnorm(1-beta))^2))*delta0,2)
      temp2=round(exp(-sqrt(1/d1*(qnorm(1-alpha)+qnorm(1-beta))^2))*delta0,2)  }

    else{stop("Choose from Superiority or Non-inferiority")}


    if (delta<=temp1 & delta>=temp2) {
      delta<-ifelse(floor(delta)==1,temp1+0.01,temp2-0.01)
    }

    if (trial=="Superiority"){d2start=ceiling((log(delta)^2/(qnorm(1-alpha)+qnorm(1-beta))^2-1/d1)^(-1))}
    else if(trial=="Non-inferiority") {d2start=ceiling(((log(delta)-log(delta0))^2/(qnorm(1-alpha)+qnorm(1-beta))^2-1/d1)^(-1))}

    #if(d2start<0 & delta==ifelse(floor(delta)==1,temp1,temp2)){
     # delta<-ifelse(floor(delta)==1,delta+0.001,delta-0.001)
      #d2start=ceiling((log(delta)^2/(qnorm(1-alpha)+qnorm(1-beta))^2-1/d1)^(-1))}

    find=findboundDF(d2=d2start,d1=d1,alpha=alpha, beta=beta,k=k,
                   option=option,param=param)
    etam=find$etam

    ctn=0
    conv<-0
    repeat{
      ctn=ctn+1

      if (trial=="Superiority"){
        temp1=round(exp(sqrt(1/d1*(etam)^2)),2)
        temp2=round(exp(-sqrt(1/d1*(etam)^2)),2)
      }
      else if (trial=="Non-inferiority"){
        temp1=round(exp(sqrt(1/d1*(etam)^2))*delta0,2)
        temp2=round(exp(-sqrt(1/d1*(etam)^2))*delta0,2)
      }

      if (delta<=temp1 & delta>=temp2) {
        delta<-ifelse(floor(delta)>=1,temp1+0.01,temp2-0.01)
      }

      if (trial=="Superiority"){d2=ceiling(1/(log(delta)^2/etam^2-1/d1))}
      else if (trial=="Non-inferiority"){d2=ceiling(1/((log(delta)-log(delta0))^2/etam^2-1/d1))}


      find=findboundDF(d2=d2,d1=d1,alpha=alpha, beta=beta,k=k,option=option,param=param)
      etam=find$etam
      #temp1=round(exp(sqrt(1/d1*(etam)^2)),3)
      #temp2=round(exp(-sqrt(1/d1*(etam)^2)),3)
      #if (delta<temp1 & delta>temp2) {
       # delta<-ifelse(floor(delta)==1,temp1+0.001,temp2-0.001)
      #}

      lb=find$lower
      ub=find$upper
      ts=find$ts
      if(ctn==1){
        q=d2
      }
      if(d2==q){
        conv<-conv+1
      }
      else{
        q<-d2
        conv<-1
      }
      if (conv>=5){break}
      if (ctn>50) {break}}
    if (conv>=5){break}
    else {
      j<-j+1
      if(j>=3){stop("Use different value of Hazard ratio")}
      delta<-ifelse(floor(delta)==1,delta+0.01,delta-0.01)}
  }

  lb=round(lb,3)
  ub=round(ub,3)
  p1<-NULL
  p2<-NULL
  l1<-NULL
  l2<-NULL
  p1[1]<-pnorm(ub[1],lower.tail = F)
  l1[1]<-pnorm(ub[1],lower.tail = F,mean =etam*sqrt(ts[length(ts)]))
  p2[1]<-pnorm(lb[1],mean = etam*sqrt(ts[1]))
  l2[1]<-pnorm(lb[1])
  for (i in 2:length(lb)){
    #p1[i]<-mvtnorm::pmvnorm(lower=c(lb[1:i-1],ub[i]),upper=c(ub[1:i-1],Inf),sigma=find$sig[2:(i+1),2:(i+1)])[1]
    p2[i]<-mvtnorm::pmvnorm(lower=c(lb[1:(i-1)],-Inf),upper=c(rep(Inf,length(lb[1:(i-1)])),lb[i]),sigma=find$sig[2:(i+1),2:(i+1)],mean = etam*sqrt(ts[1:i]))[1]
    l2[i]<-mvtnorm::pmvnorm(lower=c(lb[1:(i-1)],-Inf),upper=c(rep(Inf,length(lb[1:(i-1)])),lb[i]),sigma=find$sig[2:(i+1),2:(i+1)])[1]
  }
  k1<-pnorm(ub[1],lower.tail = F)
  k2<-NULL
  for (i in 1:length(lb)){
    #k1[i]<-pnorm(ub[i],lower.tail = F)
    k2[i]<-pnorm(lb[i],lower.tail = F)
  }

  if(intdelt!=delta){cat(red$bold(paste("The final delta used for the design is",delta)))
    cat("\n\n")
  }

  ans=list(Design=data.frame(Looks=1:length(ts),Events=round(d2*k)),Efficacy=data.frame("Information time"=round(ts[length(ts)],4),"Cumulative alpha spent"=find$alpha,"Efficacy boundary in z-score scale"=ub[1],
                                                                                        "Efficacy boundary in p-value scale"=k1,"Boundary crossing probability under H0"=p1,check.names = F),
           Futility=data.frame("Information time"=round(ts,4),
                               "Cumulative beta spent"=find$beta,"Futility boundary in z-score scale"=lb,"Futility boundary in p-value scale"=k2,"Boundary crossing probability under H0"=l2,"Boundary crossing probability under H1"=p2,check.names = F),etam=etam,d2max=d2,delta_used=delta,trial=trial)

  return(ans)}

FutDesign<-set.defaults(FutDesign,option="OBF",param=4,trial="Superiority")
#################################################################
