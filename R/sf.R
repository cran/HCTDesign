#' @title Log rank test for non-inferiority trial
#' @description Calculates the score function of the log rank test for non-inferiority trial
#' @param event event time vector from person level trial data.
#' @param status numeric vector indicating the status of event from person level trial data.
#' @param delta0 Non-inferiority margin.
#' @param group group string vector indicating the assignment of patients into control or experimental group.
#' @param experiment name of experimental group as character string.
#' @param control name of control group as character string.
#' @author Tushar Patni, Yimei Li and Jianrong Wu.
#' @return  Returns the value of score statistic.
#' @examples
#' time<-c(20,65,12,50,58,65,45,44)
#' event<-c(1,0,0,0,1,1,1,1)
#' group<-c(rep("exp",4),rep("cont",4))
#' gg<-sf(event=time,status=event,delta0=1.3,group=group,experiment="exp",control="cont")
#' @import stats
#' @importFrom Rdpack reprompt
#' @importFrom diversitree set.defaults
#' @export


sf<-function(event,status,delta0,group,experiment,control){


HR<-delta0

x1<-event[group==control]
s1<-status[group==control]

x2<-event[group==experiment]
s2<-status[group==experiment]

x<-list(x1,x2)
s<-list(s1,s2)

em<-c(length(x1),length(x2))

###LHS
num<-vector()
den<-vector()
frc<-vector()
for (i in 1:length(x1)){
  #NUmerator
  tm<-vector()
  for (j in 1:length(x2) ){
    tm[j]<-as.numeric(x2[j]>=x1[i])
  }
  num[i]<-s1[i]*sum(tm)

  #DEnominator
  sk1<-vector()
  for (o in 1:length(x1)){
    sk1[o]<-as.numeric(x1[o]>=x1[i])
  }
  sk1<-sum(sk1)

  sk2<-vector()
   for(l in 1:length(x2)){
     sk2[l]<-as.numeric(x2[l]>=x1[i])
   }
  sk2<-sum(sk2)

  den[i]<-sk1+ (HR* sk2)

  frc[i]<-num[i]/den[i]

}

#browser()

######RHS
num1<-vector()
den1<-vector()
frc1<-vector()

for (i in 1:length(x2)){
  #NUmerator
  tm<-vector()
  for (j in 1:length(x1) ){
    tm[j]<-as.numeric(x1[j]>=x2[i])
  }
  num1[i]<-s2[i]*sum(tm)

  #DEnominator
  sk1<-vector()
  for (o in 1:length(x1)){
    sk1[o]<-as.numeric(x1[o]>=x2[i])
  }
  sk1<-sum(sk1)

  sk2<-vector()
  for(l in 1:length(x2)){
    sk2[l]<-as.numeric(x2[l]>=x2[i])
  }
  sk2<-sum(sk2)

  den1[i]<-sk1+ HR* sk2

  frc1[i]<-num1[i]/den1[i]

}


w<-(HR*sum(frc))- sum(frc1)


############################################


##Calculation of variance
adj<-vector()

for (i in 1:2){

  frq<-vector()
  for (j in 1:em[i]){

    n1<-vector()
    for (e in 1:length(x1)){
      n1[e]<-as.numeric(x1[e]>=x[[i]][j])
    }

    n2<-vector()
    for (h in 1:length(x2)){
      n2[h]<-as.numeric(x2[h]>=x[[i]][j])
    }

    d1<-vector()
    for (q in 1:length(x1)){
      d1[q]<-as.numeric(x1[q]>=x[[i]][j])

    }

    d2<-vector()
    for (g in 1:length(x2)){
      d2[g]<-as.numeric(x2[g]>=x[[i]][j])
    }

    frq[j]<-(s[[i]][j]*sum(n1)*sum(n2))/(sum(d1)+(HR*sum(d2)))^2

  }

 # browser()

  adj[i]<-sum(frq)

}




sigma2<-HR*sum(adj)





ft<-w/sqrt(sigma2)


return(ft)

}









