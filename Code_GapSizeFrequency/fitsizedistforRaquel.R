################################################
fitdbhlik=function(nstem,mindbh,maxdbh,fitfcn)
# nstem is a vector of the number of trees in each size class, in order, including zeros
# mindbh is a vector of the minimum dbh of each size class, in the same order
# maxdbh is a vector of the maximum dbh of each size class
# fitfcn is a character string specifying the function to fit: "weib" for weibull, "exp" for exponential, "pow" for power
# so for example, you could pass nstem=c(10,9,10,8,3,1,0,1), mindbh=seq(9.5,16.5,by=1), maxdbh=mindbh+1
# if the data have dbh values of 10.0, 10.1, etc. with a precision of 0.1, then set up the size class boundaries to be 9.95, 10.05, 10.15, etc.
# that is, the size classes are the width of the precision, and are centered on the measured values
# fits a specified distribution fitfcn ("weib"=Weibull,"pow"=power-law)
# to a histogram of number of stems in different size classes
# using true maximum likelihood fitting
############################################
{
  totstems=sum(nstem)
  if (fitfcn=="weib" || fitfcn=="3par") {
    fit=optim(par=c(0.1,1),fn=errordbhlik,fitfcn=fitfcn,dbhmin=mindbh,dbhmax=maxdbh,n=nstem,totstems=totstems,method="SANN")
    fit=optim(par=fit$par,fn=errordbhlik,fitfcn=fitfcn,dbhmin=mindbh,dbhmax=maxdbh,n=nstem,totstems=totstems)
  }
  else if (fitfcn=="pow") {  # need to use a different function for 1-dimensional optimization
    fit=optimize(f=errordbhlik,interval=c(-1,10),fitfcn=fitfcn,dbhmin=mindbh,dbhmax=maxdbh,n=nstem,totstems=totstems)
  }
  else if (fitfcn=="exp") {  # need to use a different function for 1-dimensional optimization
    initpar=0.1
    probexpected=calcdbhhist(mindbh,maxdbh,fitfcn,initpar)
    minprob=min(probexpected)
    while(minprob==0) {
       initpar=initpar/2
       probexpected=calcdbhhist(mindbh,maxdbh,fitfcn,initpar)
       minprob=min(probexpected)
    }
    fit=optimize(f=errordbhlik,interval=c(-0.5*initpar,2*initpar),fitfcn=fitfcn,dbhmin=mindbh,dbhmax=maxdbh,n=nstem,totstems=totstems)
  }

  if (fitfcn=="weib")
    dbhfit=data.frame(par1=fit$par[1],par2=fit$par[2],loglike=fit$value,convergerror=fit$convergence,iterations=fit$counts[[1]],ksstat=NA)
  else # if (fitfcn=="pow")
    dbhfit=data.frame(par1=fit$minimum,par2=NA,loglike=fit$objective,convergerror=NA,iterations=NA,ksstat=NA)

  return(dbhfit)
} # end fitdbhlik


#############################################
errordbhlik=function(param,fitfcn,dbhmin,dbhmax,n,totstems)
# calculate likelihood of size distribution data
# given what is expected under fitfcn type dbh dist function with parameter values param
# size distribution data in the actual dbhdistribution histogram as represented by
# dbhmin, dbhmax, and n (number of stems between dbhmin and dbhmax)
{
  probexpected=calcdbhhist(dbhmin,dbhmax,fitfcn,param)

  toterr=-sum(n*log(probexpected))

  return(toterr)
} # end errordbhlik


##################################################
calcdbhhist=function(dbhmin,dbhmax,fitfcn,param)
# calculate the probability that stems occur in each size class (dbhmin-dbhmax)
# given a particular distribution type and its parameters
{
  mindbh=min(dbhmin)
  maxdbh=max(dbhmax)
  if (fitfcn=="weib") { # param[1] is the shape, param[2] is the scale
    expectedn=pweibull(dbhmax,param[1],param[2])-pweibull(dbhmin,param[1],param[2])
    totexpected=pweibull(maxdbh,param[1],param[2])-pweibull(mindbh,param[1],param[2])
  }
  else if (fitfcn=="pow") {  # param[1] is the slope of the log-log size distribution
#    if (param[1]>1) {
#        expectedn=pexp(log(dbhmax),(param[1]-1))-pexp(log(dbhmin),(param[1]-1))
#        totexpected=pexp(log(maxdbh),param[1]-1)-pexp(log(mindbh),param[1]-1)
#    }
#    else {
        expectedn=dbhmax^(1-param[1])-dbhmin^(1-param[1])
        totexpected=maxdbh^(1-param[1])-mindbh^(1-param[1])
#    }
  }
  else if (fitfcn=="exp") {  # param[1] is the parameter of a negative exponential size distributionb
    expectedn=pexp(dbhmax,param[1])-pexp(dbhmin,param[1])
    totexpected=pexp(maxdbh,param[1])-pexp(mindbh,param[1])
  }
  else if (fitfcn=="3par") {  #### Christian's integral
    j=seq(0.005,0.995,by=0.01)
    nbins=length(dbhmin)
    expectedn=rep(NA,nbins)
    for (i in 1:nbins)
      expectedn[i]=sum(prob3parfcn(dbhmin[i]+j(dbhmax[i]-dbhmin[i])),param)
    totexpected=sum(expectedn)
  }
  expectedn=expectedn/totexpected  # normalize to account for fact that no data below 10 mm, etc.
  return(expectedn)
} # end calcdbhhist

