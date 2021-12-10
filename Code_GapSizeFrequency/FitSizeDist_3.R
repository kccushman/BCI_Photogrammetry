
################################################
bootstrapdbhfit=function(truehist,bsdbhdists,mindbh,maxdbh,fitfcn,alpha=c(0.05,0.01)) {
  # This function bootstraps size distributions by subplots of dimension gridsize x gridsize meters
  # and then fits dbh distributions to each bootstrap and returns statistics of bootstrap
  # The fit is done using a straightforward maximum likelihood method:
  # calculate the probability of having a tree in a given interval given the parameters,
  # multiply the log of this probability by the number of trees in that interval,
  # and sum over all intervals.
  # fitfcn if the function to be fitted to the size distribution, and can be can be
  #    "weib"=Weibull (2-parameter) or
  #    "pow" = power-law (1-parameter because it's a true probability distribution)
  #    "exp" = negative exponential (1-parameter)
  ################################################
  nfitdiv=length(mindbh)
  nreps=dim(bsdbhdists)[2]
  # fits to the data for the whole plot (no bootstrapping)
  truefit=fitdbhlik(truehist,mindbh,maxdbh,fitfcn)
  
  catfittedn=sum(truehist)*calcdbhhist(mindbh,maxdbh,fitfcn,c(truefit$par1,truefit$par2))
  catreg=lm(log(catfittedn+1)~log(truehist+1))
  cats.r=summary(catreg)
  
  # fits to bootstrapped data
  bsfits=apply(bsdbhdists,2,fitdbhlik,mindbh=mindbh,maxdbh=maxdbh,fitfcn=fitfcn)
  bsfitpars=as.data.frame(t(matrix(as.numeric(as.data.frame(bsfits)),nrow=6,ncol=nreps)))
  names(bsfitpars)=names(bsfits[[1]])
  
  # put results into a data fraome
  dbhfits=data.frame(par1est=truefit$par1,
                     par2est=truefit$par2,
                     loglike=truefit$loglike,
                     convergerror=truefit$convergerror,
                     iterations=truefit$iterations,
                     catchisq=sum(((truehist-catfittedn)^2)/catfittedn)/nfitdiv,
                     catr2=cats.r$r.squared,
                     catadjr2=cats.r$adj.r.squared,
                     alpha1=alpha[1],
                     par1lo1=quantile(bsfitpars$par1,alpha[1]/2),
                     par1hi1=quantile(bsfitpars$par1,1-alpha[1]/2),
                     par1lo2=quantile(bsfitpars$par1,alpha[2]/2),
                     par1hi2=quantile(bsfitpars$par1,1-alpha[2]/2),
                     alpha2=alpha[2],
                     par2lo1=ifelse(fitfcn=="weib",quantile(bsfitpars$par2,alpha[1]/2),NA),
                     par2hi1=ifelse(fitfcn=="weib",quantile(bsfitpars$par2,1-alpha[1]/2),NA),
                     par2lo2=ifelse(fitfcn=="weib",quantile(bsfitpars$par2,alpha[2]/2),NA),
                     par2hi2=ifelse(fitfcn=="weib",quantile(bsfitpars$par2,1-alpha[2]/2),NA))
  return(dbhfits)
  
} # end bootstrapdbhfit