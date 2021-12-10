################################################
doallbootstrapdbhfits=function(alldata,
  nreps=100, # number of bootstrap replicates - should be 1000 or so
  alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
  fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
  # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
  filestem="test", # this is the beginning of the names of the output files
  ddiv=seq(10.05,200.05,by=0.1)) # these are the divisions, in mm, between size classes
################################################
# This function calls fits to the size distribution with bootstrapping by subplot
################################################
{
  nfcn=length(fitfcn)
  nfitdiv=length(ddiv)-1
  mindbh=ddiv[1:nfitdiv]
  maxdbh=ddiv[2:(nfitdiv+1)]
  
  allbsdbhcounts=list()
  
    # make the bootstrapped dataset
  dbhhists=makebootstrapdbhhist(alldata,ddiv, nreps)
  
  # compile 95% CIs on the numbers of stems in each size category
  allbsdbhcounts[[i]]=bootstrapdbhcounts(mindbh,maxdbh,dbhhists$truehist,
                                         dbhhists$bshistmatrix,alpha)

  statdbhfit=data.frame(nreps=nreps,n=sum(dbhhists$truehist),nfitdiv=nfitdiv,
                        mindbh=min(mindbh), maxdbh=max(maxdbh))
    
  for (j in 1:nfcn) { # do fits for all types of fits ("weib","pow","exp")
    # do the fits
    thisdbhfit=bootstrapdbhfit(dbhhists$truehist,dbhhists$bshistmatrix,
                               mindbh,maxdbh,fitfcn[j],alpha)
    names(thisdbhfit)=paste(fitfcn[j],names(thisdbhfit),sep="")
    
    # put results into a data fraome
    if (j==1)
      thesedbhfits=cbind(statdbhfit,thisdbhfit)
    else
      thesedbhfits=cbind(thesedbhfits,thisdbhfit)
  }
  alldbhfits=thesedbhfits

  outfitfn=paste(datadir,filestem,"sizedistbsfit.txt",sep="")
  write.table(alldbhfits,file=outfitfn,sep="\t",row.names=F)
  print("Writing size distribution fitted parameters & their CIs to the file")
  print(outfitfn)
  names(allbsdbhcounts)=paste(site,census,sep="")
  outcountfn=paste(datadir,filestem,"sizedistbscount.rdata",sep="")
  save(allbsdbhcounts,file=outcountfn)
  print("Writing size distribution counts & CIs for the fitting size classes to the file")
  print(outcountfn)
  return(alldbhfits)
} # end doallbootstrapdbhfits




################################################
makebootstrapdbhhist=function(alldata,ddiv=seq(10.05,200.05,by=0.1),nbootstrap)
  # make a list of dataframes giving counts of stems in different size classes according to ddiv
  ################################################
# alldata is a dataframe with columns 
# "dbh" giving the dbh in cm and 
# "quadnum" giving a quadrat number or identifier, where quadrats are the unit for bootstrapping
# "nbootstrap" is the number of bootstrap iterations to run
# bootstraps are samples of quadrats with replacement
{
  inc=alldata$dbh>=min(ddiv) & alldata$dbh<max(ddiv) & !is.na(alldata$dbh)& !is.na(alldata$quadnum) 
  
  quadsums=countbyquadndivsimp(alldata[inc,c("dbh","quadnum")],ddiv)

  truecounts=apply(quadsums,2,sum,na.rm=T)
  
  ngrid=dim(quadsums)[1]
  dbhhistmatrix=matrix(0,nrow=length(ddiv)-1,ncol=nbootstrap)
  for (j in 1:nbootstrap) {
    whichgrids=sample(ngrid,ngrid,replace=T)
    dbhhistmatrix[,j]=apply(quadsums[whichgrids,],2,sum,na.rm=T)
    if(j%%50==0) cat("finished boot ", j, "\n")
  }
  nclass <- length(ddiv)-1
  rownames(dbhhistmatrix)=paste(ddiv[1:nclass],"-",ddiv[2:(nclass+1)],sep="")
  dbhhistmatrix[is.na(dbhhistmatrix)]=0
  output=list(truecounts,dbhhistmatrix)
  names(output)=c("truehist","bshistmatrix")
  return(output)
} # end makebootstrapdbhhist


################################################
countbyquadndivsimp=function(alldata,ddiv)
  # a generic function for counting data by quadrat quadnum 
  # and size classes ddiv
  # alldata is a dataframe with columns "dbh" and "quadnum" defined as above
  # NOTE: if there are no stems in an xy quadrat, 
  # then it will not be a column in the output file
  ################################################
{
  inc=!is.na(alldata$dbh)&alldata$dbh>=min(ddiv)&alldata$dbh<=max(ddiv)&!is.na(alldata$quadnum)
  usedata=alldata[inc,]
  ndclass=length(ddiv)-1
  dbhclass=as.numeric(cut(usedata$dbh,breaks=ddiv,right=F,labels=1:ndclass))
  dn=tapply(usedata$dbh,list(usedata$quadnum,dbhclass),length)
  dmtch=match((1:ndclass),colnames(dn))
  qs=sort(unique(usedata$quadnum))
  xmtch=match(qs,rownames(dn))
  dn=dn[xmtch,dmtch]
  colnames(dn)=paste(ddiv[1:ndclass],"-",ddiv[2:(ndclass+1)],sep="")
  dn[is.na(dn)]=0
  return(dn)
} # end countbyquadndivsimp



################################################
bootstrapdbhcounts=function(mindbh,maxdbh,truehist,bshistmatrix,alpha=c(0.05,0.01)) {
  quants=apply(bshistmatrix,1,quantile,probs=c(alpha[1]/2,1-alpha[1]/2,alpha[2]/2,1-alpha[2]/2))
  quantdf=as.data.frame(t(quants))
  
  names(quantdf)=c(paste("lo",100*alpha[1],sep=""),paste("hi",100*alpha[1],sep=""),
                   paste("lo",100*alpha[2],sep=""),paste("hi",100*alpha[2],sep=""))
  bsdbhcounts=data.frame(mindbh=mindbh,maxdbh=maxdbh,n=truehist)
  bsdbhcounts=cbind(bsdbhcounts,quantdf)
  return(bsdbhcounts)
} # end bootstrapdbhcounts


