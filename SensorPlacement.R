sensorpl <- function(V,np,ttheta) {
  # Description: This function takes in input and returns the order of points selected
  #   using the method of Krause et al. (2008)
  # Input:
  #   V: the designmatrix where each row is a design point
  #   np: # of sensors to place
  #   ttheta: theta for the correlation function
  # Output:
  #   A.inds: vector of indices of the selected rows
  
  dd <- ncol(V) # dimensions of problem
  
  A.inds <- c() # Selected indices
  Abar.inds <- 1:(dim(V)[1]) # Unselected indices
  
  # This is to find first point
  deltas <- rep(NA,length(Abar.inds)) # We pick point with highest delta
  for (Abar.ind in 1:length(Abar.inds)) { # loop through all unselected and calculate delta for it
    Abar <- V # Points not selected
    yy <- Abar[Abar.ind,] # The point we are checking
    Abar.temp <- matrix(Abar[-Abar.ind,],nrow=nrow(Abar)-1,byrow=T) # All the remaining points
    #SAy <- apply(A,1,function(rrow){gaussian.corr(yy,rrow,ttheta)})
    SAbar.tempy <- apply(Abar.temp,1,function(rrow){gaussian.corr(yy,rrow,ttheta)})
      # Correlation of Abar.tempy with yy
    
    SAbarAbar.temp <- matrix(NA,length(Abar.inds)-1,length(Abar.inds)-1) # Correlation matrix for SAbar
    for (row.ind in 1:(length(Abar.inds)-1)) {
      SAbarAbar.temp[row.ind,] <- apply(Abar.temp,1,function(rrow){gaussian.corr(Abar.temp[row.ind,],rrow,ttheta)})
    }
    s2y <- 1 # I think this is what it should be, not sure actually
    delta <- 1/(s2y-t(SAbar.tempy)%*%solve(SAbarAbar.temp,SAbar.tempy))
    deltas[Abar.ind] <- delta # Update with the value of delta
  }
  #print(deltas)
  begin.with <- which.max(deltas) # Stores index to store first
  rm(yy,Abar.ind,Abar.temp,SAbar.tempy,SAbarAbar.temp,row.ind,s2y,deltas,delta,A.inds,Abar.inds)
  #   Delete everything, did this to make sure there was no contamination
  ### End of first point selection
  
  A.inds <- begin.with # selected indices
  Abar.inds <- setdiff(1:(dim(V)[1]),begin.with) # Unselected indices
  
  # Get rest of indices
  for(k in 2:np) {#while(k<np) { # place next sensor
    # get A corr matrix
    if(dd==1) # Need separate cases when 1-D vs 2+-D
      A <- matrix(V[A.inds],nrow=length(A.inds))
    else { 
      if(length(A.inds)==1)
        A <- matrix(V[A.inds,],nrow=1)
      else
        A <- V[A.inds,]
    }
    SAA <- matrix(NA,length(A.inds),length(A.inds)) # Correlation matrix for A
    for (row.ind in 1:length(A.inds)) {
      SAA[row.ind,] <- apply(A,1,function(rrow){gaussian.corr(A[row.ind,],rrow,ttheta)})
    }
    SAA.inv <- solve(SAA)
    # get Abar corr matrix BUT Abar needs to exclude y!!!
    if(dd==1)
      Abar <- matrix(V[-A.inds],nrow=length(Abar.inds),byrow=T)
    else
      Abar <- V[-A.inds,]
    
    # try new in for loop
    deltas <- rep(NA,length(Abar.inds))
    for (Abar.ind in 1:length(Abar.inds)) {
      yy <- Abar[Abar.ind,]
      Abar.temp <- matrix(Abar[-Abar.ind,],nrow=nrow(Abar)-1,byrow=T)
      SAy <- apply(A,1,function(rrow){gaussian.corr(yy,rrow,ttheta)})
      SAbar.tempy <- apply(Abar.temp,1,function(rrow){gaussian.corr(yy,rrow,ttheta)})
        
      SAbarAbar.temp <- matrix(NA,length(Abar.inds)-1,length(Abar.inds)-1)
      for (row.ind in 1:(length(Abar.inds)-1)) {
        SAbarAbar.temp[row.ind,] <- apply(Abar.temp,1,function(rrow){gaussian.corr(Abar.temp[row.ind,],rrow,ttheta)})
      }
      # option 1 is to find inverse
      #SAbarAbar.temp.inv <- solve(SAbarAbar.temp)
      #delta <- -(0-t(SAy)%*%SAA.inv%*%SAy)/(0-t(SAbar.tempy)%*%SAbarAbar.temp.inv%*%SAbar.tempy)
      # Changed it to solve system
      s2y <- 1
      delta <- (s2y-t(SAy)%*%SAA.inv%*%SAy)/(s2y-t(SAbar.tempy)%*%solve(SAbarAbar.temp,SAbar.tempy))
      deltas[Abar.ind] <- delta
    }
    # have deltas now
    
    #y.ind <- which.min(dMIs)
    ind.add.to.A <- Abar.inds[which.max(deltas)]
    if(length(ind.add.to.A)>1) ind.add.to.A <- ind.add.to.A[1]
    # add smallest to A
    A.inds <- c(A.inds,ind.add.to.A)
    Abar.inds <- Abar.inds[Abar.inds!=ind.add.to.A]
  }
  return(A.inds)
}

gaussian.corr <- function(xx1,xx2,ttheta) { # Gets the Gaussian correlation
  exp(-sum(abs(xx1-xx2)^2*ttheta))
}
if (F) {
  # 1D
  sensorpl(matrix(c(0,.2,.4,.6,.8,1),6,1),1,1)
  sensorpl(matrix(runif(6),6,1),1,1)
  set.seed(0)
  npts=10
  kpts=2
  V1 <- matrix(runif(npts),ncol=1)
  V1.inds <- sensorpl(V1,kpts,1)
  mrkr <- rep('o',npts)
  mrkr[V1.inds] <- as.character(1:length(V1.inds))
  plot(V1,col=((1:npts)%in%V2.inds + 1),pch=mrkr)
  stripchart(V1[-V1.inds,])
  stripchart(V1[V1.inds,],pch=as.character(1:length(V1.inds)),col='red',add=T,at = 1.1)
}
get.1D.sensors <- function(reps,npts,kpts,seed=0,theta=1,LHS=F) {
  if(LHS) require(lhs)
  dat.store <- list()
  dat.store.inds <- list()
  set.seed(seed)
  for(ii in 1:reps) {
    if (LHS)
      V1 <- matrix(sort(lhs::maximinLHS(npts,1)),ncol=1)
    else
      V1 <- matrix(sort(runif(npts)),ncol=1)
    V1.inds <- sensorpl(V1,kpts,theta)
    mrkr <- rep('o',npts)
    mrkr[V1.inds] <- as.character(1:length(V1.inds))
    dat.store[[ii]] <- V1
    dat.store.inds[[ii]] <- V1.inds
  }
  #plot(V1,col=((1:npts)%in%V2.inds + 1),pch=mrkr)
  stripchart(dat.store,col='white',xlab='Design points',ylab='Replicate number',main='Sensor placement in 1-D')
  for (ii in 1:reps){
    stripchart(dat.store[[ii]][-dat.store.inds[[ii]]],add=T,at=ii,pch=4);#browser()
    stripchart(dat.store[[ii]][dat.store.inds[[ii]]],add=T,at=ii,pch=4,col=2);#browser()
    for(jj in 1:(kpts+1))
      stripchart(dat.store[[ii]][dat.store.inds[[ii]][jj]],add=T,at=ii+.2*ifelse(ii==reps,-1,1),col=1+jj,pch=as.character(jj))
  }
  #stripchart(V1[V1.inds,],pch=as.character(1:length(V1.inds)),col='red',add=T,at = 1.1)
}
if (F) {
  #get.1D.sensors(5,8,2)
  get.1D.sensors(5,8,3,seed=3000,theta=1)
  get.1D.sensors(5,8,3,seed=3000,theta=10)
  get.1D.sensors(5,8,3,seed=3000,theta=10,LHS=T)
}
if (F) {
  # 2D
  set.seed(0)
  npts=30
  kpts=8
  V2 <- matrix(runif(npts*2),ncol=2)
  V2.inds <- sensorpl(V2,kpts,1)
  mrkr <- rep('o',npts)
  mrkr[V2.inds] <- as.character(1:length(V2.inds))
  plot(V2,col=((1:npts)%in%V2.inds + 1),pch=mrkr)
  plot(V2[-V2.inds,])
  points(V2[V2.inds,],pch=as.character(1:length(V2.inds)),col='red')
}
get.2D.sensors <- function(reps,npts,kpts,seed=0,theta=1,LHS=F) {
  if(LHS) require(lhs)
  set.seed(seed)
  if (LHS)
    V2 <- lhs::maximinLHS(npts,2)#matrix(runif(npts*2),ncol=2)
  else
    V2 <- matrix(runif(npts*2),ncol=2)
  V2.inds <- sensorpl(V2,kpts,theta)
  mrkr <- rep('o',npts)
  mrkr[V2.inds] <- as.character(1:length(V2.inds))
  #plot(V2,col=((1:npts)%in%V2.inds + 1),pch=mrkr)
  plot(V2[-V2.inds,],xlab='x1',ylab='x2',main='Sensor placement in 2-D')
  points(V2[V2.inds,],pch=as.character(1:length(V2.inds)),col='red')
}
if (F) {
  get.2D.sensors(1,50,8,seed=0,theta=10)
  get.2D.sensors(1,50,8,seed=0,theta=10,LHS=T)
}