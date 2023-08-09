cutoff<-function(dat=list(),threshold){
  for (i in 1:length(dat)) {
    if(as.numeric(dat[i])<as.numeric(dat[1])*threshold){
      return(i-1)
    }
    if(i==length(dat))
    {return (length(dat))}
  }
  }


