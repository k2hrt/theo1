sread2<-function(x,r,tau0=1)
{
  # Get size of phase data array
  N=length(x)
  
  # Get # data frame rows
  n=nrow(r)
  
  # Allocate # analysis points array
  num<-1:n
  
  # Fill # analysis points array
  # ith row of tau column is r[i,1]
  for(i in 1:n){
    af=(r[i,1]/0.75)/tau0    
    num[i]=((N-af)/2)*af
  }
  
  # Create output data frame
  # This is Stable32 Theo1 Read format
  # DF column not included
  # Add num column to data frame
  r<-cbind(r,num)
  # Re-arrange data frame columns
  # Now: tau,l,t,u,num
  # Want: tau,num,t,l,u
  r<-r[c(1,5,3,2,4)]
  
  # Write output to file
  # Note: Filename is hard-coded
  write.table(r, "C:\\Data\\sigma.tau", row.names=F, col.names=F)
  
  # Return output data frame
  invisible(r)
}
