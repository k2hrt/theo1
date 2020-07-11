sread<-function(x,r)
{
  # Get size of phase data array
  N=length(x)
  
  # Set size of working arrays
  # n=((N-1)/2)-1
  n=(N-1)/2
  
  # Allocate working arrays
  tau<-1:n
  min<-1:n
  nom<-1:n
  max<-1:n
  num<-1:n
 
  # Fill working arrays from data frame
  tau<-r[,1]
  min<-r[,2]
  nom<-r[,3]
  max<-r[,4]
  
  # Fill # analysis points array
  for(i in 1:n){
    af=(tau[i]/0.75)/tau0    
    num[i]=((N-af)/2)*af
  }
  
  # Create output data frame
  # This is Stable32 Theo1 Read format
  # DF column not included
  rr<-data.frame(tau, num, nom, min, max)
  
  # Write output to file
  # Note: Filename is hard-coded
  write.table(rr, "C:\\Data\\sigma.tau", row.names=F, col.names=F)
  
  # Return output data frame
  invisible(rr)
}
