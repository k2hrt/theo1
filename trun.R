trun <- function(x, tau0=1.0, ci=0.683) 
{
  # Allocate arrays
  N=length(x)
  t<-1:((N-1)/2)-1
  u<-1:((N-1)/2)-1 
  l<-1:((N-1)/2)-1
  tau<-1:((N-1)/2)-1

  # Call TheoBR()
  cf=TheoBR(x,t,u,l,tau,tau0,ci,1)

  # Print the bias correction factor
  print(cf)

  # Crudely estimate noise type using original nominal cf
  # WPM=0.4, FPM=0.6, WFM=1.0, FFM=1.71, RWFM=2.24
  if(cf<sqrt(0.4*0.6)) a=2
  else if(cf<sqrt(0.6*1.0)) a=1
  else if(cf<sqrt(1.0*1.71)) a=0
  else if(cf<sqrt(1.71*2.24))a=-1
  else a=-2
  
  # Print power law noise exponent
  print(a)

  # Create data frame for results
  r<-data.frame(tau,l,t,u)

  # Get plot y scale
  top=10^(ceiling(log10(max(t))))
  bot=10^(floor(log10(min(t))))
  
  # Plot results
  plot(tau,t,type="l",ylim=c(bot,top),main="Theo1BR Plot",xlab="Tau", ylab="Theo1BR",log="xy",col="red") 
  polygon(c(tau, rev(tau)), c(u, rev(l)), col = "grey70", border = NA)
  points(tau,t,type="l",col="red")
  grid()

  # Return data frame without printing it
  invisible(r)
} 
