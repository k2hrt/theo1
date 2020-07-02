logtics <- function(x, y)
{
  # Get log scale limits
  left = floor(log10(min(x)))  
  right = ceiling(log10(max(x)))
  bot = floor(log10(min(y)))  
  top = ceiling(log10(max(y)))
  # Draw X-Axis Tics
  for(xx in left:right)
  { 
    axis(side=1,at=(2:9)*10^xx,labels=FALSE)
  }
  # Draw Y-Axis Tics
  for(yy in bot:top)
  { 
    axis(side=2,at=(2:9)*10^yy,labels=FALSE)
  }
} 
