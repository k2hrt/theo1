/*****************************************************************************/
/*                                   TheoBR.cpp                              */
/*                  All Tau Bias Removed Theo1 Calculation for R             */
/*            Using Fast Lewis Theo1 and Taylor/Howe Fast TheoBR Methods     */
/*                                                                           */
/* This code calculates a bias-removed all tau version of Theo1 for R,       */
/* based on the Lewis fast Theo1 algorithm and the Taylor/Howe fast Theo1    */
/* method.  It is completely contained within this file.  The N point phase  */
/* data x[], and the (N-1)/2 point Theo1 t[], upper [u] and lower [l] error  */
/* bar, and Theo1 tau[] array arguments are NumericVectors.  The data        */
/* sampling interval tau0 is input as a double with a default value of 1,    */
/* the confidence factor as a double with default value 0.683, and the bias  */
/* removal flag as an integer (0 or 1) with defalut value 1.  The  x[] is    */
/* unchanged, t[] receives the Theo1BR results, u[] and l[] receive the      */
/* upper and lower Theo1 error bar values, tau[] receives the corresponding  */
/* Theo1 tau values, and the function returns a double, either the Theo1     */
/* bias factor, 0 if no correction, or an error code.  The tau values are    */
/* 0.75 * tau0 * even averaging factors.  If the last br argument is 0,      */
/* the function returns 1.0 with the raw uncorrected Theo1 values.  Chi-     */
/* squared error bars are calculated for each AF based on the # of data      */
/* points, the Theo1 EDF, and the desired CF.  A negative CF entry denotes   */
/* use of single-sided confidence limits.                                    */
/*                                                                           */ 
/* To set up:                                                                */
/*   Install Rcpp with: > library(Rcpp)                                      */
/*   Compile with: > sourceCpp(“Theo1BR.cpp”)                                */
/*   after adding full path to filename e.g., "C:\\R\\Theo1BR.cpp"           */
/*                                                                           */
/* To use:                                                                   */
/*   Load phase data into vector x (the tau0 must be known)                  */
/*   e.g., > x<-scan("C:\\Data\\phase.dat")                                  */
/*   Determine the # of phase data points with > N=length(x)                 */ 
/*   Allocate results vector t with: > t<-numeric((N-1)/2)                   */
/*   Allocate u vector tau with: > u<-numeric((N-1)/2)                       */
/*   Allocate l vector tau with: > l<-numeric((N-1)/2)                       */
/*   Allocate tau vector tau with: > tau<-numeric((N-1)/2)                   */
/*                                                                           */  
/* Call with:                                                                */
/* > Theo1BR(x, t, u, l, tau, tau0, cf, br)                                  */
/* where tau0, cf and br have defaults of 1.0, 0.683 and 1 respectively      */
/*                                                                           */
/* Display results with:                                                     */
/*   > t or create a data frame                                              */
/*   > r<-data.frame(t,u,l,tau) and display it with > r                      */
/* or plot the results with (for example, red line, no error bars):          */
/*   > plot(tau,t,type="l",main="Theo1BR Plot",xlab="Tau",ylab="Theo1BR",    */
/*   + log="xy", col="red")                                                  */
/* To add lines representing upper and lower error bounds:                   */
/*   > points(tau, u, type="l", col="blue")                                  */
/*   > points(tau, l, type="l", col="blue")                                  */
/*                                                                           */
/* To show the error bounds as a grey area on the plot:                      */
/*   > plot(tau,t,type="l",main="Theo1BR Plot",xlab="Tau",ylab="Theo1BR",    */
/*   + log="xy", col="red")                                                  */ 
/*   > polygon(c(tau, rev(tau)), c(u, rev(l)), col = "grey70", border = NA)  */
/*   > points(tau,t,type="l",col="red")                                      */
/*                                                                           */
/* Other R plot formats can be used as desired                               */
/* An R wrapper function can encapsulate the complete analysis process       */
/*                                                                           */
/* References:                                                               */
/*  [1] B. Lewis, "W Fast Algorithm for Calculation of Theo1", Submitted     */
/*      to IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency   */
/*      Control, May 2020.                                                   */
/*  [2] J.A. Taylor and D.A. Howe, "Fast Theo1BR: A Method for Long Data     */
/*      Set Stability Analysis", IEEE Transactions on Ultrasonics, Ferro-    */
/*      electrics, and Frequency Control, September 2010. I                  */
/*                                                                           */
/* Revision Record:                                                          */
/*  06/10/20  Created and running                                            */
/*  06/16/20  Added draft error bar code                                     */
/*  06/17/20  Editing and debugging                                          */
/*  06/18/20  Changed name to TheoBR()                                       */
/*  06/19/20  Added provisions for single-sided confidence limits            */
/*  06/23/20  Validated and documented                                       */
/*                                                                           *//* W.J. Riley, Hamilton Technical Services, Beaufort, SC 29907 USA © 2020    */
/* E-mail: bill@wriley.com                                   License: MIT    */
/*                                                                           */
/*****************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

/*****************************************************************************/
/* Calculate TheoBR for a phase dataset taken at time increments tau0.  The  */
/* phase data are copied to a working array X, and any linear component to   */
/* the dataset is removed as a first step.  The raw all tau Theo1 values are */
/* are calculated using the Lewis fast algorithm of Reference [1], the       */
/* corresponding tau values are written, and Theo1BR bias corrections are    */
/* made per Reference [2].  Chi-squared upper and lower error bars are       */
/* returned in the u and l arrays for the given confidence factor. The # of  */
/* Theo1 results is (N-1)/2.  This code should preferably be compiled with:  */
/* -O3 -ffast-math.                                                          */
/*                                                                           */
/* Arguments:                                                                */
/*   x    = NumericVector = Phase data input                                 */
/*   t    = NumericVector = Theo1BR output                                   */
/*   u    = NumericVector = Upper Theo1BR error bar                          */
/*   l    = NumericVector = Lower Theo1BR error bar                          */
/*   tau  = NumericVector = Tau output                                       */
/*   tau0 = double = Data sampling time input (default=1)                    */
/*   cf   = double = Confidence factor                                       */
/*                   Negative = Single-Sided (default=0.683)                 */
/*   br   = integer = Bias remove flag (0=no, 1=yes, default=1)              */
/*                                                                           */
/* Return = double = Theo1 Bias Factor                                       */
/*          or 1.0 is no bias correction                                     */
/*          or 0 if error                                                    */
/*****************************************************************************/

// [[Rcpp::export]]

// The core Theo1BR computational routine was created by
// Ben Lewis of the University of Strathclyde UK
// b.lewis@strayh.ac.uk per Reference [1] above © 2020 IEEE

// TheoBR function 
double TheoBR(NumericVector x, NumericVector t, NumericVector u,
  NumericVector l, NumericVector tau, double tau0=1.0, double cf=0.683,
  int br=1){

  // Find phase array size
  int N=x.size();

  // Copy phase data to working array
  double* X = new double[N];
  for(int i = 0; i < N; i++){
    X[i] = x[i];
  }
  
  // Initializations
  int k_max = (N-1)/2;
  int single = 0; // Flag for single-sided confidence limits
  if(cf<0){
	 single=1;
	 cf = -cf;
  }
  double* C1 = new double[N];
  double* C3 = new double[k_max+1];
  double* C4 = new double[k_max*2];

  // Preprocess by removing linear part
  double midpoint = ((double) (N-1))/2;
  long double sum1 = 0;
  long double sum2 = 0;
  for(int i = 0; i < N; i++){
    sum1 += X[i];
    sum2 += X[i]*(i-midpoint);
  }
  double a = sum1/N;
  double b = sum2/N*12/((double)N*N-1);
  for(int i = 0; i < N; i++){
    X[i] -= a + b*(i-midpoint);
  }

  // Calculate C1
  double s=0;
  for(int i = 0; i < N; i++){
    s += (X[i]*X[i]);
    C1[i] = s;
  }

  // Main loop
  C3[0] = C1[N-1];
  for(int k=1; k<=k_max; k++){
    //Calculate C2 values
    double C2_2k = 0;
    double C2_2k_1 = 0;
    for(int j = 0; j <= N-2*k-1; j++){
      C2_2k += (X[j]*X[j+2*k]);
      C2_2k_1 += (X[j]*X[j+2*k-1]);
    }
    C2_2k_1 += (X[N-2*k]*X[N-1]);

    // Update C3, C4 in place
    for(int v=0; v < k; v++){
      C3[v] -= (X[k-1-v]*X[k-1+v])
             + (X[N-k+v]*X[N-k-v]);
    }
    for(int v = 1;v<=2*k-2;v++){
      C4[v-1] -= (X[2*k-1-v]*X[2*k-1])
               + (X[2*k-2-v]*X[2*k-2])
               + (X[N-2*k]*X[N-2*k+v])
               + (X[N-2*k+1]*X[N-2*k+1+v]);   
    }
    C3[k] = C2_2k;
    C4[2*k-2] = 2*C2_2k_1 - (X[0]*X[2*k-1])
              - (X[N-2*k]*X[N-1]);
    C4[2*k-1] = 2*C2_2k;

    //Calculate un-normalised T_k from C1-C4
    double T_k = 0;
    double A0 = C1[N-1] - C1[2*k-1]
              + C1[N-2*k-1] + 2*C2_2k;
    for(int v = 1;v<=k;v++){
      double A1 = A0 - C1[v-1]
                + C1[N-1-v] - C1[2*k-v-1]
                + C1[N-1-2*k+v];
      double A2 = C3[k-v] - C4[v-1]
                - C4[2*k-v-1];
      T_k += (A1+2*A2)/v;
    }

    // For testing - check for negative Theo1
    // if(T_k<0) return 0;
    
    // Apply normalisation to get Theo1
    // t[k-1] = (T_k/(3*(double)(N-2*k)*k*k));
    t[k-1] = (T_k/(3*(double)(N-2*k)*k*k*tau0*tau0));
    if(br==0){
      // Take square root for Theo1 deviation
      t[k-1] = sqrt(t[k-1]);
    }

    // Enter tau = averaging factor * tau0
    tau[k-1] = 0.75*2*k*tau0;
  }

  //Release memory
  delete[] C1;
  delete[] C3;
  delete[] C4;
  
  // Return raw Theo1 if no bias correction
  if(br==0){
    return 1.0;
  }
  
  // Now have raw all tau Theo1 variances and their corresponding tau values.
  // Next need to determine the bias correction factor kf using n=(0.1*N/3)-3
  // AVAR/Theo1 variance ratio averages per Ref [2] (other n's have been used,
  // e.g., (N/6)-3).  The AVAR and Theo1 variance values used are for the same
  // taus, with AFs of 9+3i and 12+4i respectively where i goes from 0 to n
  // (sometimes n to 0 has been used instead).  The Theo1 variance values come
  // from the t[] array with the index above.
  double kf = 0.0;
  int n = ((0.1*N)/3)-3;
  int m; // Averaging factor
  
  // Correction factor summation loop
  for(int i=0; i<=n; i++){
    // Calculate AVAR at AF=m=9+3i
    m=9+3*i;
    double avar = 0.0;
    for(int j=0; j<N-2*m; j++){
      // avar += pow((X[j+2*m] - 2*X[j+m] + X[j]), 2);
      avar += ((X[j+2*m] - 2*X[j+m] + X[j])*(X[j+2*m] - 2*X[j+m] + X[j]));
    }
    
    // Apply scale factor
    double M=m;
    avar /= (2*M*M*(N-2*M)*tau0*tau0);
    
    // For testing
    // Rprintf("i=%d, m=%d, adev=%e, theo1 dev=%e\n", i,m,sqrt(avar),sqrt(t[5+2*i]));
    
    // Divide AVAR by corresponding Theo1 variance
    // and add that to correction factor sum
    kf += (avar / t[5+2*i]);
  }
  // Divide kf sum by # ratios
  kf /= (n+1);
  
  // For testing
  // Rprintf("kf=%e\n", kf);
  
  // For testing - Force kf=1.0
  // kf=1.0;
  
  // Apply bias corrections and convert to Theo1 deviation.
   for(int i=0; i < k_max; i++){
    t[i] *= kf;
    t[i] = sqrt(t[i]);
   }
  
  // To set error bars, must determine
  // dominant power law noise type from
  // Theo1 bias correction factor, then
  // calculate Theo1 edf, and use that to
  // determine upper and lower chi-squared
  // confidence limits at (say) 68% 
  // 1-sigma Confidence Factor

  // Use Theo1 bias factor to determine noise type
  // There is one bias factor value and one noise
  // type determination for the whole run 
  // Thresholds for noise sorting are geometric means
  // of nominal Theo1 variance bias values

  // Note: Nominal Theo1 variance bias factors in
  // original 2003 FCS paper were:
  // WPM=0.4, FPM=0.6, WFM=1.0, FFM=1.71, RWFM=2.24
  // These are not used here
  
  // Error bar determination loop
  for(int i=0; i < k_max; i++){
    // Get nominal Theo1 var bias factors for various
    // power law noise types versus averaging factor
    m=i+1;
    // double mm = 0.75*m;
    double mm = 0.75*2*m;
    double bWPM = (0.09 + (0.74/pow(mm, 0.40)));
    double bFPM = (0.14 + (0.82/pow(mm, 0.30)));
    double bWFM = 1.00;
    double bFFM = (1.87 + (-1.05/pow(mm, 0.79)));
    double bRWFM = (2.70 + (-1.53/pow(mm, 0.85)));
    
    // For testing
    // Display the nominal bias factors
    // Rprintf("mm=%d, WPM=%f, FPM=%f, WFM=%f, FFM=%f, RWFM=%f\n",
    //   mm,bWPM,bFPM,bWFM,bFFM,bRWFM);

    // Find noise type
    int a; // Power law noise exponent alpha
    if(kf<sqrt(bWPM*bFPM)) a=2;
    else if(kf<sqrt(bFPM*bWFM)) a=1;
    else if(kf<sqrt(bWFM*bFFM)) a=0;
    else if(kf<sqrt(bFFM*bRWFM))a=-1;
    else a=-2;

    // Find edf
    // N is # phase data points
    // mm is tau-s the Theo1 stride = 0.75*m
    double edf; // Theo1 EDF
    if(a==2) edf=(0.86*(N+1)*(N-((4.0/3.0)*mm)))/
      ((N-mm)*(mm/(mm+1.14)));
    if(a==1) edf=((4.798*N*N)-(6.374*N*mm)+(12.387*mm))/
      (sqrt(mm+36.6)*(N-mm)*(mm/(mm+0.3)));
    if(a==0) edf=(((4.1*N+0.8)/mm)-((3.1*N+6.5)/N))*
      (pow(mm,1.5)/(pow(mm,1.5)+5.2));
    if(a==-1) edf=(2.0*N*N-1.3*N*mm-3.5*mm)/
      (N*mm)*((mm*mm*mm)/(mm*mm*mm+2.3));
    if(a==-2) edf=(4.4*N-2)/(2.9*mm)*(((4.4*N-1)*
      (4.4*N-1)-8.6*mm*(4.4*N-1)+
      11.4*mm*mm)/((4.4*N-3)*(4.4*N-3)));
 
    edf=((((4.1*N)+0.8)/mm)-(((3.1*N)+6.5)/N))*(pow(mm,1.5)/(pow(mm,1.5)+5.2));    
        
    // Find inverse upper and lower Chi-squared values for given confidence
    // factor and edf with an embedded R calls.
    // See: https://teuder.github.io/rcpp4everyone_en/
    // 220_dpqr_functions.html#chi-squared-distribution
    // Set confidence interval and fill in upper and lower error bar arrays
    // Lower error bar equals nominal TheoBR for single-sided conf interval
    if(single){
	  double upper = R::qchisq((cf),edf,0,0);
      u[i] = sqrt(t[i]*t[i]*edf/upper);
      l[i] = t[i]; 
    }
    else{ // Double-sided confidence intervals:
      double lower = R::qchisq((1-cf)/2,edf,0,0);
      double upper = R::qchisq(1-(1-cf)/2,edf,0,0);
      u[i] = sqrt(t[i]*t[i]*edf/upper);
      l[i] = sqrt(t[i]*t[i]*edf/lower);
    }
    
    // For testing
    // if(i==4) Rprintf("tau=%f, min=%f, nom=%f, max=%f\n", tau[i],l[i],t[i],u[i]);
  }

  // Done
  return kf;
}
