#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NO_BINS 400
#define k       10

double weight(double x) {return (exp(-pow(x,2)))/sqrt(M_PI);}
double func(double x) {return pow(x,2);}

int main()
{
  double xn,xt,del,wSum,aucor,rate;
  int    bins[NO_BINS],bin_n,i,n,naccept;
  double norm,mean,sigma,integrate,integrate2,fki,RND,R,eta,normalize_value;
  FILE *fp;

  /* get number steps to take, from the keyboard */
  printf("Enter number of steps: ");
  scanf("%i",&n);
  double acorrelation[n];
  /* initialization */
  xn  = 0;                           // starting point
  del = 3;                             // 3.0 step size
  norm = 0.;
  mean=0.0;
  normalize_value=sqrt(M_PI);        //integration w(x) over [-inf,inf]=normalize_value//
  naccept=0;
  sigma=0.0;                     
  fki=0.0;
  integrate=0.0;
  integrate2=0.0;
  for(i=0; i<NO_BINS; i++) bins[i] = 0;

  /* random walker */
  for(i=0; i<n; i++) {

    /* metropolis algorithm consisting of:
     * 1- trial step 
     * 2- move or stay 
     */
     RND=( (double)(rand())/(double)(RAND_MAX));
     xt=xn+del*(2*RND-1);
     R=weight(xt)/weight(xn);
     if(R>1)
     {
     xn=xt;
     naccept++;
     }
     else
     {eta=( (double)(rand())/(double)(RAND_MAX));
     if(R>eta)
     {
     xn=xt;
     naccept++;
     }
     }
     

    /* integral and weighs calculation */
    wSum       += weight(xn);
    integrate  += func(xn);
    integrate2 += pow(func(xn),2);
    /* auto correlation */
    
    acorrelation[i]=xn;
    /* make histogram */
    //    printf("%f\n",xn);
    bin_n = (int) floor(xn*50) + 200;
    if(bin_n < NO_BINS && bin_n>=0)
      bins[bin_n]++;
  } 

  /* Integral, sigma and auto-correlation */
  mean=integrate/n;
  
  norm=integrate2/n;
  sigma=sqrt((((norm)-pow(mean,2)))/(n-1));
  rate=((double)(naccept)/(double)(n));
  for(i=0;i<n-k;i++)
  fki+=(func(acorrelation[i])*func(acorrelation[i+k]))/(n-k);  
  aucor=(fki-pow(mean,2))/(norm-pow(mean,2));
  /* print out acceptance ratio, auto-correlation, Integral and sigma*/
  printf("integration f(x)=x^2*exp(-x^2) over [-inf,inf]=%f\n",mean*normalize_value);
  printf("sigma=%1f\n",sigma);
  printf("auto-correlation=%f\n",aucor);
  printf("acceptAnce rate=%f\n",rate);

  /* write histogram into a ascii file */
  fp = fopen("datafie.dat","w");
  for(i=0; i<NO_BINS; i++) fprintf(fp,"%i\t%i\n",i-200,bins[i]);
  fclose(fp);

  return 0; 
} 

