
#include <Rcpp.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

// refer to:
// https://github.com/evanbiederstedt/poilogcpp

// Function declarations
double poilog(int x, double my, double sig);
double bipoilog(int x, int y, double my1, double my2, double sig1, double sig2, double ro);
double maxf(int x, double my, double sig);
double upper(int x, double m, double my, double sig);
double lower(int x, double m, double my, double sig);
double my_f(double z, int x, double my, double sig, double fac);
double my_f2(double z,int y,int x,double my1,double my2,double sig1,double sig2,double ro,double fac);
void my_f_vec(double *z, int n, void *p);
void my_f2_vec(double *z, int n, void *p);
NumericVector poilog2(NumericVector x, NumericVector y, double my1, double my2, double sig1, double sig2, double ro, int nrN);
NumericVector poilog1(NumericVector x, double my, double sig, int nrN);

struct My_fparams { 
   int x; 
   double sig; 
   double my; 
   double fac;
};

struct My_f2params { 
   int x;   
   int y; 
   double sig1; 
   double sig2; 
   double my1; 
   double my2;
   double ro; 
   double fac;
};

// [[Rcpp::export]]
NumericVector poilog2(NumericVector x, NumericVector y, double my1, double my2, double sig1, double sig2, double ro, int nrN){
   std::vector<double> vect;
   for (int i = 0; i < nrN; i++){
      vect.push_back(bipoilog(x[i],y[i], my1, my2, sig1, sig2, ro));
      Rcpp::checkUserInterrupt();
   }
   return Rcpp::wrap(vect);
}

// [[Rcpp::export]]
NumericVector poilog1(NumericVector x, NumericVector my, NumericVector sig, int nrN){
   std::vector<double> vect;
   for (int i = 0; i < nrN; i++){
      vect.push_back(poilog(x[i], my[i], sig[i]));
      Rcpp::checkUserInterrupt();
   }
   return Rcpp::wrap(vect);
}


double maxf(int x, double my, double sig){
   double d,z;
   z=0;
   d=100;
   while (d>0.00001) {
      if (x-1-exp(z)-1/sig*(z-my)>0) {
         z=z+d; 
      } else {
         z=z-d;
      }
      d=d/2;
   }
   return(z);
}


double upper(int x, double m, double my, double sig){
   double d,z,mf;
   mf = (x-1)*m-exp(m)-0.5/sig*((m-my)*(m-my));
   z = m+20;
   d = 10;
   while (d>0.000001) {
      if ((x-1)*z-exp(z)-0.5/sig*((z-my)*(z-my))-mf+log(1000000)>0) {
         z=z+d; 
      } else {
         z=z-d;
      }
      d=d/2;
   }
   return(z);
}


double lower(int x, double m, double my, double sig) {
   double d,z,mf;
   mf = (x-1)*m-exp(m)-0.5/sig*((m-my)*(m-my));
   z = m-20;
   d = 10;
   while (d>0.000001) {
      if ((x-1)*z-exp(z)-0.5/sig*((z-my)*(z-my))-mf+log(1000000)>0) {
         z=z-d; 
      } else {
         z=z+d;
      }
      d=d/2;
   }
   return(z);
}


double my_f(double z, int x, double my, double sig, double fac) {
   return exp(z*x-exp(z)-0.5/sig*((z-my)*(z-my))-fac);
}


double my_f2(double z,int y,int x,double my1,double my2,double sig1,double sig2,double ro,double fac){
   return( poilog(y,my2+ro*sqrt(sig2/sig1)*(z-my1),sig2*(1-ro*ro)) *
       exp(x*z-exp(z)-fac-0.5/sig1*(z-my1)*(z-my1)) );
}


void my_f_vec(double *z, int n, void *p){
   struct My_fparams *params = (struct My_fparams *)p;
   int x      = (params->x);
   double sig = (params->sig);
   double my  = (params->my);
   double fac = (params->fac);
   for (int i=0;i<n;i++) {
      z[i]=my_f(z[i],x,my,sig,fac);
   }
   return;
}


void my_f2_vec(double *z, int n, void *p){
   struct My_f2params *params = (struct My_f2params *)p;
   int x       = (params->x);
   int y       = (params->y);
   double sig1 = (params->sig1);
   double sig2 = (params->sig2);
   double my1  = (params->my1);
   double my2  = (params->my2);
   double ro   = (params->ro);
   double fac  = (params->fac);
   for (int i=0;i<n;i++) {
      z[i]=my_f2(z[i],y,x,my1,my2,sig1,sig2,ro,fac);
   }
   return;
}


double poilog(int x, double my, double sig) {
   double a, b, m, fac, val;
   double result, abserr;
   int last, neval, ier;
   int lenw;
   int limit=100;
   double reltol=0.00001;
   double abstol=0.00001;
   lenw = 4 * limit;
   int *iwork = (int *)calloc(limit, sizeof(int));
   double *work = (double *)calloc(lenw, sizeof(double));

   m = maxf(x,my,sig);
   a = lower(x,m,my,sig);
   b = upper(x,m,my,sig);
   fac = lgamma(x+1);

   struct My_fparams p = { x, sig, my, fac };

   Rdqags(my_f_vec, (void *) &p, &a, &b,
      &abstol,&reltol,&result,&abserr,&neval,&ier,
      &limit,&lenw, &last, iwork, work);

   if (ier!=0) {
      Rcpp::stop("error in integration\n");
   } 

   val = result*(1/sqrt(2*M_PI*sig));
   free(iwork);
   free(work);
   return(val);
}


double bipoilog(int x, int y, double my1, double my2, double sig1, double sig2, double ro) {
   double a, b, m, fac, val;
   double result, abserr;
   int last, neval, ier;
   int lenw;
   int limit=100;
   double reltol=0.00001;
   double abstol=0.00001;
   lenw = 4 * limit;
   int *iwork = (int *)calloc(limit, sizeof(int));
   double *work = (double *)calloc(lenw,  sizeof(double));

   m = maxf(x,my1,sig1);
   a = lower(x,m,my1,sig1);
   b = upper(x,m,my1,sig1);
   fac = lgamma(x+1);

   struct My_f2params p = { x, y, sig1, sig2, my1, my2, ro, fac};

   Rdqags(my_f2_vec, (void *) &p, &a, &b,
      &abstol,&reltol,&result,&abserr,&neval,&ier,
      &limit,&lenw, &last,iwork, work);

   if (ier!=0) {
      Rcpp::stop("error in integration\n");
   } 

   val = result*(1/sqrt(2*M_PI*sig1));
   free(iwork);
   free(work);

   return(val);
}

