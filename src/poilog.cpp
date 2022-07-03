
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
void my_f_vec(double *z, int n, void *p);
std::vector<double> poilog1(std::vector<int> x, std::vector<double> my, std::vector<double> sig);

struct My_fparams { 
   int x; 
   double sig; 
   double my; 
   double fac;
};

// [[Rcpp::export]]
std::vector<double> poilog1(std::vector<int> x, std::vector<double> my, std::vector<double> sig){
   int nrN = x.size();
   std::vector<double> vect(nrN);
   for (int i = 0; i < nrN; i++) {
      vect[i] = poilog(x[i], my[i], sig[i]);
   }
   return vect;
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

// [[Rcpp::export]]
double l_lnpois_cpp(std::vector<int> Y_obs, std::vector<double> lambda_ref, int d, double mu, double sig, double phi = 1.0) {

   int n = Y_obs.size();

   double l = 0;
   double p = 0;

   for (int i = 0; i < n; i++) {
      p = poilog(Y_obs[i], mu + log(phi * d * lambda_ref[i]), std::pow(sig, 2));
      if (p == 0) {p = 1e-15;}
      l += log(p);
   }

   return l;
}