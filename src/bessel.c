// This file contains the implementation of Bessel functions, according to the Numerical Recipes Book
// [ http://numerical.recipes/book/book.html ]

#include <math.h>

static inline double poly(const double *cof, const int n, const double x) {
    // Evaluate a polynomial
    double ans = cof[n];
    for (int i = n-1; i>=0 ; i--){
        ans = ans*x + cof[i];
    }
    return ans;
}

double BesselK0(double x){
// Returns the modified Bessel function of the second kind and zero-th order K0(x) for positive real x.
// Adapted from the routine `Doub k0` of page 280;

static const double k0pi[5] = { 1.000000000000000e+0 , 2.346487949187396e-1 , 1.187082088663404e-2 , 2.150707366040937e-4 , 1.425433617130587e-6 };
static const double k0qi[3] = { 9.847324170755358e-1 , 1.518396076767770e-2 , 8.362215678646257e-5 };
static const double  k0p[5] = { 1.159315156584126e-1 , 2.770731240515333e-1 , 2.066458134619875e-2 , 4.574734709978264e-4 , 3.454715527986737e-6 };
static const double  k0q[3] = { 9.836249671709183e-1 , 1.627693622304549e-2 , 9.809660603621949e-5 };
static const double k0pp[8] = { 1.253314137315499e+0 , 1.475731032429900e+1 , 6.123767403223466e+1 , 1.121012633939949e+2 , 9.285288485892228e+1 , 3.198289277679660e+1 , 3.595376024148513e+0 , 6.160228690102976e-2 };
static const double k0qq[8] = { 1.000000000000000e+0 , 1.189963006673403e+1 , 5.027773590829784e+1 , 9.496513373427093e+1 , 8.318077493230258e+1 , 3.181399777449301e+1 , 4.443672926432041e+0 , 1.408295601966600e-1 };

    if (x <= 1.0) {
        // Use two rational approximations
        double y = x*x;
        double term = poly(k0pi,4,y)*log(x)/poly(k0qi,2,1.-y);
        return poly(k0p,4,y)/poly(k0q,2,1.-y)-term;
    } else {
        // Rational approximation with exp(-x)/sqrt(x) factored out
        double z = 1.0/x;
        return exp(-x)*poly(k0pp,7,z)/(poly(k0qq,7,z)*sqrt(x));
    }
}