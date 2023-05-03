#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include<gsl/gsl_sf_bessel.h>
#include<gsl/gsl_sf_legendre.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_histogram.h>
#include<gsl/gsl_cdf.h>

using namespace std;

/// Gaussian RNG using GSL
double RNGGaussGSL(gsl_rng * r_gaussian, double sigma, double mean){
  double gauss_val = gsl_ran_gaussian(r_gaussian, sigma)+mean;
  return gauss_val;
}

/// uniform RNG using GSL
double RNGUniformGSL(gsl_rng * r_uniform, double lower_limit, double upper_limit){
  double uniform_val = gsl_ran_flat(r_uniform, lower_limit, upper_limit);
  return uniform_val;
}

// Breit-Wigner Distribution

double RNGBreitWigner(gsl_rng * r,
                      double events,
                      double lower,
                      double upper,
                      double Ezero,
                      double Width,
                      int    angMom,
                      int    fragMass,
                      int    flagGamma){

  double n = 1000*upper; // Number of bins
  int index = n+1;

  gsl_histogram * fH; //holder histogram
  gsl_histogram_pdf * fP; //sampled histogram


  fH = gsl_histogram_alloc(n); // define histo gram and alloc
  gsl_histogram_set_ranges_uniform (fH, lower, upper);
  double y[index] = {0.0}; // Array to hold function

  // Get Constants - copied from st_rng.cc starting at line 330
  double hBarC = 197.3269631; // MeV fm
  double Nmass = 939.565346; //  MeV
  double amu = 931.494028; // from Nist
  double nucRad = 1.4; // fm from Lane and Thomas pg 266 *1.2*
  double Const = (1.0 /3.14159265 ); // Const for BW
  // reduced mass from lane in amu
  double redMass = (((Nmass/amu)*fragMass)/(Nmass/amu + fragMass))*amu;
  // reduced mass for three particles
  // double redMass = (((2.*Nmass/amu)*fragMass)/(2.*Nmass/amu + fragMass))*amu;
  double k = sqrt(2.0 * redMass) / hBarC; // k in terms of sqrt(Energy)
  // Nuclear Radius (A1^(1/3) + A2^(1/3))
  double Radius = nucRad * (pow(fragMass,1.0/3.0) + pow(Nmass/amu, 1.0/3.0));
  double x0 = k*Radius*sqrt(Ezero);

  //Functions with the resonance enrgy inputed
  double jl0 = x0 * gsl_sf_bessel_jl (angMom,x0); // Regular Spherical Bessel * x
  double nl0 = x0 * gsl_sf_bessel_yl (angMom,x0); // Irregular Spherical Bessel * x
  double PenL0 = x0 / (nl0*nl0 + jl0*jl0); // Penetrability as a function of E,W,L

  for (int i=0; i < n; i++) {// loop over bins
    double En=i*(upper-lower)/(n)+(upper-lower)/(2*n); // Convert bins to energy in MeV
    double x = k*Radius*sqrt(En); // X = k*R with energy included

    //Call functions for Penetribility
    double jl = x * gsl_sf_bessel_jl (angMom,x); // Regular Spherical Bessel * x
    double nl = x * gsl_sf_bessel_yl (angMom,x); // Irregular Spherical Bessel * x
    double PenL = x / (nl*nl + jl*jl); // Penetrability as a function of E,W,L

    double bound = 0; // Define the boundary condition
    double shift = 0; // Define the shift function

    //zwk, check this for angMom==0
    if(angMom ==0){
      shift = 0;
      bound = 0;
      PenL = 1;
      PenL0 = 1;
    }

    if(angMom == 1) { //
      shift = - 1 / ( 1 + pow(x,2)); // From lane and thomas, i.e. x*(F'*F+G'*G)/(F*F+G*G)
      bound = -1/ ( 1 + pow(x0,2)); // Set boudary so shift = 0 at eigen value (i.e. Decay energy)
    }

    if(angMom == 2) { // Same as above for l = 2
      shift = -3*(pow(x,2)+6)/(pow(x,4)+9+3*pow(x,2));
      bound = -3*(pow(x0,2)+6)/(pow(x0,4)+9+3*pow(x0,2));
    }

    double redGamma = 0.0;
    double delta = 0.0;
    double Gamma = 0.0;
    double obsGamma = 0.0;

    if(flagGamma == 1) {  // For use with Gflag - reduced width
      redGamma = Width;
      delta = -( shift - bound ) * redGamma * redGamma; // Full shift function
      Gamma = 2.0 * redGamma * redGamma * PenL; // Width
    } else if (flagGamma == 0){ // default --observed width CRH 5/12/08
      obsGamma = Width;
      delta = - ( shift - bound)* (obsGamma / (2 * PenL0));
      Gamma = (obsGamma / PenL0)*PenL;
    }

    // Full BW with shift and pen functions included pg 322 Lane and Thomas
    double BW = Const * ( (Gamma / 2.0) / ((Ezero + delta - En)*(Ezero + delta - En) + (Gamma*Gamma) / 4.0));

    y[i]=BW; // Add value to array
    gsl_histogram_accumulate (fH,En,y[i]); // add value to histogram

  }// loop over bins

  // Create PDF from histogram
  fP = gsl_histogram_pdf_alloc(n); // init PDF
  gsl_histogram_pdf_init(fP,fH); // Make pdf with histogram dist.

  //const gsl_rng_type * T;
  //gsl_rng * rng;
  //T = gsl_rng_default;
  //rng = gsl_rng_alloc(T);

  double uniform_choice = gsl_ran_flat(r, 0, 1);
  double val = gsl_histogram_pdf_sample (fP,uniform_choice);

  return val;
}
