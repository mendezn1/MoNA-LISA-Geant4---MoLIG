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

#include <vector>
#include <Math/GenVector/Boost.h>
#include <TLorentzVector.h>
#include <typeinfo>
#include "random_number_generators.hh"

using namespace std;

vector<vector<vector<double>>> kinematics(int events, double E_o, double width, int l, int mass_fragment) {

  vector<vector<double>> neutron_lab_4vec;
  vector<vector<double>> fragment_lab_4vec;
  double amu = 931.49410242 ;// atomic mass in MeV/c^2 (NIST)
  double m_neutron = 1.00866491590 * amu  ;/// neutron amu mass from in MeV/c^2 (AME2020) (0.000000000047 - sigma)
  double m_electron = 0.000548579909065 * amu; // electron amu mass from in MeV/c^2 (NIST)
  double m_fragment = 51.963213646 * amu - 20*m_electron ;/// Ca-52 fragment mass in MeV/c^2 (AME2020)(0.000000720 - sigma)
  double m_parent = 52.968451000 * amu - 20*m_electron;/// Ca-53 fragment mass in MeV/c^2 (AME2020) (0.000047000 - sigma)
  double m_beam =  53.963029359 * amu - 21*m_electron;/// Sc-54 secondary beam mass in MeV/c^2 (AME2020) (0.000015000 - sigma)
  double m_ex_neutron = 8071.31806/1000  ;/// neutron mass excess in MeV (AME2020) (0.00044 - sigma)
  double m_ex_fragment = -34266.272/1000  ;/// Ca-52 mass excess in MeV (AME2020) (0.671 - sigma)
  double m_ex_parent = -29387.707/1000 ;/// Ca-53 mass excess in MeV (AME2020) (43.780 - sigma)
  double m_ex_beam =  -34437.934/1000 ;/// Sc-54 mass excess in MeV (AME2020) (13.973 - sigma)
  double s_n_parent = m_ex_fragment+m_ex_neutron-m_ex_parent ;/// MeV | Neutron separation energy of Ca-53 (3.195 MeV (Calculated))
  double ke_beam = 180*52.968451000 ; // Kinetic energy of the secondary Sc-54 beam
  double gamma_parent = 1+ke_beam/m_parent ; // calculated gamma values
  double beta_parent = pow(1-(pow(gamma_parent,-2.0)),0.5); // calculated beta value

  const gsl_rng_type * T;
  T = gsl_rng_default;

  gsl_rng * r_bw;
  gsl_rng * r_uniform;
  r_bw      = gsl_rng_alloc(T);
  r_uniform = gsl_rng_alloc(T);

  TLorentzVector neutron_4vec;
  TLorentzVector fragment_4vec;
  double beta_x = 0.0;
  double beta_y = 0.0;
  double beta_z = beta_parent;

  for (int ijk = 0; ijk < events; ijk++){
    // Holder vector for storing 4-vectors
    vector<double> p_holder;
    vector<double> f_holder;

    // Get Briet-Wigner Distribution
    double bw_energy = RNGBreitWigner(r_bw, 1, 0, 10, E_o, width, l, mass_fragment, 0);
    // Total Energy of Unbound Nucleus
    double e_parent = m_parent+s_n_parent+bw_energy;
    // Solve for CoM Momentum
    double p_neutron = pow(pow((pow(e_parent,2)+pow(m_neutron,2)-pow(m_fragment,2))/(2*e_parent),2)-pow(m_neutron,2), 0.5);
    // Random Theta and Phi angles
    double theta = RNGUniformGSL(r_uniform, 0., 180.);
    double phi = RNGUniformGSL(r_uniform, 0., 360.);
    // Create neutron and fragment 4-vectors
    neutron_4vec.SetPxPyPzE(p_neutron*sin(theta)*cos(phi), p_neutron*sin(theta)*sin(phi), p_neutron*cos(theta), pow(pow(m_neutron,2)+pow(p_neutron,2),0.5)); // <px,py,pz,E_tot>
    fragment_4vec.SetPxPyPzE(-p_neutron*sin(theta)*cos(phi), -p_neutron*sin(theta)*sin(phi), -p_neutron*cos(theta), pow(pow(m_fragment,2)+pow(p_neutron,2),0.5)); // <px,py,pz,E_tot>
    // Boost Fragments into Lab Frame
    neutron_4vec.Boost(beta_x,beta_y,beta_z);
    fragment_4vec.Boost(beta_x,beta_y,beta_z);
    // Push Lab 4-vector components into holder vectors
    p_holder.push_back(neutron_4vec[0]);
    p_holder.push_back(neutron_4vec[1]);
    p_holder.push_back(neutron_4vec[2]);
    p_holder.push_back(neutron_4vec[3]);
    f_holder.push_back(fragment_4vec[0]);
    f_holder.push_back(fragment_4vec[1]);
    f_holder.push_back(fragment_4vec[2]);
    f_holder.push_back(fragment_4vec[3]);
    // Push neutron and fragment 4-vectors
    neutron_lab_4vec.push_back(p_holder);
    fragment_lab_4vec.push_back(f_holder);
  }
  vector<vector<vector<double>>> data;
  data.push_back(neutron_lab_4vec);
  data.push_back(fragment_lab_4vec);
  return data;
}
