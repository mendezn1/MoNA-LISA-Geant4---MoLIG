#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <numeric>

using namespace std;

vector<double> edecay_v2(int events, vector<double> neutron_TotalEnergy, vector<double> neutron_px, vector<double> neutron_py, vector<double> neutron_pz,
                                     vector<double> fragment_TotalEnergy, vector<double> fragment_px, vector<double> fragment_py, vector<double> fragment_pz){
    vector<double> edecay_v2;

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

    for (int i=0; i<events; i++){
        double edecay_val_v2 = pow(pow(m_fragment,2)+pow(m_neutron,2)+2*(neutron_TotalEnergy[i]*fragment_TotalEnergy[i]-(neutron_px[i]*fragment_px[i]+neutron_py[i]*fragment_py[i]+neutron_pz[i]*fragment_pz[i])),0.5)-m_fragment-m_neutron;
        edecay_v2.push_back(edecay_val_v2);
    }
    return edecay_v2;
}