#include "yaml-cpp/yaml.h"
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>

#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"

// Physical parameters
// (they can be read from the LHAPDF .info file corresponding ot the FF set)
const double mc         = 1.51;     //GeV
const double mb         = 4.92;     //GeV
const double mt         = 175.;     //GeV
const double mtau       = 1.777;    //GeV
const double MZ         = 91.1876;  //GeV
const double alphasMZ   = 0.118;
const double alphatau   = 1./128.;
const double sin2thetaw = 0.23122;
const double MPI        = 0.139570; //GeV
const double MPI2       = MPI * MPI;

int main()
{
 // Perturbative order
 const int pto = 1; //NLO

 // FF set
 const string FFset = "MAPFF10NLOPIsum";

 // Read BABAR data
 YAML::Node table = YAML::LoadFile("Table1.yaml");
 const double sqrt_s = table["dependent_variables"][0]["qualifiers"][1]["value"].as<double>();
 const auto xs = table["independent_variables"][0]["values"];
 
 //Initialise APFEL
 APFEL::SetAlphaQCDRef(alphasMZ, MZ);       // Alpha QED at MZ
 APFEL::SetAlphaQEDRef(alphatau, mtau);     // Alpha QCD at mtau
 APFEL::SetSin2ThetaW(sin2thetaw);          // Weinberg angle
 APFEL::SetPoleMasses(mc, mb, mt);          // Heavy quark masses
 APFEL::SetTauMass(mtau);                   // Tau mass
 APFEL::SetZMass(MZ);                       // Z-boson mass
 APFEL::SetTimeLikeEvolution(true);         // Enable time-like evolution
 APFEL::SetAlphaEvolution("exact");         // Exact solution of the RGE
 APFEL::SetPerturbativeOrder(pto);          // Set perturbative order
 APFEL::SetProcessDIS("NC");                // Set NC process
 APFEL::SetPDFSet(FFset);                   // Set PDF set
 APFEL::SetReplica(0);                      // Set FF member (here replica 0)
 APFEL::SetNumberOfGrids(3);                // Set interpolation grids used to evaluate structure functions
 APFEL::SetGridParameters(1,100,3,1e-2);    // Grid 1: 100 points in x in [1e-2,1] with a cubic spline
 APFEL::SetGridParameters(2,60,3,5e-1);     // Grid 2: 60  points in x in [5e-1,1] with a cubic spline
 APFEL::SetGridParameters(3,50,3,8e-1);     // Grid 3: 50  points in x in [8e-1,1] with a cubic spline
 APFEL::InitializeAPFEL_DIS();              // Initialise APFEL
 APFEL::ComputeStructureFunctionsAPFEL(sqrt_s, sqrt_s);
 // N.B. The initial and final scales are equal; no evolution needed.

 std::cout << std::setw(12) << "z" << "  "
	   << std::setw(12) << "1/sig dsig/dp"
	   << std::endl;
 
 for (int i = 0; i < xs.size(); i++)
   {
     double s = sqrt_s * sqrt_s;
     double p = (xs[i]["high"].as<double>() + xs[i]["low"].as<double>()) / 2.0;
     double z = 2.0 * sqrt(MPI2 + (p*p)) / sqrt_s;
     double dpdz  = z * s / (4.0 * p);
     
     std::cout << std::setprecision(5)
	       << std::setw(12)
	       << z << "  "
	       << std::setprecision(4)
	       << std::setw(12)
	       << APFEL::GetSIATotalCrossSection(0,sqrt_s,"total") * ( APFEL::F2light(z) + APFEL::F2charm(z) ) / ( APFEL::GetSIATotalCrossSection(pto,sqrt_s,"light")+ APFEL::GetSIATotalCrossSection(pto,sqrt_s,"charm") ) / dpdz
       	       << std::setprecision(4)
	       << std::setw(12)
	       << std::endl;
   }

 return 0;
 
}
