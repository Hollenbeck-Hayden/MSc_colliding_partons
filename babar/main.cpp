#include "yaml-cpp/yaml.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"

#include "../common/common.h"
using namespace common;

// Equations are take from:
// Determination of FFs of Pions, Kaons, Protons - Bertone, Carrazza, Hartland, Nocera, Rojo (2017) 

// z kinematic cut bounds
const double Z_MAX = 0.9;
const double Z_MIN = 0.075;

double theory_xsec(double z, double sqrt_s)
{
	const double zs2 = z * sqrt_s / 2.0;
	const double p = sqrt(zs2 * zs2 - M_PION * M_PION);
	const double dpdz = z * sqrt_s * sqrt_s / (4.0 * p);

	return F2<CHARM>(z) * SIATotal<TOTAL, 0>(sqrt_s) / SIATotal<CHARM>(sqrt_s) / dpdz;
}


int main( void )
{
	const std::string& pdf_name = "MAPFF10NLOPIsum";

	APFEL_SetSIA();
	APFEL_init(pdf_name);

	// Initialize vectors
	std::vector<double> zs, z_lows, z_highs;
	std::vector<double> experimental;
	std::vector<double> uncorr_uncertainty, norm_uncertainty;

	// Parse BaBar data
	double sqrt_s = parse_datafile("Table1.yaml", zs, experimental,
		[] (YAML::Node root) {
			return root["dependent_variables"][0]["qualifiers"][1]["value"].as<double>();
		},
		[&z_lows, &z_highs, &uncorr_uncertainty, &norm_uncertainty]
		(YAML::Node x, YAML::Node y, double q, double& z, double& ex) {
			const double p_low = x["low"].as<double>();
			const double p_high = x["high"].as<double>();

			const double z_low  = 2.0 * sqrt(M_PION * M_PION + (p_low  * p_low )) / q;
			const double z_high = 2.0 * sqrt(M_PION * M_PION + (p_high * p_high)) / q;

			z = (z_high + z_low) / 2.0;
			ex = y["value"].as<double>();

			// Add uncorrelated systematic and statistical error in quadrature
			const double stat_error = y["errors"][0]["symerror"].as<double>();
			const double  sys_error = y["errors"][1]["symerror"].as<double>();
			const double quad_error = sqrt(stat_error * stat_error + sys_error * sys_error);

			// Systematic relative normalization uncertainty (convert from percentage)
			const std::string s_rnu = y["errors"][2]["symerror"].as<string>();
			const double rnu = LHAPDF::lexical_cast<double>(s_rnu.substr(0, s_rnu.length()-1)) / 100.0;

			// Append data
			z_lows.push_back(z_low);
			z_highs.push_back(z_high);
			uncorr_uncertainty.push_back(quad_error);
			norm_uncertainty.push_back(rnu);
		});

	std::cout << "sqrt_s: " << sqrt_s << std::endl;

	// Gather theory data
	std::vector<double> theoretical, theoretical_dev;
	calculate_theoretical(pdf_name, zs, sqrt_s, theoretical, theoretical_dev, theory_xsec);
	
	//theoretical = calculate_cross_sections(pdf_name, 0, zs, sqrt_s, theory_xsec);
	//theoretical_dev = std::vector<double>(theoretical.size(), 0);

	// Write results to a file for plotting
	std::ofstream outfile("table2.data");
	outfile << "# z\texp\terr\tthr" << std::endl;

	std::cout << std::setw(12) << "z" << " "
		  << std::setw(12) << "1/sig dsig/dp"
		  << std::endl;
	for (int i = 0; i < zs.size(); i++)
	{
		// Write to file
		outfile << z_lows[i] << "\t"
			<< z_highs[i] << "\t"
			<< experimental[i] << "\t" << uncorr_uncertainty[i] << "\t"
			<< theoretical[i] << "\t" << sqrt(theoretical_dev[i]) << std::endl;

		// Print to screen
		std::cout << std::setprecision(5)
			  << std::setw(12)
			  << zs[i] << " "
			  << std::setprecision(4)
			  << std::setw(12)
			  << theoretical[i]
			  << std::endl;
	}
	outfile.close();

	// Filter datasets by z_cuts
	const auto   theoretical_filtered = filter_bounds(zs,        theoretical, Z_MIN, Z_MAX);
	const auto  experimental_filtered = filter_bounds(zs,       experimental, Z_MIN, Z_MAX);
	const auto   norm_uncert_filtered = filter_bounds(zs,   norm_uncertainty, Z_MIN, Z_MAX);
	const auto uncorr_uncert_filtered = filter_bounds(zs, uncorr_uncertainty, Z_MIN, Z_MAX);

	// Compute chi squared
	double chi2 = chi_squared(theoretical_filtered, experimental_filtered,
				  [&experimental_filtered, &norm_uncert_filtered, &uncorr_uncert_filtered]
				   (int i, int j) -> double {
				   	double cov = (norm_uncert_filtered[i] * experimental_filtered[i]) *
						     (norm_uncert_filtered[j] * experimental_filtered[j]);
					if (i == j)
						cov += uncorr_uncert_filtered[i] * uncorr_uncert_filtered[i];
					return cov;
				   });

	std::cout << "Chi Squared: " << chi2 << std::endl;
	std::cout << "Chi squared per point: " << chi2 / ((double) experimental_filtered.size()) << std::endl;

	return 0;
}
