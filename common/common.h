#ifndef COMMON_H
#define COMMON_H


#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"
#include "yaml-cpp/yaml.h"

#include <string>
#include <memory>
#include <cmath>

/*
 * Common definitions and boilerplate for MAPFF 1.0 fragmentation function.
 * Be careful if using an updated FF.
 */
namespace common
{
	// Particle masses (GeV)
	const double M_CHARM	= 1.51;
	const double M_BOTTOM	= 4.92;
	const double M_TOP	= 175.;
	const double M_MUON	= 0.105;
	const double M_TAU	= 1.777;
	const double M_Z	= 91.1876;
	const double M_PION	= 0.139570;

	// Coupling initial scales
	const double ALPHA_S_MZ	= 0.118;
	const double ALPHA_TAU	= 1./128.;

	// sin^2 of theta Weinberg
	const double SIN2_TW	= 0.23122;

	// Perturbative order for MAPFF (NLO)
	constexpr static const int PTO = 1;

	/*
	 * Initialize APFEL with the given PDF set
	 */
	void APFEL_init(const std::string& pdf_set)
	{
		//APFEL::EnableWelcomeMessage(false);

		// Set parameters
		//APFEL::SetMassScheme("ZM-VFNS");
		APFEL::SetAlphaQCDRef(ALPHA_S_MZ, M_Z);
		APFEL::SetAlphaQEDRef(ALPHA_TAU, M_TAU);
		APFEL::SetSin2ThetaW(SIN2_TW);
		APFEL::SetPoleMasses(M_CHARM, M_BOTTOM, M_TOP);
		APFEL::SetTauMass(M_TAU);
		APFEL::SetZMass(M_Z);
		APFEL::SetAlphaEvolution("exact");
		APFEL::SetPerturbativeOrder(PTO);

		// Set grid parameters (cubic splines)
		APFEL::SetNumberOfGrids(3);
		APFEL::SetGridParameters(1, 100, 3, 1e-2);
		APFEL::SetGridParameters(2,  60, 3, 5e-1);
		APFEL::SetGridParameters(3,  50, 3, 8e-1);

		// Set PDF set
		APFEL::SetPDFSet(pdf_set);
		APFEL::SetReplica(0);

		// Initialize DIS
		APFEL::InitializeAPFEL_DIS();
	}

	// Set timelike neutral current for SIA
	void APFEL_SetSIA()
	{
		APFEL::SetTimeLikeEvolution(true);
		APFEL::SetProcessDIS("NC");
	}

	// Set APFEL replica
	void APFEL_Replica(int replica)
	{
		APFEL::SetReplica(replica);
	}

	// Initialize APFEL at a energy scale q (GeV)
	void APFEL_DIS_init(double q)
	{
		APFEL::ComputeStructureFunctionsAPFEL(q, q);
	}

	/*
	 * Tags for which order to compute F2 and SIA total cross sections to.
	 */
	const int LIGHT = 0;
	const int CHARM = 1;
	const int BOTTOM = 2;
	const int TOTAL = 3;

	/*
	 * Compute the sum of F2 functions from light to quark (LIGHT, CHARM, BOTTOM, TOTAL).
	 * Parameters are given via template for the compiler to optimize out.
	 */
	template<int quark>
	double F2(double z)
	{
		double total = 0.0;

		if (LIGHT  <= quark) total += APFEL::F2light(z);
		if (CHARM  <= quark) total += APFEL::F2charm(z);
		if (BOTTOM <= quark) total += APFEL::F2bottom(z);
		if (TOTAL  <= quark) total += APFEL::F2top(z);

		return total;
	}

	/*
	 * Compute the sum of FL functions from light to quark (see F2).
	 */
	template<int quark>
	double FL(double z)
	{
		double total = 0.0;

		if (LIGHT  <= quark) total += APFEL::FLlight(z);
		if (CHARM  <= quark) total += APFEL::FLcharm(z);
		if (BOTTOM <= quark) total += APFEL::FLbottom(z);
		if (TOTAL  <= quark) total += APFEL::FLtop(z);

		return total;
	}

	/*
	 * Compute the sum of SIA total cross sections from light to quark at the given perturbative order.
	 * Parameters are given via template for the compiler to optimize out. 
	 */
	template<int quark, int pto = PTO>
	double SIATotal(double q)
	{
		double total = 0.0;

		if (LIGHT  <= quark) total += APFEL::GetSIATotalCrossSection(pto, q, "light");
		if (CHARM  <= quark) total += APFEL::GetSIATotalCrossSection(pto, q, "charm");
		if (BOTTOM <= quark) total += APFEL::GetSIATotalCrossSection(pto, q, "bottom");
		if (TOTAL  <= quark) total += APFEL::GetSIATotalCrossSection(pto, q, "top");

		return total;
	}

	/*
	 * Computes cross sections for the given PDF and replica number. The cross-sections are
	 * evaluated at each z in zs, with a center of mass energy q.
	 * initialize_func is run before setting the APFEL parameters (i.e. setting SIA evaluation).
	 * cross_section(double z, double q) should compute the cross section for the given z and q.
	 */
	template<typename XSec>
	std::vector<double> calculate_cross_sections(const std::string& pdf_name, int replica,
						     const std::vector<double>& zs, double q,
						     XSec cross_section)
	{
		APFEL_Replica(replica);
		APFEL_DIS_init(q);		// Need to recompute DIS functions for every replica

		std::vector<double> theoretical;

		for (double z : zs)
			theoretical.push_back(cross_section(z, q));

		return theoretical;
	}

	/*
	 * Computes the average over a function T(int) over the range [0,size)
	 */
	template<typename T>
	double average(T t, int size)
	{
		double total = 0;
		for (int i = 0; i < size; i++)
			total += t(i);
		return total / ((double) size);
	}

	/*
	 * Gives the number of replicas for a given PDF set
	 */
	int getReplicaCount(const std::string& pdf_name)
	{
		std::unique_ptr<LHAPDF::PDFInfo> pdf_info(LHAPDF::mkPDFInfo(pdf_name, 0));
		return pdf_info->get_entry_as<int>("NumMembers")-1;
	}

	/*
	 * For a given PDF set, calculates the cross sections for the given z and sqrt_s values.
	 * The central average is given by theory_avg (should just match cross-sections directly from
	 * replica 0). The variance is given by theory_var.
	 */
	template<typename XSec>
	void calculate_theoretical(const std::string& pdf_name, const std::vector<double>& zs, double q,
				   std::vector<double>& theory_avg, std::vector<double>& theory_var,
				   XSec cross_section)
	{
		std::vector<std::vector<double>> theory_points;

		const int replica_num = getReplicaCount(pdf_name);
		for (int i = 1; i <= replica_num; i++)
			theory_points.push_back(calculate_cross_sections(pdf_name, i, zs, q, cross_section));

		theory_avg = std::vector<double>(zs.size(), 0);
		theory_var = std::vector<double>(zs.size(), 0);

		// Calculate average
		for (int i = 0; i < zs.size(); i++)
			theory_avg[i] = average([theory_points, i] (int j) {
				return theory_points[j][i];
			}, theory_points.size());

		// Calculate variance
		for (int i = 0; i < zs.size(); i++)
			theory_var[i] = average([theory_points, theory_avg, i] (int j) {
				double diff = theory_points[j][i] - theory_avg[i];
				return diff * diff;
			}, theory_points.size());
	}

	/*
	 * Determines if a is in the bounds [low, high]
	 */
	bool in_bounds(double a, double low, double high)
	{
		return (a >= low) and (a <= high);
	}

	/*
	 * Parses a YAML datafile identified by filename.
	 * z values are stored in zs
	 * Measured values corresponding to zs are stored in experimental
	 * The parse(const YAML::Node&, const YAML::Node&, double, double&, double&) -> void
	 *     parses a single line of the datafile. Extracting any errors / extra data should
	 *     be done by capturing external variables in the supplied functor.
	 */
	template<typename EnergyFunc, typename ParseFunc>
	double parse_datafile_sia(const std::string& filename, std::vector<double>& zs,
			          std::vector<double>& experimental, EnergyFunc energy, ParseFunc parse)
	{
		auto table = YAML::LoadFile(filename);

		const double q = energy(table);

		const auto xs = table["independent_variables"][0]["values"];
		const auto ys = table[  "dependent_variables"][0]["values"];

		zs           = std::vector<double>(xs.size(), 0);
		experimental = std::vector<double>(xs.size(), 0);

		for (int i = 0; i < xs.size(); i++)
			parse(xs[i], ys[i], q, zs[i], experimental[i]);

		return q;
	}

	/*
	 * Returns a vector of arr but filtered out datapoints which are outside
	 * the bounds [z_min, z_max]
	 */
	std::vector<double> filter_bounds(const std::vector<double>& zs, const std::vector<double>& arr, double z_min, double z_max)
	{
		std::vector<double> out;

		for (int i = 0; i < zs.size(); i++)
			if (in_bounds(zs[i], z_min, z_max))
				out.push_back(arr[i]);

		return out;
	}

	/*
	 * Computes the chi-squared for the given theoretical th values and experimental ex values.
	 * cov(int, int) is the covariance matrix which returns the covariance between points i and j.
	 */
	template<typename CovMat>
	double chi_squared(const std::vector<double>& th, const std::vector<double>& ex, CovMat cov)
	{
		double chi2 = 0;
		for (int i = 0; i < ex.size(); i++)
		{
			const double left = th[i] - ex[i];

			for (int j = 0; j < ex.size(); j++)
			{
				const double right = th[j] - ex[j];
				
				chi2 += left * right / cov(i,j);
			}
		}

		return chi2;
	}

	/*
	 * Checks if two doubles are equal within an epsilon neighborhood.
	 */
	const double CMP_EPSILON = 1e-6;
	bool f_eq(double a, double b)
	{
		return fabs(a - b) < CMP_EPSILON;
	}
};

#endif
