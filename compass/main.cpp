#include <apfel/apfelxx.h>
#include <apfel/SIDIS.h>
#include "LHAPDF/LHAPDF.h"
#include "../common/common.h"
#include "yaml-cpp/yaml.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <memory>
#include <array>
#include <functional>

/* Variable Naming Convention
 *	q, Q	Neutral current virtuality
 *	x	Incoming parton momentum fraction
 *	y	Energy transfer
 *	z	Outgoing parton momentum fraction
 *	M	Multiplicity
 *	s	Center of Mass energy
 *
 *	pdf	Parton Distribution Function
 *	ff	Fragmentation Function
 *
 *	F2	Second structure function (of the proton)
 *	FL	Longitudinal structure function (of the proton)
 *	C____	Coefficient functions (i.e. CL2qg is 2nd order quark-gluon
 *			longitudinal coefficient function)
 *
 *	DIS	Deep Inelastic Scattering
 *	SIDIS	Semi-Inclusive Deep Inelastic Scattering
 */

/*	General program structure
 * 
 * The purpose of the program is to calculate the pion+ multiplicity from theoretical parton
 * distribution function (PDF) and fragmentation function (FF) models, and compare them to the
 * results from COMPASS experiments. Directly computing multiplicity is extremely slow because
 * it requires a large amount of integrals and convolutions. APFEL++ is a library which speeds
 * up these computations by tabulating computations on a grid and interpolating between them.
 * We have the following grids:
 * 	x-space		PDF parameter
 * 	y-space		Only used in cross section calculation
 * 	z-space		FF parameter
 * 	Q-space		energy parameter (PDF, FF, and cross section)
 * 
 * Important objects are as follows:
 * 	TabulateObject		Uses factory functions (that we provide) to build a grid in Q-
 * 				space and x-space. If parameterized by a DoubleObject, also
 * 				builds a grid in z-space.
 *	DoubleObject		A linear combination of a terms. Each term is a product of
 *				a function in x-space and a function in z-space.
 *	Operator		An operator that has been Mellin transformed over a x- or
 *				z-space.
 *	Distribution		A function that may be evaluated at x or z (i.e. a PDF or
 *				FF respectively).
 * 
 * Using APFEL amounts to building the appropriate factory functions that let TabulateObject
 * generate grids in each variable. In general, these factory functions are parameterized by
 * the energy Q, and sometimes a quark index. 
 */

/* ----- Constants ----- */

const int gluon = 21;	// Gluon PDF index

// Quark mass thresholds
const std::vector<double> Thresholds{0, 0, 0, common::M_CHARM, common::M_BOTTOM, common::M_TOP};

// Basis of real quark indices at an energy q
std::vector<int> real_quark_basis(double const& q, std::vector<double> const& Thrs) {
	return std::vector<int>{1, 2, 3, 4};
	/*
	std::vector<int> quarks;
	for (int i = 0; i < Thrs.size(); i++)
		if (Thrs[i] <= q)
			quarks.push_back(i);
	return quarks;
	*/
}

// Basis of real and anti quark indices
std::vector<int> quark_basis(double const& q, std::vector<double> const& Thrs) {
	auto quarks = real_quark_basis(q, Thrs);
	const int s = quarks.size();
	for (int i = 0; i < s; i++)
		quarks.push_back(-quarks[i]);
	return quarks;
}

// Basis of real and anti quark indices, and the gluon index
std::vector<int> qcd_basis(double const& q, std::vector<double> const& Thrs) {
	auto quarks = quark_basis(q, Thrs);
	quarks.push_back(gluon);
	return quarks;
}

// APFEL's structure function basis (total = 0, light = 1, charm = 2, bottom = 3, top = 4)
std::vector<int> structure_basis(double const& q, std::vector<double> const& Thrs) {
	auto quarks = real_quark_basis(q, Thrs);
	std::vector<int> basis;
	for (int i = 2; i < quarks.size(); i++)
		basis.push_back(i-1);
	return basis;
}


const double rel_eps = 1e-3;	// Relative epsilon for integration accuracy
const int InterDeg = 3;		// Interpolation degree (cubic splines)

// Running couplings as a function of renormalization scale
std::function<double(double const&)> alpha;	// QED
std::function<double(double const&)> alpha_s;	// QCD

/* ----- Bound Structure ----- */

// Close interval
struct Bound
{
	const double min;
	const double max;
};

const Bound QGridBound{1., 10.};	// Bounds for Q-grids

/* ----- Calculate Integration Bounds ----- */

// Calculate Q integration bound
Bound QBound(const Bound& x_bound, const Bound& y_bound, double s)
{
	return Bound{
		sqrt(x_bound.min * y_bound.min * s),
		sqrt(x_bound.max * y_bound.max * s)
	};
}

// Calculate modified x integration bound
Bound ModXBound(const Bound& x_bound, const Bound& y_bound, double q, double s)
{
	return Bound{
		std::max(x_bound.min, q*q / (s * y_bound.max)),
		std::min(x_bound.max, q*q / (s * y_bound.min))
	};
}

/* ----- Cross Section Calculations ----- */

double SIDIS_xsec(double x, double y, double z, double q,
			const apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>>& F2,
			const apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>>& FL)
{
	return ((apfel::FourPi * alpha_s(q) * alpha_s(q)) / (q * q * q)) * ((1.+(1.-y)*(1.-y))*F2.EvaluatexzQ(x, z, q) - y*y*FL.EvaluatexzQ(x, z, q));
}

double DIS_xsec(double x, double y, double q,
			const apfel::TabulateObject<apfel::Distribution>& F2,
			const apfel::TabulateObject<apfel::Distribution>& FL)
{
	return ((apfel::FourPi * alpha_s(q) * alpha_s(q)) / (q * q * q)) * ((1.+(1.-y)*(1.-y))*F2.EvaluatexQ(x, q) - y*y*FL.EvaluatexQ(x, q));
}

/* ----- Multiplicity ----- */

// Integrate over SIDIS bounds
double multiplicity_numerator(double s, double y,
				const Bound& x_bound,
				const Bound& y_bound,
				const Bound& z_bound,
				const apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>>& F2,
				const apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>>& FL)
{
	const Bound q_bound = QBound(x_bound, y_bound, s);

	apfel::Integrator q_integral([&] (double const& q) -> double {
		const Bound mod_x_bound = ModXBound(x_bound, y_bound, q, s);

		apfel::Integrator x_integral([&] (double const& x) -> double {
			apfel::Integrator z_integral([&] (double const& z) -> double {
				
				return SIDIS_xsec(x, y, z, q, F2, FL);
			});

			return z_integral.integrate(z_bound.min, z_bound.max, rel_eps);
		});

		return x_integral.integrate(mod_x_bound.min, mod_x_bound.max, rel_eps);
	});

	return q_integral.integrate(q_bound.min, q_bound.max, rel_eps);
}

// Integrate over DIS bounds
double multiplicity_denominator(double s, double y,
				const Bound& x_bound,
				const Bound& y_bound,
				const Bound& z_bound,
				const apfel::TabulateObject<apfel::Distribution>& F2,
				const apfel::TabulateObject<apfel::Distribution>& FL)
{
	const Bound q_bound = QBound(x_bound, y_bound, s);

	apfel::Integrator q_integral([&] (double const& q) {
		const Bound mod_x_bound = ModXBound(x_bound, y_bound, q, s);

		apfel::Integrator x_integral([&] (double const& x) {
			return DIS_xsec(x, y, q, F2, FL);
		});

		return x_integral.integrate(mod_x_bound.min, mod_x_bound.max, rel_eps);
	});

	return q_integral.integrate(q_bound.min, q_bound.max, rel_eps) * (z_bound.max - z_bound.min);
}

double multiplicity(double s, double y,
			const Bound& x_bound,
			const Bound& y_bound,
			const Bound& z_bound,
			const apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>>& sidis_F2,
			const apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>>& sidis_FL,
			const apfel::TabulateObject<apfel::Distribution>& dis_F2,
			const apfel::TabulateObject<apfel::Distribution>& dis_FL)
{
	return multiplicity_numerator  (s, y, x_bound, y_bound, z_bound, sidis_F2, sidis_FL)
	     / multiplicity_denominator(s, y, x_bound, y_bound, z_bound,   dis_F2,   dis_FL);
}

/* ----- FILE IO ----- */

void add_if_missing(std::vector<double>& values, double x)
{
	for (int i = 0; i < values.size(); i++)
		if (common::f_eq(x, values[i]))
			return;
	values.push_back(x);
}

// Compute bins for the given independent variables
void write_bins_to_file(const std::string& name, const YAML::Node& values, std::ofstream& outfile)
{
	std::vector<double> bounds;

	for (int i = 0; i < values.size(); i++)
	{
		add_if_missing(bounds, values[i]["low" ].as<double>());
		add_if_missing(bounds, values[i]["high"].as<double>());
	}

	outfile << "# " << name << " bins" << std::endl;
	for (int i = 0; i < bounds.size(); i++)
		outfile << bounds[i] << " ";
	outfile << std::endl;
}


/* ----- MAIN ----- */

int main(int argc, char** argv)
{
	int ff_replica_num = 0;
	if (argc == 2)
		ff_replica_num = LHAPDF::lexical_cast<int>(argv[1]);

	// Define interpolation grids for PDFs and FFs
	const apfel::Grid pdf_grid{
		{apfel::SubGrid{100, 1e-3, InterDeg},
		 apfel::SubGrid{ 60, 1e-1, InterDeg},
		 apfel::SubGrid{ 50, 6e-1, InterDeg},
		 apfel::SubGrid{ 50, 8e-1, InterDeg}}};
	
	const apfel::Grid ff_grid{
		{apfel::SubGrid{100, 1e-2, InterDeg},
		 apfel::SubGrid{ 60, 5e-1, InterDeg},
		 apfel::SubGrid{ 50, 8e-1, InterDeg}}};
	
	// Calculate running couplings
	apfel::AlphaQED alpha_qed(common::ALPHA_TAU , common::M_TAU, {0, common::M_MUON, common::M_TAU}, Thresholds, common::PTO);
	apfel::AlphaQCD alpha_qcd(common::ALPHA_S_MZ, common::M_Z  , Thresholds, common::PTO);
	apfel::TabulateObject<double> alpha_qeds{alpha_qed, 100, 0.9, 1001, InterDeg};
	apfel::TabulateObject<double> alpha_qcds{alpha_qcd, 100, 0.9, 1001, InterDeg};

	alpha   = [&] (double const& mu) -> double { return alpha_qeds.Evaluate(mu); };
	alpha_s = [&] (double const& mu) -> double { return alpha_qcds.Evaluate(mu); };

	// Initialize PDF and FF
	const LHAPDF::PDFSet pdf_set("NNPDF31_nlo_pch_as_0118");
	const LHAPDF::PDFSet  ff_set("MAPFF10NLOPIp");

	std::shared_ptr<LHAPDF::PDF> pdf_member(pdf_set.mkPDF(0));
	std::shared_ptr<LHAPDF::PDF>  ff_member( ff_set.mkPDF(ff_replica_num));

	// Use isospin symmetry to deduce neutron PDF (swap u and d quarks)
	const auto neutron_isospin = [] (int i) -> int {
		switch(i)
		{
			case -2: return -1;
			case -1: return -2;
			case 1: return 2;
			case 2: return 1;
			default: return i;
		}
	};

	// Construct PDF for Lithium-Deuterium target (average over proton and neutron pdfs)
	const auto  proton_pdf = [&] (int i, double x, double q) -> double { return pdf_member->xfxQ(i, x, q); };
	const auto neutron_pdf = [&] (int i, double x, double q) -> double { return pdf_member->xfxQ(neutron_isospin(i), x, q); };
	const double lid_mass = 6.94 + 2.01;
	const double lid_anum = 3. + 1.;
	const auto lid_pdf = [&] (int i, double x, double q) -> double { return (lid_anum * proton_pdf(i,x,q) + (lid_mass - lid_anum) * neutron_pdf(i,x,q)) / lid_mass; };

	// Construct a map of factory functions which, given an energy scale q, create a product of distributions
	// 	First index - pdf
	// 	Second index - ff
	std::map<int, std::map<int, std::function<apfel::DoubleObject<apfel::Distribution>(double const&)>>> dists;
	for (int i : qcd_basis(QGridBound.max, Thresholds))
	{
		dists[i] = std::map<int, std::function<apfel::DoubleObject<apfel::Distribution>(double const&)>>();
		for (int j : qcd_basis(QGridBound.max, Thresholds))
		{
			dists[i][j] = [&, pdf_flavor=i, ff_flavor=j] (double const& q) -> apfel::DoubleObject<apfel::Distribution> {
				return apfel::DoubleObject<apfel::Distribution>{{
					{1.,
					//apfel::Distribution{pdf_grid, [=] (double const& x) -> double { return pdf_member->xfxQ(pdf_flavor, x, q); }},
					apfel::Distribution{pdf_grid, [=] (double const& x) -> double { return lid_pdf(pdf_flavor, x, q); }},
					apfel::Distribution{ ff_grid, [=] (double const& x) -> double { return ff_member->xfxQ(ff_flavor, x, q); }}}
				}};
			};
		}
	}

	// Initialize SIDIS coefficient functions
	auto SIDIS_Cs = apfel::InitializeSIDIS(pdf_grid, ff_grid);
	
	// Use to build structure functions from coefficient functions
	const auto Fi = [&] (apfel::DoubleObject<apfel::Operator> const& C0qq,
			     apfel::DoubleObject<apfel::Operator> const& C1qq,
			     apfel::DoubleObject<apfel::Operator> const& C1qg,
			     apfel::DoubleObject<apfel::Operator> const& C1gq)
			{
				return [&] (double const& q) -> apfel::DoubleObject<apfel::Distribution> {
					apfel::DoubleObject<apfel::Distribution> total;
					const double cp = alpha_s(q) / apfel::FourPi;
					for (int quark : quark_basis(q, Thresholds))
					{
						const double eq2 = (abs(quark) % 2 == 0) ? 4./9. : 1./9.;
						total += eq2 * ((C0qq + cp * C1qq) * dists[quark][quark](q)
								      + cp * C1qg  * dists[gluon][quark](q)
								      + cp * C1gq  * dists[quark][gluon](q));
					}
					return total;
				};
			};

	// Tabulate SIDIS structure functions
	const apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>> sidis_F2{
		Fi(SIDIS_Cs.C20qq, SIDIS_Cs.C21qq, SIDIS_Cs.C21qg, SIDIS_Cs.C21gq),
		100, QGridBound.min, QGridBound.max, InterDeg, Thresholds
	};

	const apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>> sidis_FL{
		Fi(apfel::DoubleObject<apfel::Operator>(), SIDIS_Cs.CL1qq, SIDIS_Cs.CL1qg, SIDIS_Cs.CL1gq),
		100, QGridBound.min, QGridBound.max, InterDeg, Thresholds
	};

	// Rewrite PDFS in DIS structure function basis
	const auto DIS_pdfs = [&] (double const& x, double const& Q) -> std::map<int, double> {
		std::map<int, double> results;

		results[apfel::DISNCBasis::CNS] = 0;
		for (int quark : real_quark_basis(Q, Thresholds))
			results[apfel::DISNCBasis::CNS] += lid_pdf(quark, x, Q) - lid_pdf(-quark, x, Q);

		results[apfel::DISNCBasis::CS] = 0;
		for (int quark : real_quark_basis(Q, Thresholds))
			results[apfel::DISNCBasis::CS] += lid_pdf(quark, x, Q) + lid_pdf(quark, x, Q);

		results[apfel::DISNCBasis::CG] = lid_pdf(gluon, x, Q);

		return results;
	};

	// Calculate electroweak charges
	const auto fBq = [] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };

	// Have APFEL initialize and build the structure functions for us
	const auto F2Obj = apfel::InitializeF2NCObjectsZM(pdf_grid, Thresholds);
	const auto FLObj = apfel::InitializeFLNCObjectsZM(pdf_grid, Thresholds);

	const auto F2_DIS_map = apfel::BuildStructureFunctions(F2Obj, DIS_pdfs, common::PTO, alpha_s, fBq);
	const auto FL_DIS_map = apfel::BuildStructureFunctions(FLObj, DIS_pdfs, common::PTO, alpha_s, fBq);

	// Sum over the structure basis
	const auto sum_dis_funcs = [&] (std::map<int, apfel::Observable<>> const& Obj) {
		return [&] (double const& q) -> apfel::Distribution {
			apfel::Distribution dist{pdf_grid, [] (double const& x) -> double { return 0; }};
			for (int i : structure_basis(q, Thresholds))
				dist += Obj.at(i).Evaluate(q);
			return dist;
		};
	};

	// Tabulate DIS structure functions
	const apfel::TabulateObject<apfel::Distribution> dis_F2 {sum_dis_funcs(F2_DIS_map), 100, QGridBound.min, QGridBound.max, InterDeg, Thresholds};
	const apfel::TabulateObject<apfel::Distribution> dis_FL {sum_dis_funcs(FL_DIS_map), 100, QGridBound.min, QGridBound.max, InterDeg, Thresholds};

	// Read data file
	std::cout << "Writing results..." << std::endl;
	auto table = YAML::LoadFile("data1.yaml");
	const auto xs  = table["independent_variables"][0]["values"];
	const auto ys  = table["independent_variables"][1]["values"];
	const auto q2s = table["independent_variables"][2]["values"];
	const auto zs  = table["independent_variables"][3]["values"];

	const auto Ms  = table[  "dependent_variables"][0]["values"];

	// Write data to file
	std::ofstream outfile("results" + LHAPDF::lexical_cast<std::string>(ff_replica_num) + ".data");
	
	// Write bins to the file
	write_bins_to_file("x", xs, outfile);
	write_bins_to_file("y", ys, outfile);
	write_bins_to_file("z", zs, outfile);

	// Write data points to file
	outfile << "# data x y z q2 M Mt Merr" << std::endl;
	std::cout << "x\ty\tz\tM\tMt\n";
	for (int i = 0; i < xs.size(); i++)
	{
		// Read this point's values
		const double x  =  xs[i]["value"].as<double>();
		const double y  =  ys[i]["value"].as<double>();
		const double z  =  zs[i]["value"].as<double>();
		const double q2 = q2s[i]["value"].as<double>();
		const double M  =  Ms[i]["value"].as<double>();

		const double Mstat = Ms[i]["errors"][0]["symerror"].as<double>();
		const double Msys  = Ms[i]["errors"][1]["symerror"].as<double>();

		const double Merr = sqrt(Mstat * Mstat + Msys * Msys);

		// Read integration bounds
		const Bound x_bound{xs[i]["low"].as<double>(), xs[i]["high"].as<double>()};
		const Bound y_bound{ys[i]["low"].as<double>(), ys[i]["high"].as<double>()};
		const Bound z_bound{zs[i]["low"].as<double>(), zs[i]["high"].as<double>()};

		// Calculate COM energy
		const double s = q2 / (x * y);

		// Calculate multiplicity
		const double Mt = multiplicity(s, y, x_bound, y_bound, z_bound,
					       sidis_F2, sidis_FL, dis_F2, dis_FL);

		// Write to file
		outfile << x  << "\t"
			<< y  << "\t"
			<< z  << "\t"
			<< q2 << "\t"
			<< M  << "\t"
			<< Mt << "\t"
			<< Merr << std::endl;

		// Log to the screen
		std::cout << x << "\t" << y << "\t" << z << "\t" << M << "\t" << Mt << std::endl;
	}

	return 0;
}
