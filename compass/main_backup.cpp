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
#include <chrono>
#include <functional>

const int QQ = 0;
const int QG = 1;
const int GQ = 2;

const int gluon = 21;
const std::vector<int> flavor_basis{-5, -4, -3, -2, -1, gluon, 1, 2, 3, 4, 5};
const std::vector<int> quark_basis{-5, -4, -3, -2, -1, 1, 2, 3, 4, 5};
const std::vector<int> real_quark_basis{1,2,3,4,5};
const std::vector<double> Thresholds{0, 0, 0, common::M_CHARM, common::M_BOTTOM, common::M_TOP};

std::function<double(double const&)> alpha;
std::function<double(double const&)> alpha_s;

using coeff_arr_type = std::array<std::map<int, std::function<apfel::DoubleObject<apfel::Operator>(double const&)>>, 3>;
using   pdf_map_type = std::map<int, std::map<int, std::function<apfel::DoubleObject<apfel::Distribution>(double const&)>>>;
using dis_func_type = std::function<double(double const&, double const&)>;
using evpdf_map_type = std::map<int, std::map<int, apfel::DoubleObject<apfel::Distribution>>>;
using evcoeff_arr_type = std::array<std::map<int, apfel::DoubleObject<apfel::Operator>>, 3>;

double squared_charge(int quark)
{
	if ((std::abs(quark) % 2) == 0) return -1./3.;
	else return 2./3.;
}

double Fi(double x, double z, const evcoeff_arr_type& Cs, const evpdf_map_type& Fs)
{
	double total = 0;
	for (const auto& i : quark_basis)
	{
		// convert to doubleobject<op> *= doubleobject<dist>
		double quark_total = 0;
		quark_total += (
				Cs[QQ].at(i) * Fs.at(i).at(i)
			      + Cs[QG].at(i) * Fs.at(gluon).at(i)
			      + Cs[GQ].at(i) * Fs.at(i).at(gluon)
			      ).Evaluate(x, z);

		quark_total *= squared_charge(i);
		total += quark_total;
	}

	return x * total;
}

double diff_xsec(double x, double y, double z, double q,
		 const evcoeff_arr_type& C2, const evcoeff_arr_type& CL,
		 const evpdf_map_type& pdfs)
{
	const double q2 = q * q;
	return (apfel::FourPi * alpha(q2) * alpha(q2) / (x * q2 * q)) * ( (1. + (1.-y)*(1.-y)) * Fi(x,z,C2,pdfs) - y*y * Fi(x,z,CL,pdfs));
}

const double rel_eps = 1e-1;
double multiplicity_numerator(double s, double y, 
				double x_min, double x_max,
				double y_min, double y_max,
				double z_min, double z_max,
				const coeff_arr_type& C2, const coeff_arr_type& CL,
				const pdf_map_type& pdfs)
{
	std::cout << "multiplicity numerator..." << std::endl;
	apfel::Timer t;
	const double q_min = sqrt(x_min * y_min * s);
	const double q_max = sqrt(x_max * y_max * s);

	apfel::Integrator q_integral([&] (double const& q) {
		std::cout << "\tq_integral..." << std::endl;
		apfel::Timer t1;
		x_min = std::max(x_min, q*q / s / y_max);
		x_min = std::min(x_max, q*q / s / y_min);

		evpdf_map_type evaluated_pdfs;
		for (const auto& [i, ff_maps] : pdfs)
		{
			evaluated_pdfs[i] = std::map<int, apfel::DoubleObject<apfel::Distribution>>();
			for (const auto& [j, f] : ff_maps)
			{
				evaluated_pdfs[i][j] = f(q);
			}
		}

		evcoeff_arr_type evaluated_C2, evaluated_CL;
		for (int i = 0; i < 3; i++)
		{
			evaluated_C2[i] = std::map<int, apfel::DoubleObject<apfel::Operator>>();
			for (const auto& [j, C] : C2[i])
				evaluated_C2[i][j] = C(q);

			evaluated_CL[i] = std::map<int, apfel::DoubleObject<apfel::Operator>>();
			for (const auto& [j, C] : CL[i])
				evaluated_CL[i][j] = C(q);
		}

		apfel::Integrator x_integral([&] (double const& x) {
			//std::cout << "\t\tx_integral..." << std::endl;
			//apfel::Timer t2;
			
			apfel::Integrator z_integral([&] (double const& z) {
				const double result = diff_xsec(x, y, z, q, evaluated_C2, evaluated_CL, evaluated_pdfs);
				//std::cout << "(" << x << ", " << z << ") " << result << std::endl;
				return result;
			});

			const double result = z_integral.integrate(z_min, z_max, rel_eps);
			//t2.stop();
			return result;
		});

		const double result = x_integral.integrate(x_min, x_max, rel_eps);
		t1.stop();
		return result;
	});

	const double result = q_integral.integrate(q_min, q_max, rel_eps);
	t.stop();
	return result;
}

double dis_xsec(double x, double y, double q, dis_func_type F2_dis, dis_func_type FL_dis)
{
	return (apfel::FourPi * alpha(q) * alpha(q) / (x * q * q * q)) * ((1. + (1.-y)*(1.-y))*F2_dis(x, q) - y*y*FL_dis(x,q));
}

double multiplicity_denominator(double s, double y,
				double x_min, double x_max,
				double y_min, double y_max,
				double z_min, double z_max,
				dis_func_type F2_dis, dis_func_type FL_dis)
{
	const double q_min = sqrt(x_min * y_min * s);
	const double q_max = sqrt(x_max * y_max * s);

	apfel::Integrator q_integral([&] (double const& q) {
		x_min = std::max(x_min, q*q / s / y_max);
		x_max = std::min(x_max, q*q / s / y_min);

		apfel::Integrator x_integral([&] (double const& x) {
			return dis_xsec(x, y, q, F2_dis, FL_dis);
		});

		return x_integral.integrate(x_min, x_max, rel_eps);
	});

	return q_integral.integrate(q_min, q_max, rel_eps) * (z_max - z_min);
}

double multiplicity(double s, double y,
			double x_min, double x_max,
			double y_min, double y_max,
			double z_min, double z_max,
			const coeff_arr_type& C2,
			const coeff_arr_type& CL,
			std::function<double(double const&, double const&)> F2_dis,
			std::function<double(double const&, double const&)> FL_dis,
			const pdf_map_type& pdfs)
{
	std::cout << "NUMERATOR..." << std::endl;
	const double num =   multiplicity_numerator(s, y, x_min, x_max, y_min, y_max, z_min, z_max, C2, CL, pdfs);
	//const double num = diff_xsec(x_min, y_min, z_min, sqrt(x_min * y_min * s), C2, CL, pdfs);
	//const double num = dis_xsec(x_min, y_min, sqrt(s * x_min * y_min), F2_dis, FL_dis);
	std::cout << "\tnum = " << num << std::endl;

	std::cout << "DENOMINATOR..." << std::endl;
	const double den = multiplicity_denominator(s, y, x_min, x_max, y_min, y_max, z_min, z_max, F2_dis, FL_dis);
	//const double den = 1;
	std::cout << "\tden = "<< den << std::endl;

	return num / den;
}

void add_if_missing(std::vector<double>& values, double x)
{
	for (int i = 0; i < values.size(); i++)
		if (common::f_eq(x, values[i]))
			return;
	values.push_back(x);
}

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

std::chrono::time_point<std::chrono::high_resolution_clock> start_timer()
{
	return std::chrono::high_resolution_clock::now();
}

std::chrono::time_point<std::chrono::high_resolution_clock> end_timer()
{
	return std::chrono::high_resolution_clock::now();
}

template<typename T>
void print_time_results(T start, T end)
{
	const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "Duration: " << duration.count() << "ms" << std::endl;
}

int main( void )
{
	std::cout << "Reading YAML file..." << std::endl;
	auto start = start_timer();
	YAML::Node table = YAML::LoadFile("data2.yaml");
	auto end = end_timer();
	print_time_results(start, end);


	// ------- SIDIS -----------

	const apfel::Grid pdf_grid{
		{apfel::SubGrid{100, 1e-3, 3},
		 apfel::SubGrid{ 60, 1e-1, 3},
		 apfel::SubGrid{ 50, 6e-1, 3},
		 apfel::SubGrid{ 50, 8e-1, 3}}};

	const apfel::Grid ff_grid{
		{apfel::SubGrid{100, 1e-2, 3},
		 apfel::SubGrid{ 60, 5e-1, 3},
		 apfel::SubGrid{ 50, 8e-1, 3}}};

	// Couplings
	apfel::AlphaQED alpha_qed(common::ALPHA_TAU , common::M_TAU, {0, common::M_MUON, common::M_TAU}, Thresholds, common::PTO);
	apfel::AlphaQCD alpha_qcd(common::ALPHA_S_MZ, common::M_Z  , Thresholds, common::PTO);
	apfel::TabulateObject<double> alpha_qeds{alpha_qed, 100, 0.9, 1001, 3};
	apfel::TabulateObject<double> alpha_qcds{alpha_qcd, 100, 0.9, 1001, 3};

	alpha   = [&] (double const& mu) -> double { return alpha_qeds.Evaluate(mu); };
	alpha_s = [&] (double const& mu) -> double { return alpha_qcds.Evaluate(mu); };

	// Initialize PDF
	const LHAPDF::PDFSet pdf_set("NNPDF31_nlo_pch_as_0118");
	const LHAPDF::PDF* const pdf_member = pdf_set.mkPDF(0);

	// Initialize FF
	const LHAPDF::PDFSet ff_set("MAPFF10NLOPIp");
	const LHAPDF::PDF* const ff_member = ff_set.mkPDF(0);

	// Build <PDF, FF> double object
	pdf_map_type pdfs;

	for (int i : flavor_basis)
	{
		pdfs.insert({i, std::map<int, std::function<apfel::DoubleObject<apfel::Distribution>(double const&)>>()});

		for (int j : flavor_basis)
		{
			pdfs.at(i).insert({j, [&, pdf_flavor=i, ff_flavor=j] (double const& q) -> apfel::DoubleObject<apfel::Distribution> {
				apfel::DoubleObject<apfel::Distribution> double_dist{{
					{1., 
					apfel::Distribution{pdf_grid, [=] (double const& x) -> double { return pdf_member->xfxQ(pdf_flavor, x, q);}},
					apfel::Distribution{ ff_grid, [=] (double const& x) -> double { return  ff_member->xfxQ( ff_flavor, x, q);}}}
					}};
				return double_dist;
				}});
		}
	}

	// Initialize Coefficient Functions
	coeff_arr_type C2s;
	coeff_arr_type CLs;
	auto SIDIS_Cs = apfel::InitializeSIDIS(pdf_grid, ff_grid);
	for (int flavor : quark_basis)
	{
		C2s[QQ][flavor] = [=] (double const& q) {
			const double cp = alpha_s(q) / apfel::FourPi;
			return SIDIS_Cs.C20qq + cp * SIDIS_Cs.C21qq;
		};

		C2s[QG][flavor] = [=] (double const& q) {
			const double cp = alpha_s(q) / apfel::FourPi;
			return cp * SIDIS_Cs.C21qg;
		};

		C2s[GQ][flavor] = [=] (double const& q) {
			const double cp = alpha_s(q) / apfel::FourPi;
			return cp * SIDIS_Cs.C21gq;
		};

		CLs[QQ][flavor] = [=] (double const& q) {
			const double cp = alpha_s(q) / apfel::FourPi;
			return cp * SIDIS_Cs.CL1qq;
		};

		CLs[QG][flavor] = [=] (double const& q) {
			const double cp = alpha_s(q) / apfel::FourPi;
			return cp * SIDIS_Cs.CL1qg;
		};

		CLs[GQ][flavor] = [=] (double const& q) {
			const double cp = alpha_s(q) / apfel::FourPi;
			return cp * SIDIS_Cs.CL1gq;
		};
	}

	// Initialize DIS Functions
	const auto DIS_pdfs = [&] (double const& z, double const& Q) -> std::map<int, double> {
		std::map<int, double> results;
		for (const auto& i : flavor_basis)
			results[i] = pdf_member->xfxQ(i, z, Q);
		return results;
	};

	const auto F2Obj = InitializeF2NCObjectsZM(pdf_grid, Thresholds);
	const auto FLObj = InitializeFLNCObjectsZM(pdf_grid, Thresholds);

	const auto fBq = [=] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };

	const auto F2DIS = BuildStructureFunctions(F2Obj, DIS_pdfs, common::PTO, alpha_s, fBq);
	const auto FLDIS = BuildStructureFunctions(FLObj, DIS_pdfs, common::PTO, alpha_s, fBq);

	apfel::TabulateObject<apfel::Distribution> F2DIStable{
		[&] (double const& Q) -> apfel::Distribution {
			apfel::Distribution total = F2DIS.at(1).Evaluate(Q);
			for (int i = 0; i <= 5; i++)
				total += F2DIS.at(i).Evaluate(Q);
			return total;
		}, 50, 1, 1000, 3, Thresholds
	};

	apfel::TabulateObject<apfel::Distribution> FLDIStable{
		[&] (double const& Q) -> apfel::Distribution {
			apfel::Distribution total = FLDIS.at(1).Evaluate(Q);
			for (int i = 0; i <= 5; i++)
				total += FLDIS.at(i).Evaluate(Q);
			return total;
		}, 50, 1, 1000, 3, Thresholds
	};

	const auto F2_dis = [&] (double const& x, double const& Q) { return F2DIStable.EvaluatexQ(x, Q); };
	const auto FL_dis = [&] (double const& x, double const& Q) { return FLDIStable.EvaluatexQ(x, Q); };
	
	// ------ Output ------

	std::cout << "Writing results..." << std::endl;
	start = start_timer();
	const auto xs  = table["independent_variables"][0]["values"];
	const auto ys  = table["independent_variables"][1]["values"];
	const auto q2s = table["independent_variables"][2]["values"];
	const auto zs  = table["independent_variables"][3]["values"];

	const auto Ms  = table[  "dependent_variables"][0]["values"];

	// Write data to file
	std::ofstream outfile("results.data");
	
	// Write bins to the file
	write_bins_to_file("x", xs, outfile);
	write_bins_to_file("y", ys, outfile);
	write_bins_to_file("z", zs, outfile);

	// Write data points to file
	outfile << "# data x y z q2 M Mt" << std::endl;
	for (int i = 0; i < xs.size(); i++)
	{
		const double x  =  xs[i]["value"].as<double>();
		const double y  =  ys[i]["value"].as<double>();
		const double z  =  zs[i]["value"].as<double>();
		const double q2 = q2s[i]["value"].as<double>();

		const double M  =  Ms[i]["value"].as<double>();

		const double x_min = xs[i][ "low"].as<double>();
		const double x_max = xs[i]["high"].as<double>();
		const double y_min = ys[i][ "low"].as<double>();
		const double y_max = ys[i]["high"].as<double>();
		const double z_min = zs[i][ "low"].as<double>();
		const double z_max = zs[i]["high"].as<double>();
		const double s = q2 / (x * y);

		const double dMdz = multiplicity(s, y, x_min, x_max, y_min, y_max, z_min, z_max,
						 C2s, CLs, F2_dis, FL_dis, pdfs);

		outfile << x    << "\t"
			<< y    << "\t"
			<< z    << "\t"
			<< q2   << "\t"
			<< M    << "\t"
			<< dMdz << std::endl;
	}
	end = end_timer();
	print_time_results(start, end);

	delete pdf_member, ff_member; 
}
