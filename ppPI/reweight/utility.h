#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <tuple>
#include <vector>
#include <iomanip>
#include <sstream>
#include <numeric>

#include "matrix.h"
#include "unweight.h"

const std::string experiments_filename = "experiments.txt";

using data_list = std::vector<double>;

auto within_cutoff(double cutoff)
{
	return [cutoff] (double x) -> bool { return x >= cutoff; };
}

std::tuple<int, double, std::vector<std::tuple<std::string, std::string>>> read_experiment_file(const std::string& filename = experiments_filename)
{
	std::ifstream infile(filename);
	if (!infile)
	{
		std::cerr << "Error occurred reading experiments file: " << experiments_filename << std::endl;
		std::exit(1);
	}

	int replica_count = 0;
	std::vector<std::tuple<std::string, std::string>> experiment_names;

	double cutoff = 0;

	infile >> replica_count;
	infile >> cutoff;

	std::string line;
	while (std::getline(infile, line))
	{
		if (line.empty()) continue;
		if (line[0] == '#') continue;

		std::stringstream sstream;
		sstream << line;

		std::string name, prefix;
		sstream >> name >> prefix;

		experiment_names.push_back(std::make_tuple(std::move(name), std::move(prefix)));
	}

	infile.close();

	return std::make_tuple(replica_count, cutoff, std::move(experiment_names));
}

matrix_type calculate_covariance(const data_list& stat, const std::vector<data_list>& sys)
{
	const auto covariance = [&] (int i, int j) -> double {
		double err = sys[1][i] * sys[1][j];
		if (i == j)
			err += (stat[i]*stat[i]) + (sys[0][i]*sys[0][i]);
		return err;
	};

	matrix_type cov(stat.size(), data_list(stat.size(), 0));

	for (int i = 0; i < cov.size(); i++)
		for (int j = 0; j < cov[i].size(); j++)
			cov[i][j] = covariance(i, j);
	
	return cov;
}

std::tuple<data_list, data_list, matrix_type> read_experiment(const std::string& experiment, double cutoff)
{
	std::ifstream infile("../data/" + experiment + "/" + experiment + ".txt");
	
	std::string line;
	int ndata, nsys;
	double temp;

	infile >> line;			// name
	infile >> ndata >> nsys;	// ndata nsys
	infile >> temp;			// cme
	infile >> temp >> temp;		// ymin ymax
	std::getline(infile, line);	// trailing newline
	std::getline(infile, line);	// data header

	data_list pTs(ndata, 0);
	data_list experimental_data(ndata, 0);
	data_list stat_error(ndata, 0);
	std::vector<data_list> sys_error(nsys, data_list(ndata, 0));

	for (int i = 0; i < ndata; i++)
	{
		infile >> line;				// index
		infile >> pTs[i];			// pT
		
		infile >> experimental_data[i];		// obs
		infile >> stat_error[i];		// stat

		for (int j = 0; j < nsys; j++)
			infile >> sys_error[j][i];	// sys
	}

	infile.close();

	const auto pos = std::find_if(pTs.begin(), pTs.end(), within_cutoff(cutoff));
	const auto dist = std::distance(pTs.begin(), pos);

	const auto erase = [dist] (data_list& xs) { xs.erase(xs.begin(), xs.begin()+dist); };
	erase(pTs);
	erase(experimental_data);
	erase(stat_error);

	for (auto& err : sys_error)
		erase(err);

	matrix_type covariance = calculate_covariance(stat_error, sys_error);


	return std::make_tuple(std::move(pTs), std::move(experimental_data), std::move(covariance));
}

data_list read_theory(const std::string& experiment, const std::string& experiment_prefix, int replica, double cutoff)
{
	double conversion = 1;
	std::ifstream post_infile("../theory/" + experiment + "/post_process.txt");
	std::string line;
	while (std::getline(post_infile, line))
	{
		if (line.substr(0, 4) == "conv")
		{
			std::stringstream sstream;
			sstream << line.substr(5);
			sstream >> conversion;
		}
	}
	post_infile.close();


	std::stringstream formatted_replica;
	formatted_replica << std::setw(3) << std::setfill('0') << replica;
	std::ifstream infile("../theory/" + experiment + "/results/" + experiment_prefix + formatted_replica.str() + ".res");

	data_list pTs, theory_data;
	double pT;
	while (infile >> pT)
	{
		pTs.push_back(pT);

		double obs;
		infile >> obs;
		theory_data.push_back(obs * conversion);
	}

	infile.close();

	const auto pos = std::find_if(pTs.begin(), pTs.end(), within_cutoff(cutoff));
	const auto dist = std::distance(pTs.begin(), pos);

	const auto erase = [dist] (data_list& xs) { xs.erase(xs.begin(), xs.begin()+dist); };
	erase(pTs);
	erase(theory_data);

	return theory_data;
}

double fabs(const double d) {
	if (d < 0)
		return -d;
	else
		return d;
}

double calc_chi_component(int i, const data_list& exp_data, const data_list& thr_data) {
	if (thr_data[i] >= 0)
		return exp_data[i] - thr_data[i];
	return 0;
}

double chi_squared(const data_list& x_data, const data_list& exp_data, const data_list& thr_data, const matrix_type& cov)
{
	if (x_data.size() != exp_data.size() or x_data.size() != thr_data.size() or thr_data.size() != exp_data.size())
	{
		std::cout << "ERROR: INCORRECT DATASET SIZES:" << std::endl;
		std::cout << "THEORY: " << thr_data.size() << std::endl;
		std::cout << "X: " << x_data.size() << std::endl;
		std::cout << "EXP: " << exp_data.size() << std::endl;
	}

	double chi2 = 0;

	const auto inv_cov = cholesky_inverse(cov);

	//const int width = 10;
	//std::cout << std::setw(width) << "X" << std::setw(width+1) << "EXP" << std::setw(width+1) << "THR" << std::setw(width+1) << "DIFF" << std::endl << std::scientific << std::setprecision(2);
	//std::cout << std::setw(width) << "X" << std::setw(width+1) << "SLICE" << std::endl << std::scientific << std::setprecision(2);
	for (int i = 0; i < exp_data.size(); i++)
	{
		const auto i_slice = calc_chi_component(i, exp_data, thr_data);
		//std::cout << std::setw(width) << x_data[i] << " " << std::setw(width) << exp_data[i] << " " << std::setw(width) << thr_data[i] <<  " " << std::setw(width) << exp_data[i] - thr_data[i] << std::endl;
		//std::cout << std::setw(width) << x_data[i] << std::setw(width) << i_slice << std::endl;
		for (int j = 0; j < exp_data.size(); j++)
		{
			const auto j_slice = calc_chi_component(j, exp_data, thr_data);

			chi2 += i_slice * inv_cov[i][j] * j_slice;
		}
	}

	return chi2;
}

data_list calculate_reweights(const data_list& chi2s, const int dof)
{
	data_list reweights(chi2s.size(), 0);

	const double n = static_cast<double>(dof);
	const auto weight_func = [n] (double chi2) -> double {
		return pow(chi2, (n-1.)/2.) * exp(-0.5 * chi2);
	};
	std::transform(chi2s.begin(), chi2s.end(), reweights.begin(), weight_func);

	const double norm = std::accumulate(reweights.begin(), reweights.end(), 0.0) / static_cast<double>(reweights.size());
	std::transform(reweights.begin(), reweights.end(), reweights.begin(), [norm] (double x) { return x / norm; });

	return reweights;
}


#endif
