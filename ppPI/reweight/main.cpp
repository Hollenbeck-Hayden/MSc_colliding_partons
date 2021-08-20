#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <tuple>
#include <vector>
#include <iomanip>
#include <sstream>

#include "matrix.h"
#include "unweight.h"
#include "utility.h"

#include <algorithm>
#include <numeric>

std::tuple<double, int> calculate_experiment_data(const std::vector<std::tuple<std::string, std::string>>& experiment_names, const int replica, const double cutoff)
{
	double chi2 = 0;
	int n = 0;

	for (const auto& experiment : experiment_names)
	{
		data_list pTs, experimental_data, theoretical_data;
		matrix_type covariance;
		
		std::tie(pTs, experimental_data, covariance) = read_experiment(std::get<0>(experiment), cutoff);
		theoretical_data = read_theory(std::get<0>(experiment), std::get<1>(experiment), replica, cutoff);

		n += pTs.size();
		chi2 += chi_squared(pTs, experimental_data, theoretical_data, covariance);
	}

	return std::make_tuple(chi2, n);
}

const double min_eps = 1e-6;
unsigned int calculate_Neff(const data_list& reweights)
{
	const double N = static_cast<double>(reweights.size());
	return static_cast<unsigned int>(exp(std::transform_reduce(
		reweights.begin(), reweights.end(), 0., std::plus<>(), [N] (double w) -> double {
			return w > 1e-10 ? w * log(N / w) : 0.0; }
		) / N));
}

unsigned int theta(double x)
{
	if (x < 0)
		return 0;
	else
		return 1;
}

std::vector<unsigned int> calculate_unweights(const data_list& reweights)
{
	// Calculate cummulative probabilities
	const double Nrep = static_cast<double>(reweights.size());
	data_list Pks;
	for (int i = 0; i < reweights.size(); i++)
	{
		if (i == 0)
			Pks.push_back(reweights[i] / Nrep);
		else
			Pks.push_back(Pks[i-1] + (reweights[i] / Nrep));
	}
	Pks.back() = 1.1;

	// Calculate Effective N'
	const unsigned int Neff = calculate_Neff(reweights);
	std::cout << "Neff: " << Neff << std::endl;

	// Unweight
	std::vector<unsigned int> unweights(Pks.size(), 0);
	for (int k = 0; k < Pks.size(); k++)
	{
		const double Pk_left  = (k == 0) ? 0 : Pks[k-1];
		const double Pk_right = Pks[k];

		for (int j = 0; j < Neff; j++)
		{
			const double dj = static_cast<double>(j) / static_cast<double>(Neff);
			unweights[k] += theta(dj - Pk_left) * theta(Pk_right - dj);
		}
	}

	return unweights;
}

template<typename T>
void write_to_file(const std::string& filename, const std::vector<T>& data)
{
	std::ofstream outfile(filename);
	for (const auto& d : data)
		outfile << d << std::endl;
	outfile.close();
}

void reweight_experimental_file(const std::string& path) {
	int replica_count = 0;
	double cutoff = 0;
	std::vector<std::tuple<std::string, std::string>> experiment_names;

	std::tie(replica_count, cutoff, experiment_names) = read_experiment_file(path + "/experiments.txt");
	std::cout << "Num replicas: " << replica_count << std::endl;

	data_list chi2s(replica_count, 0.);
	int n = 0;

	for (int k = 0; k < replica_count; k++)
		std::tie(chi2s[k], n) = calculate_experiment_data(experiment_names, k, cutoff);

	const auto reweights = calculate_reweights(chi2s, n);
	std::cout << "Reweights:" << std::scientific << std::endl;
	for (int k = 0; k < replica_count; k++)
		std::cout << k+1 << "\t" << reweights[k] << "\t" << chi2s[k] << std::endl;

	std::cout << "Sum: " << std::accumulate(reweights.begin(), reweights.end(), 0.) << std::endl;

	const auto unweights = calculate_unweights(reweights);
	std::cout << "Unweights:" << std::scientific << std::endl;
	for (int k = 0; k < replica_count; k++)
		std::cout << k+1 << "\t" << unweights[k] << std::endl;

	unsigned int total = 0;
	for (auto w : unweights)
		total += w;
	std::cout << "Sum: " << total << std::endl;

	write_to_file(path + "/reweights.data", reweights);
	write_to_file(path + "/unweights.data", unweights);

	//saveWeightedSet("MAPFF10NLOPIp", "testff_p", unweights, "./testff_p/");
	//saveWeightedSet("MAPFF10NLOPIm", "testff_m", unweights, "./testff_m/");
}

int main(int argc, char** argv)
{
	if (argc < 2) {
		std::cerr << "Usage: " << std::string(argv[0]) << " paths..." << std::endl;
		std::cerr << "paths is a list of experiment directory paths containing a experiments.txt file to reweight" << std::endl;

		std::exit(1);
	}

	for (int i = 1; i < argc; i++) {
		const std::string path = argv[i];
		reweight_experimental_file(path);
	}

	return 0;
}
