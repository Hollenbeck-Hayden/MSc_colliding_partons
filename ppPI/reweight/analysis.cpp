#include "utility.h"
#include <numeric>

data_list read_reweights(const std::string& filename)
{
	std::ifstream infile(filename);
	data_list reweights;
	double weight;
	while (infile >> weight)
		reweights.push_back(weight);
	infile.close();
	return reweights;
}

double probability_dist(const data_list& chi2s, const int n_points, const double alpha)
{
	const double n = static_cast<double>(n_points);
	const double alpha2 = alpha * alpha;
	return std::transform_reduce(chi2s.begin(), chi2s.end(), 0.0, std::plus<double>(),
		[n, alpha2] (double chi2) -> double {
			return pow(chi2/alpha2, (n-1.0)/2.0) * exp(-0.5 * chi2 / alpha2);
		}) / alpha;
}


int main(int argc, char** argv)
{
	if (argc < 2)
	{
		std::cerr << "Usage: " << std::string(argv[1]) << " paths..." << std::endl;
		std::exit(1);
	}

	double total_original_chi2 = 0.0;
	double total_weighted_chi2 = 0.0;
	int total_npoints = 0;

	for (int e = 1; e < argc; e++)
	{
		const std::string path(argv[e]);

		int n_replica = 0;
		double cutoff = 0;
		std::vector<std::tuple<std::string, std::string>> experiments;
		std::tie(n_replica, cutoff, experiments) = read_experiment_file(path + "/plot_experiments.txt");

		data_list reweights = read_reweights(path + "/reweights.data");
		double N_weights = std::accumulate(reweights.begin(), reweights.end(), 0.);
		std::cout << "N_weights: " << N_weights << std::endl;

		std::ofstream outfile(path + "/result_chi2s.txt");

		for (const auto& experiment : experiments)
		{
			const auto experiment_name = std::get<0>(experiment);
			const auto experiment_prefix = std::get<1>(experiment);

			data_list xs, experimental_data;
			matrix_type covariance;

			// Experimental data y
			std::tie(xs, experimental_data, covariance) = read_experiment(experiment_name, cutoff);

			// Ensemble average over replicas to get y[f]
			data_list theory_data(xs.size(), 0);
			data_list weight_data(xs.size(), 0);
			for (int i = 1; i <= n_replica; i++)
			{
				const auto theory_replica = read_theory(experiment_name, experiment_prefix, i, cutoff);
				for (int j = 0; j < theory_data.size(); j++)
				{
					theory_data[j] +=                theory_replica[j] / static_cast<double>(n_replica);
					weight_data[j] += reweights[i] * theory_replica[j] / N_weights;
				}
			}

			const double N_points = static_cast<double>(xs.size());
			total_npoints += xs.size();


			// Compute chi^2
			double original_chi2 = chi_squared(xs, experimental_data, theory_data, covariance);
			double weighted_chi2 = chi_squared(xs, experimental_data, weight_data, covariance);

			total_original_chi2 += original_chi2;
			total_weighted_chi2 += weighted_chi2;

			std::cout << "Experiment: " << experiment_name << std::endl;
			std::cout << "CUTOFF: " << cutoff << std::endl;
			std::cout << "N = " << N_points << std::endl;
			std::cout << "Original: " << original_chi2 << ", " << original_chi2 / N_points << std::endl;
			std::cout << "Weighted: " << weighted_chi2 << ", " << weighted_chi2 / N_points << std::endl;
			std::cout << "Ratio: " << weighted_chi2 * 100.0 / original_chi2 << std::endl;
			
			outfile << experiment_name << "\t" << static_cast<int>(N_points) << std::setprecision(2) << std::scientific
						   << "\t" << original_chi2 << "\t" << weighted_chi2
						   << std::fixed << std::setw(10) << std::right
						   << "\t\t" << std::setw(10) << original_chi2 / N_points << "\t" << std::setw(10) << weighted_chi2 / N_points
						   << "\t\t" << std::setw(10) << weighted_chi2 * 100.0 / original_chi2 << std::endl;
		}

		outfile << "total \t" << total_npoints << total_original_chi2 << "\t" << total_weighted_chi2
			<< "\t\t" << total_original_chi2 / static_cast<double>(total_npoints) << "\t" << total_weighted_chi2 / static_cast<double>(total_npoints)
			<< "\t\t" << total_weighted_chi2 * 100.0 / total_original_chi2 << std::endl;


		outfile.close();


		std::ofstream chi_o_file(path + "/replica_chi2_original.txt");

		int n_points = 0;
		data_list replica_chi2s(n_replica, 0);
		for (int i = 1; i <= n_replica; i++)
		{
			double replica_chi2_original = 0;

			for (const auto& experiment : experiments)
			{
				const auto experiment_name = std::get<0>(experiment);
				const auto experiment_prefix = std::get<1>(experiment);

				data_list xs, experimental_data;
				matrix_type covariance;
				std::tie(xs, experimental_data, covariance) = read_experiment(experiment_name, cutoff);
				const auto theory_data = read_theory(experiment_name, experiment_prefix, i, cutoff);
				
				if (i == 1) n_points += xs.size();
				replica_chi2_original += chi_squared(xs, experimental_data, theory_data, covariance);
			}

			replica_chi2s[i-1] = replica_chi2_original;
			chi_o_file << replica_chi2_original / static_cast<double>(n_points) << std::endl;
		}

		chi_o_file.close();

		data_list alphas, probs;
		for (double alpha = 0.4; alpha < 5; alpha += 0.01)
		{
			alphas.push_back(alpha);
			probs.push_back(probability_dist(replica_chi2s, n_points, alpha));
		}

		const auto prob_total = std::accumulate(probs.begin(), probs.end(), 0.0, std::plus<double>());
		std::transform(probs.begin(), probs.end(), probs.begin(), [prob_total] (double p) { return p / prob_total; });

		std::ofstream prob_file(path + "/probability_density.txt");

		for (int i = 0; i < alphas.size(); i++)
		{
			prob_file << alphas[i] << "\t" << probs[i] << std::endl;
		}

		prob_file.close();
	}

	return 0;
}
