#ifndef UNWEIGHT_H
#define UNWEIGHT_H

#include "yaml-cpp/yaml.h"

const std::string LHAPDF_PATH = "/usr/local/share/LHAPDF/";

void saveWeightedSet(const std::string& pdf_name,
		     const std::string& new_name,
		     const std::vector<unsigned int>& unweights,
		     const std::string& path)
{
	std::vector<std::vector<double>> x_grid, q_grid;
	std::vector<std::vector<int>> flavors;

	// Assume that all PDFs have the same x-grid, q-grid, and subgrid structure
	std::ifstream central_infile(LHAPDF_PATH + pdf_name + "/" + pdf_name + "_0000.dat");
	
	// skip over header
	for (std::string line = ""; std::getline(central_infile, line); )
	{
		// Read subgrid
		if (line == "---")
		{
			if (std::getline(central_infile, line))
			{
				// Read x grid
				std::stringstream x_line(line);
				x_grid.push_back(std::vector<double>());
				
				double x;
				while (x_line >> x)
					x_grid.back().push_back(x);

				// Read q grid
				std::getline(central_infile, line);
				std::stringstream q_line(line);
				q_grid.push_back(std::vector<double>());
				
				double q;
				while (q_line >> q)
					q_grid.back().push_back(q);

				// Read flavors
				std::getline(central_infile, line);
				std::stringstream f_line(line);
				flavors.push_back(std::vector<int>());

				int f;
				while (f_line >> f)
					flavors.back().push_back(f);
			}
			else
			{
				break;
			}
		}
	}
	central_infile.close();

	const int numGrids = x_grid.size();


	// Copy replica files based on unweighting
	int index = 1;
	for (int replica = 0; replica < unweights.size(); replica++)
	{
		std::stringstream old_path;
		old_path << LHAPDF_PATH << pdf_name << "/" << pdf_name << "_" << std::setfill('0') << std::setw(4) << replica+1 << ".dat";

		// Copy replica unweights[replica] times
		for (unsigned int n = 0; n < unweights[replica]; n++)
		{
			std::stringstream new_path;
			new_path << path << new_name << "_" << std::setfill('0') << std::setw(4) << index << ".dat";

			std::ifstream  infile(old_path.str(), std::ios::binary);
			std::ofstream outfile(new_path.str(), std::ios::binary);

			outfile << infile.rdbuf();

			infile.close();
			outfile.close();
			
			index++;
		}
	}

	const int COUNT = index;

	// Create central replica 0 file
	std::vector<std::vector<std::vector<double>>> total;
	for (int i = 0; i < numGrids; i++)
	{
		total.push_back(std::vector<std::vector<double>>());
		for (int j = 0; j < flavors[i].size(); j++)
		{
			total[i].push_back(std::vector<double>(x_grid[i].size() * q_grid[i].size(), 0));
		}
	}

	const int READ_HEADER = 0;
	const int READ_DATA = 1;

	for (int i = 1; i < index; i++)
	{
		std::stringstream replica_path;
		replica_path << path << new_name << "_" << std::setfill('0') << std::setw(4) << i << ".dat";

		std::ifstream infile(replica_path.str());
		std::string line;
		int state = READ_HEADER;
		int grid = -1;
		while (std::getline(infile, line))
		{
			if (line == "---")
			{
				state = READ_DATA;
				grid++;
				std::getline(infile, line);	// Skip x grid
				std::getline(infile, line);	// Skip q grid
				std::getline(infile, line);	// Skip flavors
				index = 0;
			}
			else if (state == READ_DATA)
			{
				std::stringstream data_line;
				data_line << line;
				double value;
				for (int f = 0; data_line >> value; f++)
				{
					total[grid][f][index] += value;
				}
				index++;
			}
		}
	}

	for (int i = 0; i < total.size(); i++)
		for (int j = 0; j < total[i].size(); j++)
			for (int k = 0; k < total[i][j].size(); k++)
				total[i][j][k] /= (double) COUNT;

	std::stringstream central_path;
	central_path << path << new_name << "_0000.dat";
	std::ofstream central_outfile(central_path.str());
	
	central_outfile << "PdfType: replica" << std::endl
			<< "Format: lhagrid1" << std::endl
			<< "---" << std::endl;

	for (int i = 0; i < total.size(); i++)
	{
		for (int j = 0; j < x_grid[i].size(); j++)
			central_outfile << std::scientific << std::setprecision(6) << x_grid[i][j] << " ";
		central_outfile << std::endl;

		for (int j = 0; j < q_grid[i].size(); j++)
			central_outfile << std::scientific << std::setprecision(6) << q_grid[i][j] << " ";
		central_outfile << std::endl;

		for (int j = 0; j < flavors[i].size(); j++)
			central_outfile << flavors[i][j] << " ";
		central_outfile << std::endl;

		for (int j = 0; j < total[i][0].size(); j++)
		{
			for (int f = 0; f < flavors[i].size(); f++)
				central_outfile << std::scientific << std::setprecision(6) << total[i][f][j] << " ";
			central_outfile << std::endl;
		}

		central_outfile << "---" << std::endl;
	}

	central_outfile.close();

	// Create .info file
	auto info_node = YAML::LoadFile(LHAPDF_PATH + pdf_name + "/" + pdf_name + ".info");
	info_node["SetDesc"] = new_name + " reweighted";
	info_node["NumMembers"] = COUNT;

	std::ofstream info_outfile(path + new_name + ".info");
	info_outfile << info_node;
	info_outfile << std::endl;
	info_outfile.close();
}

#endif
