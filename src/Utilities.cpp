#include "Utilities.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

#include "Natural_Units.hpp"
#include "Numerics.hpp"

namespace libphysica
{
using namespace libphysica::natural_units;

//1. Progress bar
void Print_Progress_Bar(double progress, int MPI_rank)
{
	if(MPI_rank == 0)
	{
		int BarLength = 50;
		std::cout << "\r";
		for(int j = 0; j < 3 * BarLength; j++)
			std::cout << " ";
		std::cout << "\r";
		if(progress >= 0.0 && progress < 0.999)
		{
			for(int i = 0; i < BarLength; i++)
			{
				if(i == BarLength / 2)
				{
					if(progress < 0.1)
					{
						std::cout << Round(100.0 * progress, 1) << "%" << std::flush;
						i += (progress < 0.01) ? 3 : 1;
					}
					else
					{
						std::cout << Round(100.0 * progress, 2) << "%" << std::flush;
						i += 2;
					}
				}
				else if(progress > 1.0 * i / BarLength)
					std::cout << "=" << std::flush;
				else
					std::cout << " " << std::flush;
			}
			std::cout << "|";	//<<std::flush;
		}
	}
}

void Print_Progress_Bar(double i, double iMax, int MPI_rank)
{
	Print_Progress_Bar(1.0 * i / iMax, MPI_rank);
}

//2. Import and export data from files
std::vector<double> Import_List(std::string filepath, double dimension, unsigned int ignored_initial_lines)
{
	std::vector<double> data;
	std::ifstream inputfile;
	inputfile.open(filepath);
	if(inputfile.good())
	{
		for(unsigned int i = 0; i < ignored_initial_lines; i++)
			inputfile.ignore(10000, '\n');
		double x;
		while(inputfile >> x)
			data.push_back(x * dimension);
		inputfile.close();
	}
	else
	{
		std::cerr << "Error in Import_Data(" << filepath << "): File does not exist." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	return data;
}

unsigned int Count_Lines(std::string filepath)
{
	unsigned int line_count = 0;
	std::ifstream file(filepath);
	if(file.is_open())
	{
		std::string line;
		while(std::getline(file, line))
			line_count++;
		file.close();
	}
	return line_count;
}

std::vector<std::vector<double>> Import_Table(std::string filepath, std::vector<double> dimensions, unsigned int ignored_initial_lines)
{
	std::vector<double> data_aux = {};
	std::ifstream inputfile;
	inputfile.open(filepath);
	if(inputfile.good())
	{
		for(unsigned int i = 0; i < ignored_initial_lines; i++)
			inputfile.ignore(10000, '\n');
		double x;
		while(inputfile >> x)
			data_aux.push_back(x);
		inputfile.close();

		unsigned int rows	 = Count_Lines(filepath) - ignored_initial_lines;
		unsigned int columns = data_aux.size() / rows;
		if(!dimensions.empty() && dimensions.size() != columns)
		{
			std::cerr << "Error in Import_Data(): Column length and dimension length do not match." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		std::vector<std::vector<double>> data(rows, std::vector<double>(columns, 0.0));
		unsigned int k = 0;
		for(unsigned int i = 0; i < rows; i++)
		{
			for(unsigned int j = 0; j < columns; j++)
			{
				double dim = dimensions.empty() ? 1.0 : dimensions[j];
				data[i][j] = data_aux[k] * dim;
				k++;
			}
		}
		return data;
	}
	else
	{
		std::cerr << "Error in Import_Data(" << filepath << "): File does not exist." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void Export_List(std::string filepath, std::vector<double> data, double dimension)
{
	std::ofstream outputfile;
	outputfile.open(filepath);
	for(unsigned int i = 0; i < data.size(); i++)
		outputfile << In_Units(data[i], dimension) << std::endl;
	outputfile.close();
}

void Export_Table(std::string filepath, const std::vector<std::vector<double>>& data, std::vector<double> dimensions)
{
	std::ofstream outputfile;
	outputfile.open(filepath);
	unsigned int lines = data.size();
	for(unsigned int line = 0; line < lines; line++)
	{
		unsigned int columns = data[line].size();
		if(!dimensions.empty() && dimensions.size() != columns)
		{
			std::cerr << "Error in Export_Data(): Column length and dimension length do not match." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		for(unsigned int column = 0; column < columns; column++)
		{
			double dim = dimensions.empty() ? 1.0 : dimensions[column];
			outputfile << In_Units(data[line][column], dim);
			if(column != columns - 1)
				outputfile << "\t";
			else if(line != lines - 1)
				outputfile << std::endl;
		}
	}
	outputfile.close();
}

void Export_Function(std::string filepath, std::function<double(double)> func, double xMin, double xMax, unsigned int steps, std::vector<double> dimensions, bool logarithmic)
{
	std::vector<std::vector<double>> data(steps, std::vector<double>(2, 0.0));
	std::vector<double> arguments = (logarithmic) ? Log_Space(xMin, xMax, steps) : Linear_Space(xMin, xMax, steps);
	for(unsigned int i = 0; i < arguments.size(); i++)
	{
		data[i][0] = arguments[i];
		data[i][1] = func(arguments[i]);
	}
	Export_Table(filepath, data, dimensions);
}

//3. Create lists with equi-distant numbers
std::vector<unsigned int> Range(unsigned int max)
{
	std::vector<unsigned int> range;
	for(unsigned int i = 0; i < max; i++)
		range.push_back(i);
	return range;
}

std::vector<int> Range(int min, int max, int stepsize)
{
	std::vector<int> range;
	if(min > max && stepsize > 0)
		for(unsigned int i = min; i > max; i -= stepsize)
			range.push_back(i);
	else
		for(unsigned int i = min; i < max; i += stepsize)
			range.push_back(i);
	return range;
}

std::vector<double> Linear_Space(double min, double max, unsigned int steps)
{
	if(steps < 2 || min == max)
		return {min};
	else
	{
		std::vector<double> result;
		double step = (max - min) / (steps - 1.0);

		for(unsigned int i = 0; i < steps; i++)
			result.push_back(min + i * step);
		return result;
	}
}

std::vector<double> Log_Space(double min, double max, unsigned int steps)
{
	if(steps < 2 || min == max)
		return {min};
	else
	{
		std::vector<double> result;
		double logmin = log(min);
		double dlog	  = log(max / min) / (steps - 1.0);

		for(unsigned int i = 0; i < steps; i++)
			result.push_back(exp(logmin + i * dlog));
		return result;
	}
}

}	// namespace libphysica
