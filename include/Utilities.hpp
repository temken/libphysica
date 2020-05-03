#ifndef __Utilities_hpp_
#define __Utilities_hpp_

#include <vector>
#include <string>
#include <functional>

//1. Progress bar
	extern void Print_Progress_Bar(double progress, int MPI_rang = 0);
	extern void Print_Progress_Bar(double i, double iMax, int MPI_rang = 0);

//2. Import and export data from files
	extern std::vector<double> Import_List(std::string filepath,double dimension=1.0);
	extern std::vector<std::vector<double>> Import_Table(std::string filepath, std::vector<double> dimensions = {});

	extern void Export_List(std::string filepath, std::vector<double> data, double dimension = 1.0);
	extern void Export_Table(std::string filepath, const std::vector<std::vector<double>>& data, std::vector<double> dimensions = {});
	extern void Export_Function(std::string filepath, std::function<double(double)>& func, double xMin, double xMax, unsigned int steps, std::vector<double> dimensions = {}, bool logarithmic = false);

//3. Create list with equi-distant numbers in log-space
	extern std::vector<double> Linear_Space(double min, double max, unsigned int steps);
	extern std::vector<double> Log_Space(double min, double max, unsigned int steps);


#endif