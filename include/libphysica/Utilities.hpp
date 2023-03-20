#ifndef __Utilities_hpp_
#define __Utilities_hpp_

#include <fstream>
#include <functional>
#include <iostream>
#include <libconfig.h++>
#include <string>
#include <vector>

namespace libphysica
{

// 1. Terminal output
extern std::string Time_Display(double seconds);
extern void Print_Progress_Bar(double progress, unsigned int MPI_rang = 0, unsigned int bar_length = 50, double time = 0.0, std::string bar_color = "Default");
extern void Print_Box(std::string str, unsigned int tabs = 0, int mpi_rank = 0, std::string box_color = "Default", std::string text_color = "Default");
extern std::string Colored_Text(std::string str, std::string color = "Red", bool bold = true, std::string background_color = "Default");

// 2. Import and export data from files
extern bool File_Exists(const std::string& file_path);

extern std::vector<double> Import_List(std::string filepath, double dimension = 1.0, unsigned int ignored_initial_lines = 0);
extern std::vector<std::vector<double>> Import_Table(std::string filepath, std::vector<double> dimensions = {}, unsigned int ignored_initial_lines = 0);

extern void Create_Folder(const std::string& path, int mpi_rank = 0, bool terminal_output = true);
extern void Export_List(std::string filepath, std::vector<double> data, double dimension = 1.0, const std::string& header = "");
extern void Export_Table(std::string filepath, const std::vector<std::vector<double>>& data, std::vector<double> dimensions = {}, const std::string& header = "");
extern void Export_Function(std::string filepath, std::function<double(double)> func, const std::vector<double>& x_list, std::vector<double> dimensions = {}, const std::string& header = "");
extern void Export_Function(std::string filepath, std::function<double(double)> func, double xMin, double xMax, unsigned int steps, std::vector<double> dimensions = {}, bool logarithmic = false, const std::string& header = "");

// 3. Create lists with equi-distant numbers
extern std::vector<int> Range(int max);
extern std::vector<int> Range(int min, int max, int stepsize = 1);
extern std::vector<double> Linear_Space(double min, double max, unsigned int steps);
extern std::vector<double> Log_Space(double min, double max, unsigned int steps);

// 4. Dual stream class to write onto terminal and a log file simultaneously.
class Logger
{
  private:
	std::ofstream logfile_stream;

  public:
	Logger(const std::string& file_name)
	{
		logfile_stream.open(file_name);
	};

	template <class T>
	Logger& operator<<(const T& x)
	{
		std::cout << x;
		logfile_stream << x;
		return *this;
	}
	Logger& operator<<(std::ostream& (*os)(std::ostream&) )
	{
		std::cout << os;
		logfile_stream << os;
		return *this;
	}
};

// 5. Configuration class using libconfig
class Configuration
{
  protected:
	libconfig::Config config;
	std::string cfg_file;

	virtual void Read_Config_File();

	virtual void Initialize_Result_Folder(int MPI_rank = 0);
	virtual void Copy_Config_File(int MPI_rank = 0);

	virtual void Initialize_Parameters() {};

  public:
	std::string ID, results_path;
	explicit Configuration(std::string cfg_filename, int MPI_rank = 0);

	virtual void Print_Summary(int MPI_rank = 0) {};
};

// 6. Other utilities
extern std::vector<int> Workload_Distribution(unsigned int workers, unsigned int tasks);
extern unsigned int Locate_Closest_Location(const std::vector<double>& sorted_list, double target);
extern void Check_For_Error(bool error_condition, std::string function_name, std::string error_message);
extern void Check_For_Warning(bool warning_condition, std::string function_name, std::string warning_message);

}	// namespace libphysica

#endif