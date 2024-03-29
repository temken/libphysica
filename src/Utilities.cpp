#include "libphysica/Utilities.hpp"

#include <algorithm>
#include <cmath>
#include <sys/stat.h>	 //required to create a folder
#include <sys/types.h>	 // required for stat.h

#include "libphysica/List_Manipulations.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "version.hpp"

namespace libphysica
{
using namespace libphysica::natural_units;
using namespace libconfig;

// 1. Terminal output
std::string Time_Display(double seconds_total)
{
	std::vector<std::string> units_strings = {"y", "w", "d", "h", "m", "s", "ms"};
	std::vector<double> units			   = {year, week, day, hr, minute, sec, milli * sec};
	std::vector<int> times;
	std::vector<std::string> time_strings;
	for(unsigned int i = 0; i < units.size(); i++)
	{
		times.push_back(std::floor(seconds_total * sec / units[i]));
		seconds_total -= 1.0 * times.back() * units[i] / sec;
		time_strings.push_back(std::to_string(times.back()));
		if(time_strings.back().length() == 1)
			time_strings.back() = "0" + time_strings.back();
		if(time_strings.back().length() == 2 && i == units.size() - 1)
			time_strings.back() = "0" + time_strings.back();
	}
	unsigned int i;
	for(i = 0; i < 4; i++)
		if(times[i] > 0)
			break;
	std::string separator = ":";
	return "[" + time_strings[i] + units_strings[i] + separator + time_strings[i + 1] + units_strings[i + 1] + separator + time_strings[i + 2] + units_strings[i + 2] + "]";
}

void Print_Progress_Bar(double progress, unsigned int MPI_rank, unsigned int bar_length, double time, std::string bar_color)
{
	if(MPI_rank == 0 && progress >= 0.0 && progress <= 1.0)
	{
		std::cout << "\r";
		for(unsigned int j = 0; j < 2 * bar_length; j++)
			std::cout << " ";
		std::cout << "\r";
		for(unsigned int i = 0; i < bar_length; i++)
		{
			if(i == bar_length / 2)
			{
				if(progress < 0.1 || progress == 1.0)
				{
					std::cout << Round(100.0 * progress, 1) << "%" << std::flush;
					i += (progress < 0.01 || progress == 1.0) ? 3 : 1;
				}
				else
				{
					std::cout << Round(100.0 * progress, 2) << "%" << std::flush;
					i += 2;
				}
			}
			else if(progress > 1.0 * i / bar_length)
				std::cout << Formatted_String("█", bar_color) << std::flush;
			else
				std::cout << Formatted_String("░", bar_color) << std::flush;
		}
		if(time > 0.0 && progress > 1.0e-3)
		{
			double t = (progress > 0.9999) ? time : (1.0 - progress) * time / progress;
			std::cout << " " << Time_Display(t) << std::flush;
		}
	}
}

void Print_Box(std::string str, unsigned int tabs, int mpi_rank, std::string box_color, std::string text_color)
{
	if(mpi_rank == 0)
	{
		std::string box_string_1, box_string_2 = "";
		unsigned int length = str.length() + 2;
		for(unsigned int i = 0; i < tabs; i++)
			box_string_1 += "\t";
		box_string_1 += "╔";
		for(unsigned int i = 0; i < length; i++)
			box_string_1 += "═";
		box_string_1 += "╗\n";
		for(unsigned int i = 0; i < tabs; i++)
			box_string_1 += "\t";
		box_string_1 += "║ ";
		box_string_2 += " ║\n";
		for(unsigned int i = 0; i < tabs; i++)
			box_string_2 += "\t";
		box_string_2 += "╚";
		for(unsigned int i = 0; i < length; i++)
			box_string_2 += "═";
		box_string_2 += "╝\n";

		std::cout << Formatted_String(box_string_1, box_color, true) << Formatted_String(str, text_color, true) << Formatted_String(box_string_2, box_color, true) << std::endl;
	}
}

const std::vector<std::string> colors				  = {"Default", "Black", "Red", "Green", "Yellow", "Blue", "Magenta", "Cyan", "White"};
const std::vector<std::string> color_codes			  = {"0", "30", "31", "32", "33", "34", "35", "36", "37"};
const std::vector<std::string> background_color_codes = {"49", "40", "41", "42", "43", "44", "45", "46", "47"};
extern std::string Formatted_String(std::string str, std::string color, bool bold, bool underlined, std::string background_color)
{
	if(color == "Default" && !bold)
		return str;
	else if(List_Contains(colors, color) == false || List_Contains(colors, background_color) == false)
	{
		std::cerr << Formatted_String("Warning", "Yellow", true) << ": in libphysica::Formatted_String(): Unknown color " << color << " or background color " << background_color << "." << std::endl;
		return str;
	}
	else
	{
		int i_color						  = Find_Indices(colors, color)[0];
		int i_background_color			  = Find_Indices(colors, background_color)[0];
		std::string color_code			  = color_codes[i_color];
		std::string background_color_code = background_color_codes[i_background_color];
		std::string bold_code			  = (bold) ? "1" : "0";
		std::string underlined_code		  = (underlined) ? "4" : "0";
		return "\033[" + bold_code + ";" + underlined_code + ";" + color_code + ";" + background_color_code + "m" + str + "\033[0m";
	}
}

// 2. Import and export data from files
bool File_Exists(const std::string& file_path)
{
	struct stat buffer;
	return (stat(file_path.c_str(), &buffer) == 0);
}

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
		std::cerr << "Error in libphysica::Import_Data(" << filepath << "): File does not exist." << std::endl;
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
			std::cerr << "Error in libphysica::Import_Data(): Column length and dimension length do not match." << std::endl;
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
		std::cerr << "Error in libphysica::Import_Data(" << filepath << "): File does not exist." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void Create_Folder(const std::string& path, int mpi_rank, bool terminal_output)
{
	if(mpi_rank == 0)
	{
		mode_t nMode = 0733;   // UNIX style permissions
		int nError	 = 0;
#if defined(_WIN32)
		nError = _mkdir(path.c_str());	 // can be used on Windows
#else
		nError = mkdir(path.c_str(), nMode);   // can be used on non-Windows
#endif
		if(nError != 0 && terminal_output)
			std::cerr << "\nWarning in libphysica::Create_Folder(): The folder " << path << " exists already." << std::endl
					  << std::endl;
	}
}

void Export_List(std::string filepath, std::vector<double> data, double dimension, const std::string& header)
{
	std::ofstream outputfile;
	outputfile.open(filepath);
	if(header.length() > 0)
		outputfile << header << std::endl;
	for(unsigned int i = 0; i < data.size(); i++)
		outputfile << In_Units(data[i], dimension) << std::endl;
	outputfile.close();
}

void Export_Table(std::string filepath, const std::vector<std::vector<double>>& data, std::vector<double> dimensions, const std::string& header)
{
	std::ofstream outputfile;
	outputfile.open(filepath);
	if(header.length() > 0)
		outputfile << header << std::endl;
	unsigned int lines = data.size();
	for(unsigned int line = 0; line < lines; line++)
	{
		unsigned int columns = data[line].size();
		if(!dimensions.empty() && dimensions.size() != columns)
		{
			std::cerr << "Error in libphysica::Export_Data(): Column length and dimension length do not match." << std::endl;
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

void Export_Function(std::string filepath, std::function<double(double)> func, const std::vector<double>& x_list, std::vector<double> dimensions, const std::string& header)
{
	std::vector<std::vector<double>> data;
	for(auto& x : x_list)
		data.push_back({x, func(x)});
	Export_Table(filepath, data, dimensions, header);
}

void Export_Function(std::string filepath, std::function<double(double)> func, double xMin, double xMax, unsigned int steps, std::vector<double> dimensions, bool logarithmic, const std::string& header)
{
	std::vector<double> x_list = (logarithmic) ? Log_Space(xMin, xMax, steps) : Linear_Space(xMin, xMax, steps);
	Export_Function(filepath, func, x_list, dimensions, header);
}

// 3. Create lists with equi-distant numbers
std::vector<int> Range(int max)
{
	return Range(0, max, 1);
}

std::vector<int> Range(int min, int max, int stepsize)
{
	std::vector<int> range;
	if(min > max && stepsize > 0)
		for(int i = min; i > max; i -= stepsize)
			range.push_back(i);
	else
		for(int i = min; i < max; i += stepsize)
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

// 5. Configuration class using libconfig
void Configuration::Read_Config_File()
{
	try
	{
		config.readFile(cfg_file.c_str());
	}
	catch(const FileIOException& fioex)
	{
		std::cerr << "Error in libphysica::Configuration::Read_Config_File(): I/O error while reading configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	catch(const ParseException& pex)
	{
		std::cerr << "Error in libphysica::Configuration::Read_Config_File(): Configurate file parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void Configuration::Initialize_Result_Folder(int MPI_rank)
{
	try
	{
		ID = config.lookup("ID").c_str();
	}
	catch(const SettingNotFoundException& nfex)
	{
		std::cerr << "Error in libphysica::Configuration::Initialize_Result_Folder(): No 'ID' setting in configuration file." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	// 1. Create the /results/ folder if necessary
	std::string results_folder = TOP_LEVEL_DIR "results/";
	Create_Folder(results_folder, MPI_rank, false);
	// 2. Create a /result/<ID>/ folder for result files.
	results_path = results_folder + ID + "/";
	Create_Folder(results_path, MPI_rank);
	Copy_Config_File(MPI_rank);
}

void Configuration::Copy_Config_File(int MPI_rank)
{
	if(MPI_rank == 0)
	{
		std::ifstream inFile;
		std::ofstream outFile;
		inFile.open(cfg_file);
		outFile.open(TOP_LEVEL_DIR "results/" + ID + "/" + ID + ".cfg");
		outFile << "// " << TOP_LEVEL_PROJECT_NAME << "-v" << TOP_LEVEL_PROJECT_VERSION << "\tgit:" << TOP_LEVEL_GIT_BRANCH << "/" << TOP_LEVEL_GIT_COMMIT_HASH << std::endl;
		outFile << inFile.rdbuf();
		inFile.close();
		outFile.close();
	}
}

Configuration::Configuration(std::string cfg_filename, int MPI_rank)
: cfg_file(cfg_filename), results_path("./")
{

	// 1. Read the cfg file.
	Read_Config_File();

	// 2. Find the run ID, create a folder and copy the cfg file.
	Initialize_Result_Folder(MPI_rank);
}

// 6. Other utilities
std::vector<int> Workload_Distribution(unsigned int workers, unsigned int tasks)
{
	int tasks_per_worker = tasks / workers;
	std::vector<int> index_list(workers + 1, 0);
	for(unsigned int i = 0; i < workers; i++)
		index_list[i + 1] = index_list[i] + tasks_per_worker;
	// Distribute the remainder on the workers starting at the end of the list.
	int remainder = tasks % workers;
	for(int i = 0; i < remainder; i++)
		index_list[workers - i] += (remainder - i);
	return index_list;
}

unsigned int Locate_Closest_Location(const std::vector<double>& sorted_list, double target)
{
	if(std::is_sorted(std::begin(sorted_list), std::end(sorted_list)) == false)
	{
		std::cerr << "Error in libphysica::Locate_Closest_Location(): The list is not sorted." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	auto const it = std::upper_bound(sorted_list.begin(), sorted_list.end(), target);
	if(it == sorted_list.end())
		return sorted_list.size() - 1;
	else
	{
		unsigned int index = std::distance(sorted_list.begin(), it);
		if(index == 0)
			return 0;
		else if(index == sorted_list.size())
			return sorted_list.size() - 1;
		else
		{
			double diff1 = std::fabs(sorted_list[index - 1] - target);
			double diff2 = std::fabs(sorted_list[index] - target);
			if(diff1 < diff2)
				return index - 1;
			else
				return index;
		}
	}
}

void Check_For_Error(bool error_condition, std::string function_name, std::string error_message)
{
	if(error_condition)
	{
		std::cerr << Formatted_String("Error", "Red", true) << " in " << function_name << ": " << error_message << std::endl;
		std::exit(EXIT_FAILURE);
	}
}
void Check_For_Warning(bool warning_condition, std::string function_name, std::string warning_message)
{
	if(warning_condition)
		std::cerr << Formatted_String("Warning", "Yellow", true) << " in " << function_name << ": " << warning_message << std::endl;
}

}	// namespace libphysica
