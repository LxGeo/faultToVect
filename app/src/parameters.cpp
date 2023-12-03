#include "parameters.h"

namespace LxGeo
{
	namespace faultToVector
	{
		Parameters::Parameters(int argc, char* argv[])
		{
			init();
			parse(argc, argv);
			
		}


		Parameters::~Parameters()
		{
		}


		bool Parameters::initialized()
		{
			bool is_initialized = input_raster_path.empty() | input_thin_raster_path.empty();
			if (!is_initialized) help();

			return is_initialized;
		}


		void Parameters::init()
		{
			printed_help = false;

			input_raster_path.clear();
			input_thin_raster_path.clear();

			output_basename = "result";
			output_shapefile = "result.shp";
			temp_dir = "temp_dir";
			overwrite_output = false;
			min_heat_value = 0;
			ignore_thin = false;
			ignore_multi_thresh_thin = false;
			max_n_level_neighbour = 2;
			simplification_factor = 1.5f;

		}


		void Parameters::help()
		{
			if (printed_help) return;

			std::cout << "faultToVector.exe [args]" << std::endl
				<< std::endl
				<< "where [args] are : " << std::endl
				<< "  [-h | --help] -> print this help message" << std::endl
				<< std::endl
				<< "** Basic parameters : " << std::endl
				<< std::endl
				<< "  [-i] [input_raster_path] -> provide path of input prediction raster" << std::endl
				<< "  [-o] [basename] -> specify basename of output file" << std::endl
				<< "  [-max_n_level_neighbour] [value] -> specify maximum N level for neighbourhood creation" << std::endl
				<< "  [-min_heat_value] [-mhv] [value] -> specify minimum pixel value to include in vectorization (not included!). Default: 0" << std::endl
				<< "  [-simplification_factor] [-sf] [value] -> specify maximum allowed error for pixel vectorization. Default: 1.5" << std::endl
				<< "  [--overwrite_output] -> flag to overwrite output if exists" << std::endl
				<< "  [--ignore_thin] -> flag to ignore thinning step!" << std::endl
				<< "  [--ignore_multi_thresh_thin] -> flag to ignore multi thresh thinning step!" << std::endl
				<< std::endl
				<< "Version compiled on : " << __DATE__ << std::endl;

			printed_help = true;
		}


		void Parameters::parse(int argc, char* argv[])
		{
			if (argc == 1) {
				help();
				return;
			}

			std::list<std::string> unknown_args;

			size_t r = 1;
			while (r < argc) {
				std::string arg = argv[r];
				if (arg == "-h" || arg == "--help") {
					help();
					return;
				}
				else if (arg == "-i" && r + 1 < argc) {
					input_raster_path = argv[r + 1];
					r += 2;
				}
				else if (arg == "-ith" && r + 1 < argc) {
					input_thin_raster_path = argv[r + 1];
					r += 2;
				}
				else if ((arg == "-o" || arg == "--output") && r + 1 < argc) {
					std::string f = argv[r + 1];
					std::string extension = (f.size() > 4 ? f.substr(f.size() - 4, 4) : "");
					if (extension == ".shp" || extension == ".SHP") {
						output_shapefile = f;
						output_basename = f.substr(0, f.size() - 4);
					}
					else {
						std::cout << "Warning : invalid output filename. Writing result in result.shp" << std::endl;
					}
					r += 2;

				}
				else if ((arg == "-min_heat_value" || arg == "-mhv") && r + 1 < argc) {
					min_heat_value = atoi(argv[r + 1]);
					r += 2;
				}
				else if ((arg == "-max_n_level_neighbour" || arg == "-mnln") && r + 1 < argc) {
					max_n_level_neighbour = atoi(argv[r + 1]);
					r += 2;
				}
				else if ((arg == "-simplification_factor" || arg == "-sf") && r + 1 < argc) {
					simplification_factor = atof(argv[r + 1]);
					r += 2;
				}
				else if (arg == "--overwrite_output") {
					overwrite_output = true;
					r += 1;
				}
				else if (arg == "--ignore_thin") {
					ignore_thin = true;
					r += 1;
				}
				else if (arg == "--ignore_multi_thresh_thin") {
					ignore_multi_thresh_thin = true;
					r += 1;
				}
				else {
					unknown_args.push_back(arg);
					r += 1;
				}
			}

			if (!unknown_args.empty()) {
				std::cout << "There were unknown arguments in command line call :" << std::endl << '\t';
				for (const std::string& arg : unknown_args) std::cout << arg << " ";
				std::cout << std::endl;
				help();
			}
		}

		Parameters* params = nullptr;
	}
}