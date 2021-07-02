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
			bool is_initialized = !input_raster_path.empty();
			if (!is_initialized) help();

			return is_initialized;
		}


		void Parameters::init()
		{
			printed_help = false;

			input_raster_path.clear();

			output_basename = "result";
			output_shapefile = "result.shp";
			temp_dir = "temp_dir";
			overwrite_output = false;

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
				<< "  [--overwrite_output] -> flag to overwrite output if exists" << std::endl
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
				else if (arg == "--overwrite_output" && r + 1 < argc) {
					overwrite_output = true;
					r += 2;
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