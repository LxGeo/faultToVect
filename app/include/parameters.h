#pragma once
#include "defs.h"
#include <boost/filesystem.hpp>


namespace LxGeo
{
	namespace faultToVector
	{
		class Parameters
		{
		public:
			Parameters(int argc, char* argv[]);

			~Parameters();

			bool initialized();

		protected:
			void init();

			void help();

			void parse(int argc, char* argv[]);

		public:
			bool printed_help;

			std::string input_raster_path;

			std::string output_basename;
			std::string output_shapefile;
			std::string temp_dir;

			int min_heat_value;
			short int max_n_level_neighbour;

			bool overwrite_output;
			bool ignore_thin;
			bool ignore_multi_thresh_thin;
			
			double simplification_factor;


		};

		extern Parameters* params;
	}
}
