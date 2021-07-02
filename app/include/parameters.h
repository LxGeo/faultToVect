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

			bool overwrite_output;


		};

		extern Parameters* params;
	}
}
