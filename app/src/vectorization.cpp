#include "vectorization.h"
#include "io_raster.h"
#include "defs.h"
#include "parameters.h"


namespace LxGeo
{

	using namespace IO_DATA;
	namespace faultToVector
	{

		bool Vectorization::pre_check() {

			BOOST_LOG_TRIVIAL(info) << "Pre check input parameters";
			RasterIO temp_raster = RasterIO();
			if (!temp_raster.load_raster(params->input_raster_path.c_str()))
				return false;
			
			if (temp_raster.band_count != 1) {
				BOOST_LOG_TRIVIAL(fatal) << "Input raster bands count is different from 1! Check input data!";
				return false;
			}

			//output dirs creation
			boost::filesystem::path output_path(params->output_shapefile);
			boost::filesystem::path output_parent_dirname = output_path.parent_path();
			boost::filesystem::path output_temp_path = output_parent_dirname / params->temp_dir;
			params->temp_dir = output_temp_path.string();
			if (boost::filesystem::exists(output_parent_dirname) || boost::filesystem::create_directory(output_parent_dirname))
			{
				BOOST_LOG_TRIVIAL(info) << fmt::format("Directory Created: {}", output_parent_dirname.string());
			}
			else {
				BOOST_LOG_TRIVIAL(fatal) << fmt::format("Cannot create output directory: {}!", output_parent_dirname.string());
				return false;
			}
			if (boost::filesystem::exists(output_temp_path) || boost::filesystem::create_directory(output_temp_path))
			{
				BOOST_LOG_TRIVIAL(info) << fmt::format("Directory Created: {}", output_temp_path.string());
			}
			else {
				BOOST_LOG_TRIVIAL(fatal) << fmt::format("Cannot create temporary output directory: {}!", output_temp_path.string());
				return false;
			}

			if (boost::filesystem::exists(output_path) && !params->overwrite_output) {
				BOOST_LOG_TRIVIAL(fatal) << fmt::format("output shapefile already exists: {}!", output_path.string());
				BOOST_LOG_TRIVIAL(fatal) << fmt::format("Add --overwrite_output !", output_path.string());
				return false;
			}
			

			return true;

		}

		void Vectorization::run() {

		}


	}
}