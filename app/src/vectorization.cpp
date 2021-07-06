#include "vectorization.h"
#include "io_raster.h"
#include "defs.h"
#include "parameters.h"
#include "raster_transformer.h"
#include <boost/log/trivial.hpp>
#include "pixgraph.h"


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

			// load raster
			RasterIO pred_raster = RasterIO();
			pred_raster.load_raster(params->input_raster_path, GA_ReadOnly, false);

			// Apply thinning on grayscale
			RasterIO thinned_raster = RasterIO(pred_raster);
			if (!(params->ignore_thin))
				thin_raster(pred_raster, thinned_raster);
			else thinned_raster.raster_data = pred_raster.raster_data;

			/*boost::filesystem::path thinned_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("thinned.tif");
			thinned_raster.write_raster(thinned_output_path.string().c_str(), true);*/
			
			thinned_raster.raster_data.convertTo(thinned_raster.raster_data, CV_8U);
			assert(thinned_raster.raster_data.type() == CV_8U);
			PixGraph pix_graph = PixGraph(thinned_raster.raster_data, thinned_raster.raster_X_size, thinned_raster.raster_Y_size, thinned_raster.get_origin_point(), true);
			pix_graph.Init_graph_2();
			pix_graph.compute_pixels_angles();

		}


	}
}