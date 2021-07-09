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
			matrix thinned_zeros = matrix::zeros(pred_raster.raster_Y_size, pred_raster.raster_X_size, pred_raster.raster_data.type());
			RasterIO thinned_raster = RasterIO(pred_raster, thinned_zeros);
			if (!(params->ignore_thin)){
				thin_raster(pred_raster, thinned_raster);
				boost::filesystem::path thinned_output_path = boost::filesystem::path(params->temp_dir) /
					boost::filesystem::path("thinned.tif");
				thinned_raster.write_raster(thinned_output_path.string().c_str(), true);
			}
			else thinned_raster.raster_data = pred_raster.raster_data;
			
			thinned_raster.raster_data.convertTo(thinned_raster.raster_data, CV_8U);
			assert(thinned_raster.raster_data.type() == CV_8U);
			PixGraph pix_graph = PixGraph(thinned_raster.raster_data, thinned_raster.get_pixel_width(), thinned_raster.get_pixel_height(), thinned_raster.get_origin_point(),pred_raster.geotransform , true);
			pix_graph.Init_graph_2();
			pix_graph.compute_pixels_angles();
			pix_graph.compute_pixels_angles_homogenity();
			pix_graph.generate_free_segments();
			pix_graph.write_free_segments_shapefile(params->output_shapefile, pred_raster.spatial_refrence->Clone());

			/*save angles matrix*/
			matrix angle_matrix = matrix(pred_raster.raster_Y_size, pred_raster.raster_X_size, CV_64F, cv::Scalar(0.f));
			pix_graph.transform_relative_vector_to_matrix(*pix_graph.get_all_vertcies_angle(), angle_matrix);
			RasterIO angles_raster = RasterIO(pred_raster, angle_matrix);
			boost::filesystem::path angles_raster_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("angles.tif");
			angles_raster.write_raster(angles_raster_output_path.string().c_str(), true);

			/*save angles homogenity matrix*/
			matrix angle_homogenity_matrix = matrix(pred_raster.raster_Y_size, pred_raster.raster_X_size, CV_64F, cv::Scalar(0.f));
			pix_graph.transform_relative_vector_to_matrix(*pix_graph.get_all_vertcies_angle_homegenity(), angle_homogenity_matrix);
			RasterIO angles_homogenity_raster = RasterIO(pred_raster, angle_homogenity_matrix);
			boost::filesystem::path angles_homogenity_raster_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("angles_homogenity.tif");
			angles_homogenity_raster.write_raster(angles_homogenity_raster_output_path.string().c_str(), true);

		}


	}
}