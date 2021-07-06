#include "raster_transformer.h"
#include "defs.h"
#include "io_raster.h"

#include "parameters.h"
namespace LxGeo
{

	using namespace IO_DATA;
	namespace faultToVector
	{

		void thin_raster(RasterIO& raster_to_thin, RasterIO& thinned_raster, int min_value) {
			
			assert(raster_to_thin.raster_size == thinned_raster.raster_size &&
				"raster to thin size is different from thinned raster!");

			matrix binary_image;
			cv::threshold(raster_to_thin.raster_data, binary_image, min_value, 255, cv::THRESH_BINARY);
			binary_image.convertTo(binary_image, CV_8U);

			/*RasterIO bin_raster(raster_to_thin, binary_image);
			boost::filesystem::path binary_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("binary.tif");
			bin_raster.write_raster(binary_output_path.string().c_str(), true);*/

			matrix thinned_binary;// = matrix::zeros(binary_image.rows, binary_image.cols, CV_8UC1);
			cv::ximgproc::thinning(binary_image, thinned_binary, cv::ximgproc::THINNING_ZHANGSUEN);
			

			/*RasterIO bin_thin_raster(raster_to_thin, thinned_binary);
			boost::filesystem::path binary_thin_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("binary_thinned.tif");
			bin_thin_raster.write_raster(binary_thin_output_path.string().c_str(), true);*/

			
			//raster_to_thin.raster_data.copyTo(thinned_raster.raster_data, thinned_binary);
			thinned_binary.convertTo(thinned_binary, raster_to_thin.raster_data.type());
			cv::bitwise_and(thinned_binary, raster_to_thin.raster_data, thinned_raster.raster_data);
			
			/*boost::filesystem::path temp_thin_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("temp_thinned.tif");
			thinned_raster.write_raster(temp_thin_output_path.string().c_str(), true);*/

		}

	}
}