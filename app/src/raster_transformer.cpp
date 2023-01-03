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
			cv::ximgproc::thinning(binary_image, thinned_binary, cv::ximgproc::THINNING_GUOHALL);
			
            int morph_size = 1;

            // Create structuring element
            matrix element = cv::getStructuringElement(
                cv::MORPH_RECT,
                cv::Size(2 * morph_size + 1,
                    2 * morph_size + 1),
                cv::Point(morph_size, morph_size));
            cv::morphologyEx(thinned_binary, thinned_binary, cv::MORPH_CLOSE, element);
            cv::ximgproc::thinning(thinned_binary, thinned_binary, cv::ximgproc::THINNING_GUOHALL);

            /*RasterIO closed_bin_raster(raster_to_thin, binary_image);
            boost::filesystem::path closed_binary_output_path = boost::filesystem::path(params->temp_dir) /
                boost::filesystem::path("binary_closed.tif");
            closed_bin_raster.write_raster(closed_binary_output_path.string().c_str(), true);*/

			RasterIO bin_thin_raster(raster_to_thin, thinned_binary);
			boost::filesystem::path binary_thin_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("binary_thinned.tif");
			bin_thin_raster.write_raster(binary_thin_output_path.string().c_str(), true);
            			
			raster_to_thin.raster_data.copyTo(thinned_raster.raster_data, thinned_binary);
			//thinned_binary.convertTo(thinned_binary, raster_to_thin.raster_data.type());
			//cv::bitwise_and(thinned_binary, raster_to_thin.raster_data, thinned_raster.raster_data);
			
			/*boost::filesystem::path temp_thin_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("temp_thinned.tif");
			thinned_raster.write_raster(temp_thin_output_path.string().c_str(), true);*/

		}

		

		matrix multi_threshold_thinning(matrix& in_probability_matrix, std::list<double> threshold_list, std::list<double> thresh_weigths) {

			assert(threshold_list.size() > 1, "Theshhold values list must have at least 2 values!");

			if (thresh_weigths.size() == 0) thresh_weigths = std::list<double>(threshold_list.size(), 1);

			assert(threshold_list.size() == thresh_weigths.size(), "threshold values and weights size do not match!");

			matrix aggregated(in_probability_matrix.size(), CV_32F, cv::Scalar(0));
			double weights_sum = 0;

			auto thresh_value_iter = threshold_list.begin();
			auto thresh_weight_iter = thresh_weigths.begin();

			for (; thresh_value_iter != threshold_list.end(); ++thresh_value_iter, ++thresh_weight_iter) {
				matrix c_binary_map, thinned_binary;
				threshold(in_probability_matrix, c_binary_map, *thresh_value_iter, 255, cv::THRESH_BINARY);
				c_binary_map.convertTo(c_binary_map, CV_8U);
				cv::ximgproc::thinning(c_binary_map, thinned_binary, cv::ximgproc::THINNING_ZHANGSUEN);

				
				/*RasterIO ref_raster = RasterIO(params->input_raster_path);
				RasterIO multi_thresh_thinned_raster = RasterIO(ref_raster, thinned_binary);
				std::string multi_thresh_thinned_path = (boost::filesystem::path(params->output_shapefile).parent_path() / "th_binary.tif").string();
				multi_thresh_thinned_raster.write_raster(multi_thresh_thinned_path, true);
				*/
				cv::accumulate(thinned_binary * (*thresh_weight_iter)/255, aggregated);
				weights_sum += *thresh_weight_iter;
			}
			aggregated.convertTo(aggregated, CV_8U);
			return aggregated;

		}


		void fill_small_components(matrix& in_binary, int min_size, int connectivity) {

			matrix labels, stats, centroids;

			cv::connectedComponentsWithStats(in_binary, labels, stats, centroids, connectivity, CV_32S);
			int x, y, w, h, area;

			auto invert_pixels = [&](int& comp_id) {
				for (int c_col = x; c_col <= x+w; c_col++) {
					for (int c_row = y; c_row <= y+h; c_row++) {
						if (labels.at<int>(c_row, c_col) == comp_id)
							in_binary.at<uchar>(c_row, c_col) = 255 - in_binary.at<uchar>(c_row, c_col);
					}
				}
			};

			for (int i = 0; i < stats.rows; i++)
			{
				x = stats.at<int>(cv::Point(0, i));
				y = stats.at<int>(cv::Point(1, i));
				w = stats.at<int>(cv::Point(2, i));
				h = stats.at<int>(cv::Point(3, i));
				area = stats.at<int>(cv::Point(4, i));
				if (area <= min_size) {
					invert_pixels(i);
				}
			}

		}

	}
}