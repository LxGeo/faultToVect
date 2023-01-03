#include "vectorization.h"
#include "io_raster.h"
#include "defs.h"
#include "parameters.h"
#include "raster_transformer.h"
#include <boost/log/trivial.hpp>
#include "pixgraph.h"
#include "vPixGraph/v_pix_graph.h"
#include "io_shapefile.h"
#include "affine_geometry/affine_transformer.h"


namespace LxGeo
{

	using namespace IO_DATA;
	namespace faultToVector
	{

		bool Vectorization::pre_check() {

			std::cout << "Pre check input parameters" << std::endl;
			RasterIO temp_raster = RasterIO();

			if (!params->input_raster_path.empty()) {				
				if (!temp_raster.load_raster(params->input_raster_path.c_str()))
					return false;
			}
			else if (!params->input_thin_raster_path.empty()) {
				if (!temp_raster.load_raster(params->input_thin_raster_path.c_str()))
					return false;
				apply_thin = false;
			}
			else {
				std::cout << "Should at least provide input raster or thin input raster!" << std::endl;
				return false;
			}
			
			if (temp_raster.band_count != 1) {
				std::cout << "Input raster bands count is different from 1! Check input data!" << std::endl;
				return false;
			}

			//output dirs creation
			boost::filesystem::path output_path(params->output_shapefile);
			boost::filesystem::path output_parent_dirname = output_path.parent_path();
			boost::filesystem::path output_temp_path = output_parent_dirname; // params->temp_dir;
			params->temp_dir = output_temp_path.string();
			if (boost::filesystem::exists(output_parent_dirname) || boost::filesystem::create_directories(output_parent_dirname))
			{
				std::cout << "Directory Created: " << output_parent_dirname.string() << std::endl;
			}
			else {
				std::cout << "Cannot create output directory: " << output_parent_dirname.string() << std::endl;
				return false;
			}
			if (boost::filesystem::exists(output_temp_path) || boost::filesystem::create_directories(output_temp_path))
			{
				std::cout << "Directory Created: " << output_temp_path.string() << std::endl;
			}
			else {
				std::cout << "Cannot create temporary output directory: " << output_temp_path.string() << std::endl;
				return false;
			}

			if (boost::filesystem::exists(output_path) && !params->overwrite_output) {
				std::cout << "output shapefile already exists: " << output_path.string() << std::endl;
				std::cout << "Add --overwrite_output !" << output_path.string() << std::endl;
				return false;
			}
			

			return true;

		}

		void Vectorization::run() {

			RasterIO thinned_raster;

			if (apply_thin) {
				// load preds raster
				RasterIO pred_raster = RasterIO();
				pred_raster.load_raster(params->input_raster_path, GA_ReadOnly, false);

				// binarize
				matrix thresholded_binary_matrix;
				threshold(pred_raster.raster_data, thresholded_binary_matrix, 0.7, 255, cv::THRESH_BINARY);
				thresholded_binary_matrix.convertTo(thresholded_binary_matrix, CV_8U);
				RasterIO thresholded_binary_raster = RasterIO(pred_raster, thresholded_binary_matrix);
				boost::filesystem::path thresholded_binary_output_path = boost::filesystem::path(params->temp_dir) /
					boost::filesystem::path("thresholded_binary.tif");
				thresholded_binary_raster.write_raster(thresholded_binary_output_path.string().c_str(), true);

				// fix artefacts (fill holes and close openings)
				matrix fixed_binary_matrix = 255 - thresholded_binary_matrix.clone();
				int min_size = 0, connectivity = 8;
				if (min_size>0)
					fill_small_components(fixed_binary_matrix, 2, 4); // small holes removal
				fill_small_components(fixed_binary_matrix, min_size, connectivity);
				fixed_binary_matrix = 255 - fixed_binary_matrix;
				fill_small_components(fixed_binary_matrix, min_size, connectivity);
				boost::filesystem::path fixed_binary_output_path = boost::filesystem::path(params->temp_dir) /
					boost::filesystem::path("fixed_binary.tif");
				RasterIO fixed_binary_raster = RasterIO(pred_raster, fixed_binary_matrix);
				fixed_binary_raster.write_raster(fixed_binary_output_path.string().c_str(), true);

				// Apply thinning on binary
				matrix thinned_binary;
				cv::ximgproc::thinning(fixed_binary_matrix, thinned_binary, cv::ximgproc::THINNING_GUOHALL);
				thinned_raster = RasterIO(pred_raster, thinned_binary);
				boost::filesystem::path thinned_binary_output_path = boost::filesystem::path(params->temp_dir) /
					boost::filesystem::path("thinned_binary.tif");
				thinned_raster.write_raster(thinned_binary_output_path.string().c_str(), true);
			}
			else {
				thinned_raster.load_raster(params->input_thin_raster_path, GA_ReadOnly, false);
			}
						
			SkeletonTracer ST = SkeletonTracer(thinned_raster.raster_data, 0);
			std::vector<LineString_with_attributes> out_edges = ST.export_edge_graph_as_LSwithAttr();
			boost::filesystem::path graph_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("graph.shp");
			LineStringShapfileIO shp = LineStringShapfileIO(graph_output_path.string(), nullptr);
			shp.write_linestring_shapefile(out_edges);

			std::vector<LineString_with_attributes> all_traced;
			auto sub_graphs = ST.connected_components_subgraphs();
			size_t sub_gr_id = 0;
			for (auto& c_sub_graph : sub_graphs) {
				std::vector<LineString_with_attributes> jl = ST.export_attachedLineString_as_LSwithAttr(ST.extract_junction_lines_iterative(c_sub_graph, 2, 6), sub_gr_id);
				all_traced.insert(all_traced.end(), jl.begin(), jl.end());
				std::cout << "Cnt: " << all_traced.size() << std::endl;
				sub_gr_id++;
			}
			auto coords_transformet_mat = thinned_raster.get_matrix_transformer(thinned_raster.get_pixel_width() / 2, thinned_raster.get_pixel_height() / 2);
			for (auto& lswa : all_traced)
				lswa.set_definition(affine_transform_geometry<Boost_LineString_2, Boost_LineString_2>(lswa.get_definition(), coords_transformet_mat));
			boost::filesystem::path traced_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("traced_or2_c6.shp");
			LineStringShapfileIO shp_traced = LineStringShapfileIO(traced_output_path.string(), thinned_raster.spatial_refrence);
			shp_traced.write_linestring_shapefile(all_traced);			
			
			return;
			PixGraph pix_graph = PixGraph(thinned_raster.raster_data, thinned_raster.get_pixel_width(), thinned_raster.get_pixel_height(), thinned_raster.get_origin_point(), thinned_raster.geotransform , true);
			pix_graph.Init_graph_2();
			pix_graph.compute_pixels_angles();
			pix_graph.compute_pixels_angles_homogenity2();

			/*save angles matrix*/
			matrix angle_matrix = matrix(thinned_raster.raster_Y_size, thinned_raster.raster_X_size, CV_64F, cv::Scalar(0.f));
			pix_graph.transform_relative_vector_to_matrix(*pix_graph.get_all_vertcies_angle(), angle_matrix);
			RasterIO angles_raster = RasterIO(thinned_raster, angle_matrix);
			boost::filesystem::path angles_raster_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("angles.tif");
			angles_raster.write_raster(angles_raster_output_path.string().c_str(), true);

			/*save angles homogenity matrix*/
			matrix angle_homogenity_matrix = matrix(thinned_raster.raster_Y_size, thinned_raster.raster_X_size, CV_64F, cv::Scalar(0.f));
			pix_graph.transform_relative_vector_to_matrix(*pix_graph.get_all_vertcies_angle_homegenity(), angle_homogenity_matrix);
			RasterIO angles_homogenity_raster = RasterIO(thinned_raster, angle_homogenity_matrix);
			boost::filesystem::path angles_homogenity_raster_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("angles_homogenity.tif");
			angles_homogenity_raster.write_raster(angles_homogenity_raster_output_path.string().c_str(), true);

			pix_graph.compute_connected_components_iter();
			/*save connected components matrix*/
			matrix connected_comp_matrix = matrix(thinned_raster.raster_Y_size, thinned_raster.raster_X_size, CV_32SC1, cv::Scalar(0));
			pix_graph.transform_relative_vector_to_matrix(*pix_graph.get_all_vertcies_connected_comp_labels(), connected_comp_matrix);
			RasterIO connected_comp_raster = RasterIO(thinned_raster, connected_comp_matrix);
			boost::filesystem::path connected_comp_raster_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("connected_comp.tif");
			connected_comp_raster.write_raster(connected_comp_raster_output_path.string().c_str(), true);

			//pix_graph.generate_unidirectonal_compnents();
			pix_graph.generate_free_segments();
			pix_graph.write_free_segments_shapefile(params->output_shapefile, thinned_raster.spatial_refrence);
		}


	}
}