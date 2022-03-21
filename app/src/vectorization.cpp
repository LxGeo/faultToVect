#include "vectorization.h"
#include "io_raster.h"
#include "defs.h"
#include "parameters.h"
#include "raster_transformer.h"
#include <boost/log/trivial.hpp>
#include "pixgraph.h"
#include "vPixGraph/v_pix_graph.h"
#include "io_shapefile.h"


namespace LxGeo
{

	using namespace IO_DATA;
	namespace faultToVector
	{

		bool Vectorization::pre_check() {

			std::cout << "Pre check input parameters" << std::endl;
			RasterIO temp_raster = RasterIO();
			if (!temp_raster.load_raster(params->input_raster_path.c_str()))
				return false;
			
			if (temp_raster.band_count != 1) {
				std::cout << "Input raster bands count is different from 1! Check input data!" << std::endl;
				return false;
			}

			//output dirs creation
			boost::filesystem::path output_path(params->output_shapefile);
			boost::filesystem::path output_parent_dirname = output_path.parent_path();
			boost::filesystem::path output_temp_path = output_parent_dirname; // params->temp_dir;
			params->temp_dir = output_temp_path.string();
			if (boost::filesystem::exists(output_parent_dirname) || boost::filesystem::create_directory(output_parent_dirname))
			{
				std::cout << "Directory Created: " << output_parent_dirname.c_str() << std::endl;
			}
			else {
				std::cout << "Cannot create output directory: {}!" << output_parent_dirname.c_str() << std::endl;
				return false;
			}
			if (boost::filesystem::exists(output_temp_path) || boost::filesystem::create_directory(output_temp_path))
			{
				std::cout << "Directory Created: {}" << output_temp_path.c_str() << std::endl;
			}
			else {
				std::cout << "Cannot create temporary output directory: {}!" << output_temp_path.c_str() << std::endl;
				return false;
			}

			if (boost::filesystem::exists(output_path) && !params->overwrite_output) {
				std::cout << "output shapefile already exists: {}!" << output_path.c_str() << std::endl;
				std::cout << "Add --overwrite_output !" << output_path.c_str() << std::endl;
				return false;
			}
			

			return true;

		}

		void Vectorization::run() {

			// load raster
			RasterIO pred_raster = RasterIO();
			pred_raster.load_raster(params->input_raster_path, GA_ReadOnly, false);

			// multithresh thinning
			if (!(params->ignore_multi_thresh_thin)) {
				std::vector<double> thresh_linespace= { 0.5,0.6,0.7,0.8,0.85,0.9,0.92,0.95,0.99 };// = numcpp::linspace(0.1, 1.0, 10);
				std::list<double> thresh_values(thresh_linespace.begin(), thresh_linespace.end());//{ 0.5,0.6,0.7,0.8,0.85,0.9,0.92,0.95,0.99 };
				matrix multi_thresh_thinned_pred = multi_threshold_thinning(pred_raster.raster_data, thresh_values);
				RasterIO multi_thresh_thinned_raster = RasterIO(pred_raster, multi_thresh_thinned_pred);
				std::string multi_thresh_thinned_path = (boost::filesystem::path(params->temp_dir) / "multi_thresh_thinned.tif").string();
				multi_thresh_thinned_raster.write_raster(multi_thresh_thinned_path, true);
				pred_raster.raster_data = multi_thresh_thinned_pred;
			};

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
				std::vector<LineString_with_attributes> jl = ST.export_attachedLineString_as_LSwithAttr(ST.extract_junction_lines(c_sub_graph, 2, 4), sub_gr_id);
				all_traced.insert(all_traced.end(), jl.begin(), jl.end());
				sub_gr_id++;
				if (sub_gr_id == 1) break;
			}
			boost::filesystem::path traced_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("traced.shp");
			LineStringShapfileIO shp_traced = LineStringShapfileIO(traced_output_path.string(), nullptr);
			shp_traced.write_linestring_shapefile(all_traced);
			
			return;
			PixGraph pix_graph = PixGraph(thinned_raster.raster_data, thinned_raster.get_pixel_width(), thinned_raster.get_pixel_height(), thinned_raster.get_origin_point(),pred_raster.geotransform , true);
			pix_graph.Init_graph_2();
			pix_graph.compute_pixels_angles();
			pix_graph.compute_pixels_angles_homogenity2();

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

			pix_graph.compute_connected_components_iter();
			/*save connected components matrix*/
			matrix connected_comp_matrix = matrix(pred_raster.raster_Y_size, pred_raster.raster_X_size, CV_32SC1, cv::Scalar(0));
			pix_graph.transform_relative_vector_to_matrix(*pix_graph.get_all_vertcies_connected_comp_labels(), connected_comp_matrix);
			RasterIO connected_comp_raster = RasterIO(pred_raster, connected_comp_matrix);
			boost::filesystem::path connected_comp_raster_output_path = boost::filesystem::path(params->temp_dir) /
				boost::filesystem::path("connected_comp.tif");
			connected_comp_raster.write_raster(connected_comp_raster_output_path.string().c_str(), true);

			//pix_graph.generate_unidirectonal_compnents();
			pix_graph.generate_free_segments();
			pix_graph.write_free_segments_shapefile(params->output_shapefile, pred_raster.spatial_refrence);
		}


	}
}