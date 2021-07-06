#pragma once
#include "defs.h"
#include "parameters.h"

namespace LxGeo
{

	namespace faultToVector
	{
		//template <typename mat_type, typename std::enable_if<std::is_arithmetic<mat_c_type>::value>::type* = nullptr>
		//template <typename mat_c_type>
		class PixGraph
		{
		public:
			PixGraph(){};

			PixGraph(matrix& _image_matrix, double _pix_width, double _pix_height, Boost_Point_2 _origin_point, bool _eight_connection) :
				image_matrix(_image_matrix),pix_width(_pix_width), pix_height(_pix_height), origin_point(_origin_point), eight_connection(_eight_connection) //, min_heat_value(_min_heat_value)
			{
				assert(_image_matrix.type() == CV_8U);
				if (eight_connection) {
					pixel_connection_pairs = { {-1, 0}, {1, 0}, {0, -1}, {0, 1}, 
					{-1, -1}, {1, 1}, {-1, 1}, {1, -1}, };
				}
				else {
					pixel_connection_pairs = { {-1, 0}, {1, 0}, {0, -1}, {0, 1} };
				}
				// padding input matirx
				int pad_value = 1;
				min_heat_value = static_cast<uchar>(params->min_heat_value);
				cv::copyMakeBorder(_image_matrix, image_matrix, pad_value, pad_value, pad_value, pad_value, cv::BORDER_CONSTANT, params->min_heat_value);
			};
			~PixGraph() {
			};

			void PixGraph::Init_graph();

			void PixGraph::Init_graph_2();

			size_t PixGraph::get_node_id_from_position(size_t c_row, size_t c_col, size_t max_row);

			bool PixGraph::add_neighbours(size_t reference_vertex_id, std::vector<size_t>& neighbours_indices, bool increment_level=true);

			void PixGraph::N_level_neighbours_iterate();

			void PixGraph::compute_pixels_angles();


		public:
			//PixGraphStrategy* pxg_strategy;
			// input matrix to vectorize
			matrix image_matrix;
		private:
			// minimum value to include (default 0)
			uchar min_heat_value=0;
			double pix_width=1;
			double pix_height=1;
			bool eight_connection;
			Boost_Point_2 origin_point = Boost_Point_2(0,0);
			std::vector<Boost_Point_2> all_vertcies;
			std::vector<uchar> all_vertcies_heat;
			// map of N level neighbours where key: vertex id & value: pair of [ pair of [ ( neighbors vertex id ) as vector; ( neighbors vertex id ) as set ] ; (indices of levels) as vector]
			std::map<size_t, std::pair< std::pair<std::vector<size_t>, std::set<size_t>> , std::vector<short int> > > vertecies_N_level_neighbors;

			std::vector<double> all_vertcies_angles;
			std::vector<double> all_vertcies_angle_homegenity;
			std::vector<std::pair<short int, short int>> pixel_connection_pairs;
			Boost_RTree_2 points_tree;
		};

		/*class PixGraphStrategy{};

		template <typename mat_c_type>
		class IPixGraphStrategy : public PixGraphStrategy
		{
		public:
			mat_c_type min_heat_value;
			std::vector<mat_c_type> all_vertcies_heat;
		public:
			IPixGraphStrategy() {}
			~IPixGraphStrategy() {}			
			void get_matrix_at(matrix& mat, size_t i, size_t j, mat_c_type& pix_value) {
				pix_value = mat.at<value_type>(i, j);
			}

			void set_matrix_at(matrix& mat, size_t i, size_t j, mat_c_type& pix_value) {
				mat.at<value_type>(i, j) = pix_value;
			}

			mat_c_type get_min_heat_value() { return static_cast<value_type>(min_heat_value); }
		};*/

	}

}