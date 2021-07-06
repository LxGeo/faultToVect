#include "pixgraph.h"

namespace LxGeo
{

	namespace faultToVector
	{
		void PixGraph::Init_graph(){
			
			// count possible pixel
			size_t max_pixels_cnt = size_t(image_matrix.rows) * size_t(image_matrix.cols);
			// reserve vectors
			all_vertcies.reserve(max_pixels_cnt);
			size_t current_nodes_count = 0;
			// Iterate through all inner pixels
			for (size_t current_row = 1; current_row < image_matrix.rows; ++current_row) {
				for (size_t current_col = 1; current_col < image_matrix.cols; ++current_col) {
					
					uchar c_pix_heat = image_matrix.at<uchar>(current_row, current_col);
					if (c_pix_heat > min_heat_value) {
						current_nodes_count++;
						size_t current_node_id = get_node_id_from_position(current_row, current_col, image_matrix.rows);
						all_vertcies.push_back(Boost_Point_2(current_row-1, current_col-1));						
						for (auto connection_pair : pixel_connection_pairs) {
							size_t c_neighbour_row = current_row + connection_pair.first;
							size_t c_neighbour_col = current_col + connection_pair.second;
							uchar c_neighbour_heat = image_matrix.at<uchar>(c_neighbour_row, c_neighbour_col);
							if (!(c_neighbour_heat > min_heat_value)) continue;
							size_t neigbour_node_id = get_node_id_from_position(c_neighbour_row, c_neighbour_col, image_matrix.rows);
							// add edge if not existant in set
							auto set_addition = vertecies_N_level_neighbors[current_node_id].first.second.insert(neigbour_node_id);
							if (set_addition.second)
								vertecies_N_level_neighbors[current_node_id].first.first.push_back(neigbour_node_id);
						}
					}
				}
			}
			all_vertcies.shrink_to_fit();
			all_vertcies_angles.reserve(all_vertcies.size());
			all_vertcies_angle_homegenity.reserve(all_vertcies.size());

		}

		void PixGraph::Init_graph_2(){
		
			// count possible pixel
			size_t max_pixels_cnt = size_t(image_matrix.rows) * size_t(image_matrix.cols);
			// reserve vectors
			all_vertcies.reserve(max_pixels_cnt);
			// Iterate through all inner pixels
			for (size_t current_row = 1; current_row < image_matrix.rows; ++current_row) {
				for (size_t current_col = 1; current_col < image_matrix.cols; ++current_col) {
				
					uchar c_pix_heat = image_matrix.at<uchar>(current_row, current_col);
					if (c_pix_heat > min_heat_value) {
						all_vertcies.push_back(Boost_Point_2(current_row - 1, current_col - 1));
						all_vertcies_heat.push_back(c_pix_heat);
					}

				}
			}
			all_vertcies.shrink_to_fit();
			// make rtree for points
			points_tree = Boost_RTree_2();
			for (size_t vertex_idx = 0; vertex_idx < all_vertcies.size(); ++vertex_idx)
				points_tree.insert(Boost_Value_2(Boost_Box_2(all_vertcies[vertex_idx], all_vertcies[vertex_idx]), vertex_idx));

			// add first n level neighbour
			for (size_t c_vertex_idx = 0; c_vertex_idx < all_vertcies.size(); ++c_vertex_idx) {
				
				Boost_Point_2 c_point = all_vertcies[c_vertex_idx];

				// search for neighbours touching // meaning distance <= 1 if 4 connection or distance <= sqrt(2) if 8 connection
				double buff_value = (eight_connection) ? 1.f : sqrt(2);
				Boost_Point_2 low_left_pt = Boost_Point_2(bg::get<0>(c_point) - buff_value, bg::get<1>(c_point) - buff_value);
				Boost_Point_2 top_right_pt = Boost_Point_2(bg::get<0>(c_point) + buff_value, bg::get<1>(c_point) + buff_value);

				std::vector<Boost_Value_2> possible_neighbors;
				Boost_Box_2 query(low_left_pt, top_right_pt);
				points_tree.query(bgi::intersects(query), std::back_inserter(possible_neighbors));

				std::vector<size_t> correct_neighbors;
				correct_neighbors.reserve(possible_neighbors.size());
				for (int i = 0; i < possible_neighbors.size(); i++) {

					size_t possible_i_index = possible_neighbors[i].second;
					Boost_Point_2 possible_i_point = all_vertcies[possible_i_index];

					double diff_distance = std::abs<double>(bg::distance(c_point, possible_i_point));
					if (diff_distance <= buff_value) {
						correct_neighbors.push_back(possible_i_index);
					}
				}
				correct_neighbors.shrink_to_fit();

				add_neighbours(c_vertex_idx, correct_neighbors);

			}

			// add remaining n level neighbours
			for (size_t n_level_fill_iteration = 2; n_level_fill_iteration <= params->max_n_level_neighbour; ++n_level_fill_iteration) {
				N_level_neighbours_iterate();
			}
		}

		size_t PixGraph::get_node_id_from_position(size_t c_row, size_t c_col, size_t max_row) {
			return max_row * c_row + c_col;
		}

		bool PixGraph::add_neighbours(size_t reference_vertex_id, std::vector<size_t>& neighbours_indices, bool increment_level) {

			bool all_not_added = true;
			for (size_t correct_neighbour_idx : neighbours_indices) {
				if (correct_neighbour_idx == reference_vertex_id) continue;
				// add edge if not existant in set
				auto set_addition = vertecies_N_level_neighbors[reference_vertex_id].first.second.insert(correct_neighbour_idx);
				all_not_added *= !set_addition.second;
				if (set_addition.second)
					vertecies_N_level_neighbors[reference_vertex_id].first.first.push_back(correct_neighbour_idx);
			}

			if (increment_level && !all_not_added)
				vertecies_N_level_neighbors[reference_vertex_id].second.push_back(vertecies_N_level_neighbors[reference_vertex_id].first.first.size());

			return all_not_added;
		}

		void PixGraph::N_level_neighbours_iterate() {
			for (size_t c_vertex_idx = 0; c_vertex_idx < all_vertcies.size(); ++c_vertex_idx) {
				if (vertecies_N_level_neighbors[c_vertex_idx].first.first.empty())
					continue;

				size_t last_neighbours_end_idx = vertecies_N_level_neighbors[c_vertex_idx].first.first.size();
				size_t last_neighbours_start_idx = last_neighbours_end_idx - vertecies_N_level_neighbors[c_vertex_idx].second[vertecies_N_level_neighbors[c_vertex_idx].second.size()-1];
				std::vector<size_t> last_added_neighbours(vertecies_N_level_neighbors[c_vertex_idx].first.first.begin()+ last_neighbours_start_idx,
					vertecies_N_level_neighbors[c_vertex_idx].first.first.end());

				// bool turns to false if at least one element is added
				bool all_not_added = true;
				for (auto c_last_added_neighbour : last_added_neighbours) {

					std::vector<size_t> c_last_added_first_neighbours(
						vertecies_N_level_neighbors[c_last_added_neighbour].first.first.begin(),
						vertecies_N_level_neighbors[c_last_added_neighbour].first.first.begin() + vertecies_N_level_neighbors[c_last_added_neighbour].second[0]-1 // check this
					);
					all_not_added *= add_neighbours(c_vertex_idx, c_last_added_first_neighbours, false);
				}
				// add level index only if new neighbours created
				if (!all_not_added)
					vertecies_N_level_neighbors[c_vertex_idx].second.push_back(vertecies_N_level_neighbors[c_vertex_idx].first.first.size());

			}
		}

		void PixGraph::compute_pixels_angles() {

			for (size_t c_vertex_idx = 0; c_vertex_idx < all_vertcies.size(); ++c_vertex_idx) {
				if (vertecies_N_level_neighbors[c_vertex_idx].first.first.empty())
				{
					all_vertcies_angles[c_vertex_idx] = double(0);
					continue;
				}

				std::vector<size_t> neighbours_indices = vertecies_N_level_neighbors[c_vertex_idx].first.first;
				//weighted angles notebook

			}

		};
	
	}
}