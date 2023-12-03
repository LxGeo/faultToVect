#include "pixgraph.h"
#include "defs.h"

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
			std::cout << "Loading Raster values" << std::endl;
			tqdm bar;
			for (size_t current_row = 1; current_row < image_matrix.rows; ++current_row) {
				bar.progress(current_row, image_matrix.rows);
				for (size_t current_col = 1; current_col < image_matrix.cols; ++current_col) {
				
					uchar c_pix_heat = image_matrix.at<uchar>(current_row, current_col);
					if (c_pix_heat > min_heat_value) {
						all_vertcies.push_back(Boost_Point_2(current_row - 0.5, current_col - 0.5));
						all_vertcies_heat.push_back(c_pix_heat);
					}

				}
			}
			bar.finish();
			all_vertcies.shrink_to_fit();
			// make rtree for points
			std::cout << "Creating points tree" << std::endl;
			bar.reset();
			points_tree = Boost_RTree_2();
			for (size_t vertex_idx = 0; vertex_idx < all_vertcies.size(); ++vertex_idx) {
				bar.progress(vertex_idx, all_vertcies.size());
				points_tree.insert(Boost_Value_2(Boost_Box_2(all_vertcies[vertex_idx], all_vertcies[vertex_idx]), vertex_idx));
			}
			bar.finish();
						
			// add first n level neighbour
			std::cout << "Creating first level neighbour" << std::endl;
			bar.reset();
			for (size_t c_vertex_idx = 0; c_vertex_idx < all_vertcies.size(); ++c_vertex_idx) {
				bar.progress(c_vertex_idx, all_vertcies.size());
				Boost_Point_2 c_point = all_vertcies[c_vertex_idx];

				// search for neighbours touching // meaning distance <= 1 if 4 connection or distance <= sqrt(2) if 8 connection
				double buff_value = (eight_connection) ? sqrt(2) : 1.f;
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
			bar.finish();
			// add remaining n level neighbours
			for (size_t n_level_fill_iteration = 2; n_level_fill_iteration <= params->max_n_level_neighbour; ++n_level_fill_iteration) {
				std::cout << "Running " << n_level_fill_iteration << " N level" << std::endl;
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
			tqdm bar;
			for (size_t c_vertex_idx = 0; c_vertex_idx < all_vertcies.size(); ++c_vertex_idx) {
				bar.progress(c_vertex_idx, all_vertcies.size());
				if (vertecies_N_level_neighbors[c_vertex_idx].first.first.empty())
					continue;

				size_t last_neighbours_end_idx = vertecies_N_level_neighbors[c_vertex_idx].first.first.size();
				size_t last_neighbours_start_idx = last_neighbours_end_idx - vertecies_N_level_neighbors[c_vertex_idx].second[vertecies_N_level_neighbors[c_vertex_idx].second.size()-1];
				std::vector<size_t> last_added_neighbours(vertecies_N_level_neighbors[c_vertex_idx].first.first.begin()+ last_neighbours_start_idx,
					vertecies_N_level_neighbors[c_vertex_idx].first.first.end());

				// bool turns to false if at least one element is added
				bool all_not_added = true;
				for (auto c_last_added_neighbour : last_added_neighbours) {
					auto c_first_neighbours = get_vertex_first_neighbours(c_last_added_neighbour);
					std::vector<size_t> c_last_added_first_neighbours(
						c_first_neighbours.first,
						c_first_neighbours.second// check this
					);
					all_not_added *= add_neighbours(c_vertex_idx, c_last_added_first_neighbours, false);
				}
				// add level index only if new neighbours created
				if (!all_not_added)
					vertecies_N_level_neighbors[c_vertex_idx].second.push_back(vertecies_N_level_neighbors[c_vertex_idx].first.first.size());

			}
			bar.finish();
		}

		std::pair<std::vector<size_t>::iterator, std::vector<size_t>::iterator> PixGraph::get_vertex_first_neighbours(size_t vertex_idx) {

			/*return {
				vertecies_N_level_neighbors[vertex_idx].first.first.begin(),
				vertecies_N_level_neighbors[vertex_idx].first.first.begin() + vertecies_N_level_neighbors[vertex_idx].second[0]
			};*/
			
			auto st = vertecies_N_level_neighbors[vertex_idx].first.first.begin();
			auto ed = vertecies_N_level_neighbors[vertex_idx].first.first.begin() + vertecies_N_level_neighbors[vertex_idx].second[0];
			return { st, ed };
		}

		void PixGraph::compute_pixels_angles() {

			all_vertcies_angles.reserve(all_vertcies.size());
			for (size_t c_vertex_idx = 0; c_vertex_idx < all_vertcies.size(); ++c_vertex_idx) {
				if (vertecies_N_level_neighbors[c_vertex_idx].first.first.empty())
				{
					all_vertcies_angles.push_back( M_PI);
					continue;
				}

				std::vector<size_t>& neighbours_indices = vertecies_N_level_neighbors[c_vertex_idx].first.first;
				std::vector<short int>& neighbours_levels = vertecies_N_level_neighbors[c_vertex_idx].second;
				std::vector<short int> flattened_neighbours_label;
				flattened_neighbours_label.reserve(neighbours_indices.size());
				short int start_label = 1;
				for (short int c_label : neighbours_levels) {
					for (short int i = 0; i < c_label; i++) flattened_neighbours_label.push_back(start_label);
					start_label++;
				}
				//weighted angles notebook

				std::vector<uchar> neigbours_heat = get_vector_by_indices(all_vertcies_heat, neighbours_indices);
				std::vector<Boost_Point_2> neighbours_points = get_vector_by_indices(all_vertcies, neighbours_indices);
				//transform Boost to cgal points
				std::vector<Inexact_Point_2> neighbours_points_cgal;
				container_transform_B2C_Points(neighbours_points, neighbours_points_cgal);

				// pair of central_pixel (point; heat_value)
				Inexact_Point_2 central_pixel_cgal_point = transform_B2C_Point(all_vertcies[c_vertex_idx]);
				std::pair<Inexact_Point_2, uchar> central_pixel_pair(central_pixel_cgal_point, all_vertcies_heat[c_vertex_idx]);
				Inexact_Vector_2 common_vector = get_best_fitting_vector(central_pixel_pair, neighbours_points_cgal, neigbours_heat, flattened_neighbours_label);
				// maybe normalize vector before cross product
				common_vector = common_vector / std::sqrt(common_vector.squared_length());
				//common_vector = Inexact_Vector_2(x_pixel_sign * common_vector.x(), y_pixel_sign * common_vector.y());
				double cross_p = refrence_axis_vector * common_vector;
				double angle_value = acos(cross_p);
				//double angle_min = asin(sqrt(sin(angle_value) * sin(angle_value)));
				all_vertcies_angles.push_back(angle_value);
			}

		};

		Inexact_Vector_2 PixGraph::get_best_fitting_vector(std::pair<Inexact_Point_2, uchar> central_pix_pair,
			std::vector<Inexact_Point_2> neighbours_points,
			std::vector<uchar> neigbours_heat, std::vector<short int>& level_labels) {

			Inexact_Point_2 central_pix_point = central_pix_pair.first;
			Inexact_Vector_2 common_vector= Inexact_Vector_2(central_pix_point, central_pix_point);
			for (size_t c_neighbour_index = 0; c_neighbour_index < neighbours_points.size(); ++c_neighbour_index) {
				size_t c_level_label = level_labels[c_neighbour_index];
				Inexact_Vector_2 c_vector = Inexact_Vector_2(central_pix_point, neighbours_points[c_neighbour_index]) ;
				size_t c_weight = neigbours_heat[c_neighbour_index];
				c_vector = (common_vector * c_vector) > 0 ? c_vector : -c_vector;
				common_vector += c_vector * c_weight/ c_level_label;
			}
			return common_vector;
		}

		void PixGraph::compute_pixels_angles_homogenity(){

			all_vertcies_angle_homegenity.reserve(all_vertcies.size());
			for (size_t c_vertex_idx = 0; c_vertex_idx < all_vertcies.size(); ++c_vertex_idx) {
				if (vertecies_N_level_neighbors[c_vertex_idx].first.first.size()<2)
				{
					all_vertcies_angle_homegenity.push_back(double(M_PI));
					continue;
				}

				double c_pix_angle = all_vertcies_angles[c_vertex_idx];
				uchar  c_pix_heat = all_vertcies_heat[c_vertex_idx];

				std::vector<size_t> neighbours_indices = vertecies_N_level_neighbors[c_vertex_idx].first.first;
				std::vector<uchar> neigbours_heat = get_vector_by_indices(all_vertcies_heat, neighbours_indices);
				std::vector<double> neighbours_angles = get_vector_by_indices(all_vertcies_angles, neighbours_indices);

				double sum = std::accumulate(neighbours_angles.begin(), neighbours_angles.end(), 0.0);
				//double mean = fmod(sum / neighbours_angles.size()+1, M_PI);

				std::vector<double> diff(neighbours_angles.size());
				std::transform(neighbours_angles.begin(), neighbours_angles.end(), diff.begin(),
					[c_pix_angle](double x) {return sqrt(abs(sin(x - c_pix_angle) * sin(c_pix_angle - x))); });//[mean](double x) { return x - mean; });
				double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
				double stdev = std::sqrt(sq_sum / neighbours_angles.size());
				all_vertcies_angle_homegenity.push_back(1+stdev);

			}
		}

		void PixGraph::transform_relative_vector_to_matrix(std::vector<double>& respective_vector, matrix& output_matrix) {

			for (size_t vertex_idx = 0; vertex_idx < all_vertcies.size(); ++vertex_idx) {

				Boost_Point_2 c_point = all_vertcies[vertex_idx];
				double c_value = respective_vector[vertex_idx];

				output_matrix.at<double>(size_t(bg::get<0>(c_point)), size_t(bg::get<1>(c_point))) = c_value;

			}
		}

		void PixGraph::transform_relative_vector_to_matrix(std::vector<size_t>& respective_vector, matrix& output_matrix) {

			for (size_t vertex_idx = 0; vertex_idx < all_vertcies.size(); ++vertex_idx) {

				Boost_Point_2 c_point = all_vertcies[vertex_idx];
				int c_value = static_cast<int>(respective_vector[vertex_idx]);

				output_matrix.at<int>(size_t(bg::get<0>(c_point)), size_t(bg::get<1>(c_point))) = c_value;

			}
		}

		void PixGraph::compute_connected_components() {

			connected_components_label = std::vector<size_t>(all_vertcies.size(), 0);
			std::vector<bool> visited(all_vertcies.size(), false);

			std::function<void(size_t, std::set<size_t>&)> dfsUtil = [&](size_t src, std::set<size_t>& s) -> void {
				s.insert(src);
				visited[src] = true;
				if (!vertecies_N_level_neighbors[src].second.empty()) {
					auto neighbours_iter = get_vertex_first_neighbours(src);
					for (auto it = neighbours_iter.first; it != neighbours_iter.second; it++) {
						if (!visited[*it]) {
							dfsUtil(*it, s);
						}
					}
				}
			};

			auto dfs = [&](size_t v) {
				//visited[v] =false;
				size_t c_label = 1;
				tqdm bar;
				for (size_t i = 0; i < v; i++) {
					bar.progress(i, all_vertcies.size());
					if (!visited[i]) {
						std::set<size_t> s;
						dfsUtil(i, s);
						for (auto j = s.begin(); j != s.end(); j++) {
							connected_components_label[*j] = c_label;
						}
						c_label++;
					}
				}
				bar.finish();
			};
			dfs(all_vertcies.size());
		}

		void PixGraph::compute_connected_components_iter() {

			connected_components_label = std::vector<size_t>(all_vertcies.size(), 0);
			std::vector<bool> visited(all_vertcies.size(), false);

			auto dfs = [&](size_t v, size_t c_label) {

				// create a queue for doing BFS
				std::queue<size_t> q;
				// mark the source vertex as discovered
				visited[v] = true;
				// enqueue source vertex
				q.push(v);
				// loop till queue is empty
				while (!q.empty())
				{
					// dequeue front node and print it
					v = q.front();
					q.pop();
					connected_components_label[v] = c_label;

					// do for every edge `v —> u`
					if (!vertecies_N_level_neighbors[v].second.empty())
					{
						//auto neighbours_iter = get_vertex_first_neighbours(v);
						//for (auto it = neighbours_iter.first; it != neighbours_iter.second; it++)
						for (auto it : vertecies_N_level_neighbors[v].first.first)
						{
							if (!visited[it])
							{
								// mark it as discovered and enqueue it
								visited[it] = true;
								q.push(it);
							}
						}
					}
				}
								
			};			

			size_t c_label = 1;
			tqdm bar;
			
			for (int i = 0; i < all_vertcies.size(); i++)
			{
				bar.progress(i, all_vertcies.size());
				if (visited[i] == false)
				{
					// start BFS traversal from vertex `i`
					dfs(i, c_label);
					c_label++;
				}
			}
			bar.finish();

		}

		void PixGraph::compute_pixels_angles_homogenity2() {

			all_vertcies_angle_homegenity.clear();
			all_vertcies_angle_homegenity.reserve(all_vertcies.size());
			for (size_t c_vertex_idx = 0; c_vertex_idx < all_vertcies.size(); ++c_vertex_idx) {

				std::vector<size_t>& neighbours_indices = vertecies_N_level_neighbors[c_vertex_idx].first.first;
				std::vector<short int>& neighbours_levels = vertecies_N_level_neighbors[c_vertex_idx].second;
				if (neighbours_indices.size() < 2)
				{
					all_vertcies_angle_homegenity.push_back(double(M_PI));
					continue;
				}

				
				double c_homogenity = 0.0f;
				int last_end_level_idx = 0;
				int level = 1;
				for (size_t c_start_level_idx : neighbours_levels){
					std::vector<size_t> subvector = { neighbours_indices.begin() + last_end_level_idx, neighbours_indices.begin() + c_start_level_idx };
					for (auto c_comb : get_all_combinations(subvector, 2)){
						size_t pt_0_idx = c_comb[0], pt_1_idx = c_comb[1];
						Boost_Point_2& pt_m = all_vertcies[c_vertex_idx], pt_0 = all_vertcies[pt_0_idx], pt_1 = all_vertcies[pt_1_idx];
						double c_angle = angle3p(pt_m, pt_0, pt_1);
						c_angle = std::min<double>(c_angle, M_PI - c_angle);
						c_homogenity += c_angle * level;
					}

					level += 1;
					last_end_level_idx = c_start_level_idx;
				}
				c_homogenity /= (level - 1);
				all_vertcies_angle_homegenity.push_back(c_homogenity);

			};
		}

		void PixGraph::generate_free_segments(){
			
			vectorized_vertcies = std::vector<bool>(all_vertcies.size(), false);
			size_t vectorized_vertcies_count = 0;
			free_segments.clear();
			free_segments_weight.clear();
			free_segments.reserve(all_vertcies.size() / 10);
			free_segments_weight.reserve(all_vertcies.size() / 10);

			std::set<size_t> recommended_start_points;
			std::cout << "Generating free segments" << std::endl;
			tqdm bar;
			while (vectorized_vertcies_count < all_vertcies.size()) {
				bar.progress(vectorized_vertcies_count, all_vertcies.size());
				//BOOST_LOG_TRIVIAL(debug) << "segment count " << free_segments.size();
				size_t iteration_count = 0;

				// get start point (best homogenity_value)
				size_t start_point_index;
				for (auto c_possible_start : recommended_start_points) {
					if (vectorized_vertcies[c_possible_start] == true)
						recommended_start_points.erase(c_possible_start);
				}
				if (recommended_start_points.empty())
					start_point_index = std::min_element(all_vertcies_angle_homegenity.begin(), all_vertcies_angle_homegenity.end()) - all_vertcies_angle_homegenity.begin();
				else
				{
					start_point_index = *std::min_element(
						recommended_start_points.begin(), recommended_start_points.end(),
						[&](const auto& a, const auto& b) { return all_vertcies_angle_homegenity[a] < all_vertcies_angle_homegenity[b]; });
					recommended_start_points.erase(start_point_index);
				}
				assert(vectorized_vertcies[start_point_index] == false);


				Inexact_Point_2 central_point = transform_B2C_Point(all_vertcies[start_point_index]);
				std::set<size_t> fitted_points = { start_point_index };
				Inexact_Line_2* last_common_line = NULL;//(central_point, refrence_axis_vector);
				bool can_expand = true;
				if (vertecies_N_level_neighbors[start_point_index].second.empty()) {
					can_expand = false;
					last_common_line = &Inexact_Line_2(central_point, refrence_axis_vector);
				}
				
				while (can_expand) {
					iteration_count++;
					std::set<size_t> next_level_neighbours;
					for (size_t fitted_pt : fitted_points) {
						auto c_first_neighbours = get_vertex_first_neighbours(fitted_pt);
						next_level_neighbours.insert(c_first_neighbours.first, c_first_neighbours.second);
					}
					std::set<size_t> new_neighbours; // set difference between next_level_neighbours and fitted_points
					std::set_difference(next_level_neighbours.begin(), next_level_neighbours.end(),
						fitted_points.begin(), fitted_points.end(),
						std::inserter(new_neighbours, new_neighbours.end())
					);
					// remove already vectorized from new_neighbours
					auto new_neighbour_it = new_neighbours.begin();
					while (new_neighbour_it != new_neighbours.end())
					{
						if (vectorized_vertcies[*new_neighbour_it]) {
							new_neighbour_it = new_neighbours.erase(new_neighbour_it);
						}
						else {
							++new_neighbour_it;
						}
					}
					if (new_neighbours.empty()) break; // breaks if no more new neighbours

					// try fitting vector
					std::vector<size_t> all_possible_points_indices(new_neighbours.begin(), new_neighbours.end());
					all_possible_points_indices.insert(all_possible_points_indices.end(), fitted_points.begin(), fitted_points.end());

					std::vector<uchar> all_possible_points_heat = get_vector_by_indices(all_vertcies_heat, all_possible_points_indices);
					std::vector<Boost_Point_2> all_possible_points = get_vector_by_indices(all_vertcies, all_possible_points_indices);
					std::vector<Inexact_Point_2> all_possible_points_cgal;
					container_transform_B2C_Points(all_possible_points, all_possible_points_cgal);


					auto get_centroid = [&](std::vector<Inexact_Point_2> pts, std::vector<uchar> pts_heat) {
						double mx=0, my=0;
						size_t heat_sum = 0;
						for (size_t c_p_idx = 0; c_p_idx < pts.size();++c_p_idx) {
							const auto& c_p = pts[c_p_idx];
							const auto& c_p_h = pts_heat[c_p_idx];
							mx += c_p_h * c_p.x(); my += c_p_h*c_p.y();
							heat_sum += c_p_h;
						}
						return Inexact_Point_2(mx / heat_sum, my / heat_sum);
					};
					//update central point
					if (true) {
						central_point = get_centroid(all_possible_points_cgal, all_possible_points_heat);
					}

					std::pair<Inexact_Point_2, uchar> central_pixel_pair(central_point, all_vertcies_heat[start_point_index]);

					std::vector<short int> level_labels(all_possible_points_indices.size(), 1);
					Inexact_Vector_2 common_vector = get_best_fitting_vector(central_pixel_pair, all_possible_points_cgal, all_possible_points_heat, level_labels);
					Inexact_Line_2 common_line(central_point, common_vector);

					for (size_t already_fitted_point_idx : fitted_points) {
						double fit_err = std::sqrt(CGAL::squared_distance(transform_B2C_Point(all_vertcies[already_fitted_point_idx]), common_line));
						if (fit_err > params->simplification_factor && last_common_line != NULL) // maybe add minimum iteration
						{
							// use last fitted line
							common_line = *last_common_line;
						}
					}

					bool should_break = false;
					for (size_t possible_point_idx : all_possible_points_indices) {
						double fit_err = std::sqrt(CGAL::squared_distance(transform_B2C_Point(all_vertcies[possible_point_idx]), common_line));
						if (fit_err < params->simplification_factor) {
							fitted_points.insert(possible_point_idx);
						}
						else {							
							should_break = true;
						}
					}
					last_common_line = &common_line; // maybe out of scope ptr refrencing
					if (should_break) break;

				}

				// creating segment from fitted line & fitted points
				float segment_weight = 0;
				double min_x = all_vertcies[*fitted_points.begin()].get<0>();
				double max_x = all_vertcies[*fitted_points.begin()].get<0>();
				double min_y = all_vertcies[*fitted_points.begin()].get<1>();
				double max_y = all_vertcies[*fitted_points.begin()].get<1>();
				for (size_t fitted_pt : fitted_points) {
					min_x = std::min<double>(min_x, all_vertcies[fitted_pt].get<0>());
					max_x = std::max<double>(max_x, all_vertcies[fitted_pt].get<0>());
					min_y = std::min<double>(min_y, all_vertcies[fitted_pt].get<1>());
					max_y = std::max<double>(max_y, all_vertcies[fitted_pt].get<1>());
					vectorized_vertcies[fitted_pt] = true;
					vectorized_vertcies_count++;
					all_vertcies_angle_homegenity[fitted_pt] = DBL_MAX;
					segment_weight += all_vertcies_heat[fitted_pt];
				}

				//normalizing segment weight
				segment_weight /= fitted_points.size();

				min_x -= 0.5;
				max_x += 0.5; 
				min_y -= 0.5;
				max_y += 0.5;
				Inexact_Iso_rectangle_2 bbox(Inexact_Point_2(min_x, min_y), Inexact_Point_2(max_x, max_y));

				if (fitted_points.size() <= 1) {
					//BOOST_LOG_TRIVIAL(info) << "fitted points not enough! count= " << fitted_points.size();
					continue;
				}

				//Segment creation
				const auto intersection_result = CGAL::intersection(bbox, *last_common_line);
				if (intersection_result) {
					if (const Inexact_Segment_2* s = boost::get<Inexact_Segment_2>(&*intersection_result)) {
						free_segments.push_back(*s);
						free_segments_weight.push_back(segment_weight);
						free_segments_resp_component.push_back(connected_components_label[*fitted_points.begin()]);
						// set recommended point as closest to segmented line
						std::list<Boost_Value_2> possible_next_start_point;
						Inexact_Point_2 seg_s = s->source(), seg_e = s->target();
						bgi::query(points_tree, bgi::nearest(transform_C2B_Point(seg_s), 4), std::back_inserter(possible_next_start_point));
						bgi::query(points_tree, bgi::nearest(transform_C2B_Point(seg_e), 4), std::back_inserter(possible_next_start_point));
						for (auto c_possible_start : possible_next_start_point) {
							if (vectorized_vertcies[c_possible_start.second]==false)
								recommended_start_points.insert(c_possible_start.second);
						}
					}
					else {
						BOOST_LOG_TRIVIAL(debug) << "Error generating free segment!";
					}
					
				}
			}
			bar.finish();
			BOOST_LOG_TRIVIAL(info) << "Generated free segments count: " << free_segments.size();

		}

		void PixGraph::generate_free_polylines(size_t component_label) {

			vectorized_vertcies = std::vector<bool>(all_vertcies.size(), false);
			size_t vectorized_vertcies_count = 0;
			free_segments.clear();
			free_segments_weight.clear();
			free_segments.reserve(all_vertcies.size() / 10);
			free_segments_weight.reserve(all_vertcies.size() / 10);

			std::cout << "Generating free polylines" << std::endl;
			tqdm bar;
			while (vectorized_vertcies_count < all_vertcies.size()) {
				bar.progress(vectorized_vertcies_count, all_vertcies.size());
				size_t iteration_count = 0;


			}
			bar.finish();
			BOOST_LOG_TRIVIAL(info) << "Generated free polylines count: " << free_segments.size();

		}

		void PixGraph::write_free_segments_shapefile(const std::string& output_filename, OGRSpatialReference* source_srs)
		{
			if (free_segments.empty()) {
				std::cout << "Warning : empty list of segments. No output written." << std::endl;
				return;
			}

			GDALDataset* source_dataset = NULL;
			GDALDataset* target_dataset = NULL;

			try {
				const std::string driver_name = "ESRI Shapefile";

				GDALDriver* driver = GetGDALDriverManager()->GetDriverByName(driver_name.c_str());
				if (driver == NULL) {
					throw std::logic_error("Error : ESRI Shapefile driver not available.");
				}

				// Step 2.
				// Writes target file
				target_dataset = driver->Create(output_filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
				if (target_dataset == NULL) {
					throw std::logic_error("Error : creation of output file failed.");
				}

				OGRLayer* target_layer = target_dataset->CreateLayer("lines", source_srs, wkbLineString, NULL);
				if (target_layer == NULL) {
					throw std::logic_error("Error : layer creation failed.");
				}

				OGRFieldDefn o_field_id("ID", OFTInteger);

				if (target_layer->CreateField(&o_field_id) != OGRERR_NONE) {
					throw std::logic_error("Error : field creation failed.");
				}

				OGRFieldDefn o_field_segment_weight("seg_w", OFTReal);

				if (target_layer->CreateField(&o_field_segment_weight) != OGRERR_NONE) {
					throw std::logic_error("Error : field creation failed.");
				}

				OGRFieldDefn o_field_cluster_label("cluster", OFTInteger);

				if (target_layer->CreateField(&o_field_cluster_label) != OGRERR_NONE) {
					throw std::logic_error("Error : field creation failed.");
				}

				for (size_t i = 0; i < free_segments.size(); ++i) {
					Inexact_Segment_2 S = free_segments[i];
					float c_segment_weight = free_segments_weight[i];
					size_t c_segment_cluster = free_segments_resp_component[i];
					OGRLineString ogr_linestring;

					const Inexact_Point_2& S1 = S.source();
					double dfX =
						geotransform[0]
						+ S1.y() * geotransform[1]
						+ S1.x() * geotransform[2];
					double dfY =
						geotransform[3]
						+ S1.y() * geotransform[4]
						+ S1.x() * geotransform[5];
					//double s1_refrenced_x = origin_point.get<0>() + S1.x() * abs(pix_width);
					//double s1_refrenced_y = origin_point.get<1>() + S1.y() * abs(pix_height);
					ogr_linestring.addPoint(&OGRPoint(dfX, dfY));
					const Inexact_Point_2& S2 = S.target();
					dfX =
						geotransform[0]
						+ S2.y() * geotransform[1]
						+ S2.x() * geotransform[2];
					dfY =
						geotransform[3]
						+ S2.y() * geotransform[4]
						+ S2.x() * geotransform[5];
					//double s2_refrenced_x = origin_point.get<0>() + S2.x() * abs(pix_width);
					//double s2_refrenced_y = origin_point.get<1>() + S2.y() * abs(pix_height);
					ogr_linestring.addPoint(&OGRPoint(dfX, dfY));

					OGRFeature* feature;
					feature = OGRFeature::CreateFeature(target_layer->GetLayerDefn());

					feature->SetGeometry(&ogr_linestring);
					feature->SetField("ID", int(i));
					feature->SetField("seg_w", float(c_segment_weight));
					feature->SetField("cluster", int(c_segment_cluster));

					// Writes new feature
					OGRErr error = target_layer->CreateFeature(feature);
					if (error != OGRERR_NONE) std::cout << "Error code : " << int(error) << std::endl;
					OGRFeature::DestroyFeature(feature);
				}

			}
			catch (std::exception& e) {
				std::cout << e.what() << std::endl;
			}

			if (source_dataset != NULL) GDALClose(source_dataset);
			if (target_dataset != NULL) GDALClose(target_dataset);
		}

	
	}
}