#pragma once
#include "defs.h"
#include "orientations/common_orientations.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/make_shared.hpp>
#include <boost/range/iterator_range.hpp>
#include "geometries_with_attributes/linestring_with_attributes.h"
#include "design_pattern/insertion_order_set.h"

namespace LxGeo
{

	using namespace GeometryFactoryShared;

	namespace faultToVector
	{

		struct PixelData {
			int row;
			int col;
			double val;
			bool vectorized;
		};
		const orientationPair getPixelsOrientation(PixelData& p1, PixelData& p2) {
			return getOrientation(p1.col, p1.row, p2.col, p2.row);
		};

		typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::no_property, boost::no_property> SkeletonPixGraph;

		typedef boost::graph_traits<SkeletonPixGraph>::vertex_iterator vertex_iterator;
		typedef boost::graph_traits<SkeletonPixGraph>::vertex_descriptor vertex_descriptor;
		typedef boost::graph_traits<SkeletonPixGraph>::edge_descriptor edge_descriptor;
		typedef boost::graph_traits<SkeletonPixGraph>::adjacency_iterator adjacency_iterator;
		typedef boost::graph_traits<SkeletonPixGraph>::out_edge_iterator out_edge_iterator;

		typedef boost::filtered_graph<SkeletonPixGraph, std::function<bool(edge_descriptor)>, std::function<bool(vertex_descriptor)> > ComponentGraph;
		typedef boost::shared_ptr<std::vector<size_t>> vertex_component_map;

		enum vertexDegreeType {
			UNIQUE = 0,
			LEAF = 1,
			MIDDLE =2,
			JUNCTION = 3,
		};

		struct simplifiedPt {
			int x;
			int y;
			bool ignored;
			simplifiedPt(int _x, int _y, bool _ignored) :x(_x), y(_y), ignored(_ignored) {};
		};

		// lambda used to transform a container of simplified points 
		auto connect_simplified_pts = [](std::vector<simplifiedPt>& simplified_pts)->std::list<Boost_Point_2> {

			auto start_pt = *simplified_pts.begin();
			std::list<Boost_Point_2> simplified_line; simplified_line.push_back(Boost_Point_2(start_pt.x, start_pt.y));
			for (auto c_pt_iter = simplified_pts.begin() + 1; c_pt_iter != simplified_pts.end(); ++c_pt_iter) {
				if (c_pt_iter->ignored) continue;
				else {
					auto next_pt_iter = c_pt_iter + 1;
					if (next_pt_iter == simplified_pts.end())
						simplified_line.push_back(Boost_Point_2(c_pt_iter->x, c_pt_iter->y));
					else
						simplified_line.push_back(Boost_Point_2(
							double(next_pt_iter->x - c_pt_iter->x) / 2 + c_pt_iter->x,
							double(next_pt_iter->y - c_pt_iter->y) / 2 + c_pt_iter->y
						));
				}
			}
			return simplified_line;
		};

		auto getDegreeType = [](int degree_cnt)->vertexDegreeType { return vertexDegreeType(std::min<int>(degree_cnt, vertexDegreeType::JUNCTION)); };

		class SkeletonTracer
		{
		public:

			//SkeletonTracer() {};

			SkeletonTracer(matrix& _image_matrix, double _min_heat_value) :
				image_matrix(_image_matrix), min_heat_value(_min_heat_value)
			{
				pixels_count = static_cast<size_t>(cv::countNonZero(_image_matrix));
				sp_graph = SkeletonPixGraph(pixels_count);
				init_vertex_data_map();
				init_edges();
			};

			void init_vertex_data_map() {
				
				points_tree = Boost_RTree_2_linear_points();

				size_t c_vertex_idx = 0;
				for (int c_row = 0; c_row < image_matrix.rows; ++c_row) {
					for (int c_col = 0; c_col < image_matrix.cols; ++c_col) {

						double c_pixel_val = static_cast<double>(image_matrix.at<uchar>(c_row, c_col));
						if (c_pixel_val <= min_heat_value) continue;

						PixelData c_pixel_data{ c_row, c_col, c_pixel_val, false };
						pixel_data_map[c_vertex_idx] = c_pixel_data;

						points_tree.insert(Boost_Value_2_point(Boost_Point_2(c_row, c_col), c_vertex_idx));
						c_vertex_idx++;

					}
				}

			}

			void init_edges() {

				struct HaveConnection
				{
					PixelData const& _ref;
					std::map<size_t, PixelData>& _pixel_data_map;

					HaveConnection(PixelData const& ref, std::map<size_t, PixelData>& pixel_data_map)
						: _ref(ref), _pixel_data_map(pixel_data_map)
					{}

					bool operator()(Boost_Value_2_point const& tar_val) const
					{
						PixelData& tar = _pixel_data_map[tar_val.second];
						return (std::abs<int>(_ref.row - tar.row) <= 1)
							&& (std::abs<int>(_ref.col - tar.col) <= 1);
					}
				};

				for (size_t c_vertex_idx = 0; c_vertex_idx < pixels_count; ++c_vertex_idx) {

					PixelData& c_pixel_data = pixel_data_map[c_vertex_idx];
					Boost_Box_2 box(Boost_Point_2(c_pixel_data.row -1, c_pixel_data.col - 1), Boost_Point_2(c_pixel_data.row + 1, c_pixel_data.col +1));

					std::list<Boost_Value_2_point> possible_neighbors;
					points_tree.query(bgi::intersects(box) && bgi::satisfies(HaveConnection(c_pixel_data, pixel_data_map)),
						std::back_inserter(possible_neighbors));

					std::set<orientationPair> neighbours_orientations;
					std::transform(possible_neighbors.begin(), possible_neighbors.end(), std::inserter(neighbours_orientations, neighbours_orientations.begin()),
						[&](auto& c_neigh_value)->const orientationPair { return getPixelsOrientation(c_pixel_data, pixel_data_map[c_neigh_value.second]); }
					);

					possible_neighbors.erase(
						std::remove_if(possible_neighbors.begin(), possible_neighbors.end(),
							[&](const Boost_Value_2_point& o) {
								auto c_orientation = getPixelsOrientation(c_pixel_data, pixel_data_map[o.second]);
								if (c_orientation.symbol == orientationSymbol::C) return true;
								if (c_orientation.main_orientation()) return false;
								for (auto related_main_orientation_symbol : secondary_main_orientations_relation.at(c_orientation.symbol)) {
									auto related_main_orientation = orientations_map.at(related_main_orientation_symbol);
									if (neighbours_orientations.find(related_main_orientation) != neighbours_orientations.end())
										return true;
								}
								return false;
							}),
						possible_neighbors.end());

					for (auto& c_possible_neigh : possible_neighbors) {
						auto& respective_idx = c_possible_neigh.second;
						auto& rescpective_pixel_data = pixel_data_map[respective_idx];
						if (std::abs<int>(c_pixel_data.col - rescpective_pixel_data.col) == 1
							|| std::abs<int>(c_pixel_data.row - rescpective_pixel_data.row) == 1
							) {
							if (boost::edge(respective_idx, c_vertex_idx, sp_graph).second) continue;
							boost::add_edge(c_vertex_idx, respective_idx, sp_graph);
						}
					}

				}

			}

			std::vector<LineString_with_attributes> export_edge_graph_as_LSwithAttr() {
				std::vector<LineString_with_attributes> edges_linestrings;
				edges_linestrings.reserve( boost::num_edges(sp_graph) );
				auto es = boost::edges(sp_graph);
				for (auto eit = es.first; eit != es.second; ++eit) {
					vertex_descriptor source_ = boost::source(*eit, sp_graph), target_ = boost::target(*eit, sp_graph);
					PixelData& source_data = pixel_data_map[source_], target_data = pixel_data_map[target_];
					Boost_Point_2 source_rep_pt(source_data.row, source_data.col), target_rep_pt(target_data.row, target_data.col);

					LineString_with_attributes c_edge_container(Boost_LineString_2({ source_rep_pt, target_rep_pt }));
					c_edge_container.set_string_attribute("source", std::to_string(size_t(source_)));
					c_edge_container.set_string_attribute("target", std::to_string(size_t(target_)));
					edges_linestrings.push_back(c_edge_container);
				}
				return edges_linestrings;
			}

			std::vector<ComponentGraph> connected_components_subgraphs()
			{
				vertex_component_map mapping = boost::make_shared<std::vector<size_t>>(boost::num_vertices(sp_graph));
				size_t num = boost::connected_components(sp_graph, mapping->data());

				std::vector<ComponentGraph> component_graphs;

				for (size_t i = 0; i < num; i++)
					component_graphs.emplace_back(sp_graph,
						[mapping, i, this](edge_descriptor e) {
							return mapping->at(boost::source(e, sp_graph)) == i
								|| mapping->at(boost::target(e, sp_graph)) == i;
						},
						[mapping, i](vertex_descriptor v) {
							return mapping->at(v) == i;
						});

				return component_graphs;
			}

			struct attachedLineString {
				int start_extremity_id;
				int end_extremity_id;
				Boost_LineString_2 path_linestring;
				attachedLineString(int se, int ee, std::list<Boost_Point_2>& pp) : start_extremity_id(se), end_extremity_id(ee), path_linestring(pp.begin(), pp.end()) {};
			};

			std::list<attachedLineString> extract_junction_lines(ComponentGraph& component_graph, int max_orientations_count = 2, int max_succesive_count = 4) {
				
				std::list<attachedLineString> extracted_lines;
				
				ComponentGraph::vertex_descriptor start_pt=NULL;
				ComponentGraph::vertex_iterator v, vend;
				for (boost::tie(v, vend) = boost::vertices(component_graph); v != vend; ++v) {
					int c_vertex_degree = boost::degree(*v, component_graph);
					if (getDegreeType(c_vertex_degree) == vertexDegreeType::LEAF || getDegreeType(c_vertex_degree) == vertexDegreeType::JUNCTION) {
						start_pt = *v;
						break;
					}
				}
				if (boost::degree(start_pt, component_graph) == 0)
					return extracted_lines;
							

				// algo related datastructres
				std::set<orientationPair> traced_orientations;
				orientationPair last_orientation = orientations_map.at(orientationSymbol::C);
				int stability_count = 0;
				bool is_straight = false;
				int start_extremity = NULL, end_extremity = NULL;

				// DFS lambda for tracing
				std::function<void(ComponentGraph::vertex_descriptor, ComponentGraph::vertex_descriptor, std::vector<simplifiedPt>&)> dfs 
					= [&](ComponentGraph::vertex_descriptor v0, ComponentGraph::vertex_descriptor v1, std::vector<simplifiedPt>& simplified_pts) -> void {
					
					vertexDegreeType target_degree_type = getDegreeType(boost::degree(v1, component_graph));
					// boolean representing if v1 is simplified
					bool simplified_pt = !(target_degree_type == vertexDegreeType::JUNCTION || target_degree_type == vertexDegreeType::LEAF); //true

					PixelData& v0_pd = pixel_data_map[v0]; PixelData& v1_pd = pixel_data_map[v1];
					assert(!v1_pd.vectorized); 
					if (target_degree_type!=vertexDegreeType::JUNCTION) v1_pd.vectorized = true;
					orientationPair c_orientation = getPixelsOrientation(v0_pd, v1_pd);

					if (traced_orientations.find(c_orientation) == traced_orientations.end()) {// a new orientation is intorduced
						if (traced_orientations.size() == max_orientations_count) {//orientation break
							(simplified_pts.end()-1)->ignored = false;
							//clear respective structres
							traced_orientations.clear(); //traced_orientations.insert(c_orientation);
							stability_count = 0;
						}
						else {
							if (is_straight) {//straight line break
								(simplified_pts.end() - 1)->ignored = false;
								traced_orientations.clear(); //traced_orientations.insert(c_orientation);
								stability_count = 0;
								is_straight = false;
							}
							else {// add new orientation
								traced_orientations.insert(c_orientation);
							}
						}
					}
					else {// old orientation encountred
						stability_count += (c_orientation.symbol == last_orientation.symbol);
						if (stability_count >= max_succesive_count) { // max consecutive orientation check
							if (traced_orientations.size() <= 1) {// straight line case
								is_straight = true;
							}
							else {// succesive point break
								(simplified_pts.end() - max_succesive_count)->ignored = false; // remove simplification of N-succession point
								is_straight = true;
								traced_orientations.erase(c_orientation);
								stability_count = 0;
							}
						}
					}

					simplified_pts.emplace_back(v1_pd.row, v1_pd.col, simplified_pt);
					last_orientation = c_orientation;
					
					if (target_degree_type == vertexDegreeType::JUNCTION) {
						// add traced lines with respective extremeties
						end_extremity = v1;
						extracted_lines.emplace_back(start_extremity, end_extremity, connect_simplified_pts(simplified_pts));
						ComponentGraph::adjacency_iterator v, vend;
						for (boost::tie(v, vend) = adjacent_vertices(v1, component_graph); v != vend; ++v) {
							if (*v != v0 && !pixel_data_map[*v].vectorized) {
								start_extremity = v1;
								std::vector<simplifiedPt> c_trace_pts; c_trace_pts.emplace_back(pixel_data_map[start_extremity].row, pixel_data_map[start_extremity].col, false);
								dfs(v1, *v, c_trace_pts);
							}
						}
					}
					if (target_degree_type == vertexDegreeType::LEAF) {
						extracted_lines.emplace_back(start_extremity, end_extremity, connect_simplified_pts(simplified_pts));
					}
					if (target_degree_type == vertexDegreeType::MIDDLE) {
						// get correct next neighbour
						ComponentGraph::vertex_descriptor v2;
						ComponentGraph::adjacency_iterator v, vend;
						for (boost::tie(v, vend) = adjacent_vertices(v1, component_graph); v != vend; ++v) {
							if (*v!=v0) {
								v2 = *v;
								break;
							}
						}
						// continue tracing
						dfs(v1, v2, simplified_pts);
					}

				};

				ComponentGraph::adjacency_iterator vadj, vadj_end;
				for (boost::tie(vadj, vadj_end) = adjacent_vertices(start_pt, component_graph); vadj != vadj_end; ++vadj) {
					ComponentGraph::vertex_descriptor vadj_desc = *vadj;
					if (!pixel_data_map[vadj_desc].vectorized) {
						start_extremity = start_pt;
						std::vector<simplifiedPt> c_trace_pts; c_trace_pts.emplace_back(pixel_data_map[start_extremity].row, pixel_data_map[start_extremity].col, false);
						dfs(start_pt, *vadj, c_trace_pts);
					}
				}

				return extracted_lines;
			}

			std::list<attachedLineString> extract_junction_lines_iterative(ComponentGraph& component_graph, int max_orientations_count = 2, int max_succesive_count = 4) {

				std::list<attachedLineString> extracted_lines;

				insertion_order_set<ComponentGraph::vertex_descriptor> start_pts;
				ComponentGraph::vertex_descriptor sample_leaf=-1;
				ComponentGraph::vertex_iterator v, vend;
				for (boost::tie(v, vend) = boost::vertices(component_graph); v != vend; ++v) {
					int c_vertex_degree = boost::degree(*v, component_graph);
					if (getDegreeType(c_vertex_degree) == vertexDegreeType::LEAF | getDegreeType(c_vertex_degree) == vertexDegreeType::JUNCTION) {
						start_pts.push_back(*v);
						break;
					}
				}												
				
				int start_extremity = -1, end_extremity = -1;

				std::function<void(ComponentGraph::vertex_descriptor, ComponentGraph::vertex_descriptor, std::vector<simplifiedPt>&)> dfs_iter
					= [&](ComponentGraph::vertex_descriptor v0, ComponentGraph::vertex_descriptor v1, std::vector<simplifiedPt>& simplified_pts) {

					// algo related datastructres
					std::set<orientationPair> traced_orientations;
					orientationPair last_orientation = orientations_map.at(orientationSymbol::C);
					int stability_count = 0;
					bool is_straight = false;
					//self breakable loop
					while (true) {
						vertexDegreeType target_degree_type = getDegreeType(boost::degree(v1, component_graph));
						// boolean representing if v1 is simplified
						bool simplified_pt = (target_degree_type == vertexDegreeType::MIDDLE);
						PixelData& v0_pd = pixel_data_map[v0]; PixelData& v1_pd = pixel_data_map[v1];
						v1_pd.vectorized = true;
						orientationPair c_orientation = getPixelsOrientation(v0_pd, v1_pd);

						if (traced_orientations.find(c_orientation) == traced_orientations.end()) {// a new orientation is intorduced
							if (traced_orientations.size() == max_orientations_count) {//orientation break
								(simplified_pts.end() - 1)->ignored = false;
								//clear respective structres
								traced_orientations.clear(); //traced_orientations.insert(c_orientation);
								stability_count = 0;
							}
							else {
								if (is_straight) {//straight line break
									(simplified_pts.end() - 1)->ignored = false;
									traced_orientations.clear(); //traced_orientations.insert(c_orientation);
									stability_count = 0;
									is_straight = false;
								}
								else {// add new orientation
									traced_orientations.insert(c_orientation);
								}
							}
						}
						else {// old orientation encountred
							stability_count += (c_orientation.symbol == last_orientation.symbol);
							if (stability_count >= max_succesive_count) { // max consecutive orientation check
								if (traced_orientations.size() <= 1) {// straight line case
									is_straight = true;
								}
								else {// succesive point break
									(simplified_pts.end() - max_succesive_count)->ignored = false; // remove simplification of N-succession point
									is_straight = true;
									traced_orientations.erase(c_orientation);
									stability_count = 0;
								}
							}
						}

						simplified_pts.emplace_back(v1_pd.col, v1_pd.row, simplified_pt);

						last_orientation = c_orientation;
						if (target_degree_type == vertexDegreeType::JUNCTION) {
							// add traced lines with respective extremeties
							end_extremity = v1;
							extracted_lines.emplace_back(start_extremity, end_extremity, connect_simplified_pts(simplified_pts));
							start_pts.push_back(end_extremity);
							break;
						}
						if (target_degree_type == vertexDegreeType::LEAF) {
							end_extremity = v1;
							extracted_lines.emplace_back(start_extremity, end_extremity, connect_simplified_pts(simplified_pts));
							break;
						}
						if (target_degree_type == vertexDegreeType::MIDDLE) {
							// get correct next neighbour
							ComponentGraph::vertex_descriptor v2;
							ComponentGraph::adjacency_iterator v, vend;
							for (boost::tie(v, vend) = adjacent_vertices(v1, component_graph); v != vend; ++v) {
								if (*v != v0) {
									v2 = *v;
									break;
								}
							}
							// continue tracing
							v0 = v1; v1 = v2;
						}
					
					}

				};

				// trace all start points
				auto start_pt_iterator = start_pts.begin();
				while(start_pt_iterator != start_pts.end()) {
					auto c_start_pt = *start_pt_iterator;
					ComponentGraph::adjacency_iterator vadj, vadj_end;
					for (boost::tie(vadj, vadj_end) = adjacent_vertices(c_start_pt, component_graph); vadj != vadj_end; ++vadj) {
						ComponentGraph::vertex_descriptor vadj_desc = *vadj;
						if (pixel_data_map[vadj_desc].vectorized) continue;
						start_extremity = c_start_pt;
						std::vector<simplifiedPt> simplified_pts; simplified_pts.emplace_back(pixel_data_map[c_start_pt].col, pixel_data_map[c_start_pt].row, false);
						dfs_iter(c_start_pt, *vadj, simplified_pts);						
					}
					start_pt_iterator = std::next(start_pt_iterator);
				}
				return extracted_lines;

			}

			std::vector<LineString_with_attributes> export_attachedLineString_as_LSwithAttr(std::list<attachedLineString>& input_attached_lines,
				size_t comp_id) {
				std::vector<LineString_with_attributes> out_linestrings;
				out_linestrings.reserve(input_attached_lines.size());

				for (attachedLineString& c_attached_line : input_attached_lines) {
					LineString_with_attributes c_edge_container(c_attached_line.path_linestring);
					c_edge_container.set_string_attribute("st_ext", std::to_string(c_attached_line.start_extremity_id));
					c_edge_container.set_string_attribute("end_ext", std::to_string(c_attached_line.end_extremity_id));
					c_edge_container.set_int_attribute("comp", static_cast<int>(comp_id));
					out_linestrings.push_back(c_edge_container);
				}
				return out_linestrings;

			};

			std::vector<LineString_with_attributes> export_component_edge_graph_as_LSwithAttr(ComponentGraph& comp_graph) {
				std::vector<LineString_with_attributes> edges_linestrings;
				edges_linestrings.reserve(boost::num_edges(comp_graph));
				auto es = boost::edges(comp_graph);
				for (auto eit = es.first; eit != es.second; ++eit) {
					ComponentGraph::vertex_descriptor source_ = boost::source(*eit, sp_graph), target_ = boost::target(*eit, sp_graph);
					PixelData& source_data = pixel_data_map[source_], target_data = pixel_data_map[target_];
					Boost_Point_2 source_rep_pt(source_data.row, source_data.col), target_rep_pt(target_data.row, target_data.col);

					LineString_with_attributes c_edge_container(Boost_LineString_2({ source_rep_pt, target_rep_pt }));
					c_edge_container.set_string_attribute("source", std::to_string(size_t(source_)));
					c_edge_container.set_string_attribute("target", std::to_string(size_t(target_)));
					edges_linestrings.push_back(c_edge_container);
				}
				return edges_linestrings;
			}

			~SkeletonTracer() {};

		public:			
			// input matrix to vectorize
			matrix image_matrix;
			double min_heat_value;
		private:
			size_t pixels_count;
			SkeletonPixGraph sp_graph;
			Boost_RTree_2_linear_points points_tree;
			std::map<size_t, PixelData> pixel_data_map;

		};


	}
}