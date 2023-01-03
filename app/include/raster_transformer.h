#pragma once
#include "defs.h"
#include "io_raster.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace faultToVector
	{
		/*
		Apply thinning on grayscale raster, by computing the binary raster -> apply thinning on binary -> mask grayscale by thinned binary image
		*/
		void thin_raster(RasterIO& raster_to_thin, RasterIO& thinned_raster, int min_value = 0);

		/*
		* Apply thinning on different binary thresholded probability map -> aggregate results
		*/
		// {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.92,0.95,0.99}
		matrix multi_threshold_thinning(matrix& in_probability_matrix, std::list<double> threshold_list, std::list<double> thresh_weigths = {});

		/*
		* Invert small holes values, either background or foreground
		*/
		void fill_small_components(matrix& in_binary, int min_size, int connectivity);


	}
}