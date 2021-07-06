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

	}
}