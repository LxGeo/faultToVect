#pragma once
#include "defs.h"
#include "numcpp/line_space.h"

namespace LxGeo
{
	namespace faultToVector
	{

		/**
		*  A Vectorization class to manage running required steps to generate final vectorization shapefile.
		*/
		class Vectorization
		{

		public:
			Vectorization() {
				apply_thin = true;
			};

			~Vectorization() {};

			/**
			*  A method used to run all steps of vectorization.
			*/
			virtual void run();

			/**
			*  A method to check requirements before running vectorization steps.
			* Example: -Checking input_raster exsitance, checking output_path overwrite, check algorithm parameters ...
			* @return an bool indicating if can run algorithm
			*/
			bool pre_check();

			bool apply_thin;

		};
	}
}