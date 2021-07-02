#include "gdal.h"
#include "parameters.h"
#include "vectorization.h"


using namespace LxGeo::faultToVector;

int main(int argc, char* argv[])
{
	clock_t t_begin = clock();
	GDALAllRegister();

	// Reads command-line parameters

	params = new Parameters(argc, argv);
	if (!params->initialized()) {
		delete params;
		return 1;
	}

	// Runs process

	Vectorization* V = new Vectorization();
	if (V->pre_check())
		V->run();

	// Quits

	delete V;
	delete params;

	clock_t t_end = clock();
	//std::cout << "** Elapsed time : " << double(t_end - t_begin) / CLOCKS_PER_SEC << " s." << std::endl;

	return 0;
}