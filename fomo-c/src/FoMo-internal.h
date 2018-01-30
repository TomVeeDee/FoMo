#include <vector>
#include <string>

namespace FoMo
{
	double readgoftfromchianti(const std::string chiantifile);
	DataCube readgoftfromchianti(const std::string chiantifile, std::string & ion, double & lambda0, double & atweight);
	GoftCube emissionfromdatacube(DataCube, std::string, std::string, const FoMoObservationType);
	
#ifdef HAVE_CGAL_DELAUNAY_TRIANGULATION_2_H	
	FoMo::RenderCube RenderWithCGAL(FoMo::DataCube datacube, FoMo::GoftCube goftcube, FoMoObservationType observationtype, 
	const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec, const std::string outfile);
	
	FoMo::RenderCube RenderWithCGAL2D(FoMo::DataCube datacube, FoMo::GoftCube goftcube, FoMoObservationType observationtype, 
	const int x_pixel, const int y_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, const std::string outfile);
#endif
	
	FoMo::RenderCube RenderWithNearestNeighbour(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec, const std::string outfile);
	
	FoMo::RenderCube RenderWithThomson(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel,
	std::vector<double> lvec, std::vector<double> bvec, const std::string outfile);
	
	FoMo::RenderCube RenderWithProjection(FoMo::GoftCube goftcube, const int x_pixel, const int y_pixel, const int z_pixel, const int lambda_pixel, const double lambda_width,
	std::vector<double> lvec, std::vector<double> bvec, const std::string outfile);
}
