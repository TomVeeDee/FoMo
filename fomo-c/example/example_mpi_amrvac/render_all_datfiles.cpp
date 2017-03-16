#include "FoMo.h"
#include "FoMo-amrvac.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/replace.hpp>

namespace po = boost::program_options;
using namespace std;

// next two routines copied from http://stackoverflow.com/questions/3300419/file-name-matching-with-wildcard
void EscapeRegex(string &regex)
{
    boost::replace_all(regex, "\\", "\\\\");
    boost::replace_all(regex, "^", "\\^");
    boost::replace_all(regex, ".", "\\.");
    boost::replace_all(regex, "$", "\\$");
    boost::replace_all(regex, "|", "\\|");
    boost::replace_all(regex, "(", "\\(");
    boost::replace_all(regex, ")", "\\)");
    boost::replace_all(regex, "[", "\\[");
    boost::replace_all(regex, "]", "\\]");
    boost::replace_all(regex, "*", "\\*");
    boost::replace_all(regex, "+", "\\+");
    boost::replace_all(regex, "?", "\\?");
    boost::replace_all(regex, "/", "\\/");
}

bool MatchTextWithWildcards(const string &text, string wildcardPattern, bool caseSensitive = true)
{
    // Escape all regex special chars
    EscapeRegex(wildcardPattern);

    // Convert chars '*?' back to their regex equivalents
    boost::replace_all(wildcardPattern, "\\?", ".");
    boost::replace_all(wildcardPattern, "\\*", ".*");

    boost::regex pattern(wildcardPattern, caseSensitive ? boost::regex::normal : boost::regex::icase);

    return regex_match(text, pattern);
}

int main(int argc, char* argv[])
{
	string amrvac_version, compstring, parstring;
	int gamma_eqparposition;
	double n_unit, Teunit, L_unit;
	
	double pi=4*atan(1.);
	vector<double> bangles,langles;

	// parse options to the program
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("file,f", po::value<string>(&compstring)->default_value("\\*.dat"),"search pattern for .dat files (e.g. directory/prefix\\*.dat")
		("par,p", po::value<string>(&parstring)->default_value("amrvac.par"),"name of amrvac.par file")
		("version,V", po::value<string>(&amrvac_version)->default_value("gitlab"), "set amrvac version (gitlab or old)")
		("gamma,g", po::value<int>(&gamma_eqparposition)->default_value(0), "set position of gamma in eqpar vector")
		("n_unit,n", po::value<double>(&n_unit)->default_value(1.),"set value for de-normalisation of density to number density (code units to cm^-3)")
		("T_unit,T", po::value<double>(&Teunit)->default_value(1.),"set value for de-normalisation of temperature (code units to K)")
		("L_unit,L", po::value<double>(&L_unit)->default_value(1.),"set value for de-normalisation of distances (code units to Mm)")
		("langles,l", po::value<vector<double>>(&langles)->multitoken(),"set l angles to rotate around z-axis (radians), multiple values possible")
		("bangles,b", po::value<vector<double>>(&bangles)->multitoken(),"set b angles to rotate around y-axis (radians), multiple values possible")
		;
		
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	
	if (vm.count("help")) {
		cout << desc << "\n";
		return 1;
	}

	if (bangles.size()==0 || langles.size()==0)
	{
		if (langles.size()==0) cout << "No l angles specified! No rendering will be performed." << endl;
		if (bangles.size()==0) cout << "No b angles specified! No rendering will be performed." << endl;
	}
	
	// search for files that contain string compstring
	// first find the path of the files
	// split compstring into path + "/" + regex
	string pathstring;
	size_t lastslash = compstring.find_last_of("/");
	if (lastslash == string::npos)
	{
		// there is no slash in the compstring
		pathstring=".";
		compstring="./"+compstring;
	}
	else
	{
		pathstring = compstring.substr(0,lastslash);
	}
	
	boost::filesystem::directory_iterator end_itr; // default construction yields past-the-end
	vector<string> filelist;
	for ( boost::filesystem::directory_iterator itr( pathstring ); itr != end_itr; ++itr )
	{
		if ( boost::filesystem::is_regular_file(itr->status()) && MatchTextWithWildcards(itr->path().string(),compstring) )
		{
			string filename=itr->path().string();
			filelist.push_back( filename );
		}
	}
	int nframes=filelist.size();
	if (nframes==0) cout << "No files found." << endl;
	
	// Initialize the FoMo object
	FoMo::FoMoObject Object;
	
	for (int t=0; t<nframes; t++)
	{
		string filename=filelist[t];
		cout << "Doing file " << t+1 << " of " << nframes << ": read from " << filename << endl << flush;
		
		Object = read_amrvac_dat_file(filename.c_str(), parstring.c_str(), amrvac_version, gamma_eqparposition, n_unit, Teunit, L_unit);
		
		// data is in structure, now start the rendering
		int x_pixel=149;
		int y_pixel=64;
		int z_pixel=500;
		int lambda_pixel=100;
		double lambda_width=200000;
		Object.setresolution(x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width);
		Object.setrendermethod("NearestNeighbour");
		Object.setobservationtype(FoMo::Spectroscopic);
		Object.setchiantifile("../chiantitables/goft_table_fe_12_0194_abco.dat");
		stringstream ss;
		ss << filename << ".fomo";
		Object.setoutfile(ss.str());
		Object.render(langles,bangles);
	}
}
