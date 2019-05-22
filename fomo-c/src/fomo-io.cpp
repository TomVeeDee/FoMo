#include "FoMo.h"
#include "FoMo-internal.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

// escape double quotes from the shell:
// https://www.daniweb.com/programming/software-development/threads/348802/passing-string-as-d-compiler-option
// https://bytes.com/topic/c/answers/443326-portable-way-pass-string-c-macro
#define XSTR(x) #x
#define STR(x) XSTR(x)
#define FOMO_VER STR(FOMOVERSION)

/**
 * @brief This writes out the contents of the GoftCube.
 * 
 * This member function writes out the contents of the GoftCube to the hard disk, in the filename given by the argument. 
 * The file can then be read into IDL with the commands provide under the idl/ directory, here in particular readgoftcube.pro. The python/ directory contains routines to read the data into python. 
 * Normally, one would use it to post-process the forward modelling results (when using with a RenderCube), but it is also very useful for debugging purposes when using it with a GoftCube before the rendering.
 * @param filename This parameter specifies which filename the data needs to be written to.
 */
void writegoftcube_txt(FoMo::GoftCube goftcube, const std::string filename)
{
        // write out goftcube to file "filename"
        int commrank;
        std::string space=" ";
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
        commrank = 0;
#endif
        if (commrank==0)
        {
                std::ofstream out(filename,std::ios::ate);
                if (out.is_open())
                {
                        int nvars = goftcube.readnvars();
                        int dim = goftcube.readdim();
                        int ng = goftcube.readngrid();
						out << "FoMo-" << FOMO_VER << std::endl;
                        out << dim << std::endl;
                        out << ng << std::endl;
                        out << nvars << std::endl;
						std::vector<std::string> unitvec=goftcube.readunit();
						for (int u=0; u<dim+nvars; u++)
						{
							out << unitvec.at(u) << "\t";
						}
						out << std::endl;
                        out << goftcube.readchiantifile() << std::endl;
                        out << goftcube.readabundfile() << std::endl;
						FoMo::tgrid grid=goftcube.readgrid();
						FoMo::tvars vars;
						for (int i=0; i<nvars; i++)
						{
							vars.push_back(goftcube.readvar(i));
						}
						out << std::setprecision(8);
                        for (int j=0; j<ng; j++)
                        {
                                for (int i=0; i<dim; i++)
                                {
                                        out << grid[i][j] << space;
                                }
                                for (int i=0; i<nvars; i++)
                                {
                                        out << vars[i][j] << space;
                                }
                                out << std::endl;
                        }
                }
                else std::cerr << "Unable to write to " << filename << std::endl;
                out.close();
        }
}

/**
 * @brief This writes out the contents of the GoftCube.
 * 
 * This member function writes out the contents of the GoftCube to the hard disk, in the filename given by the argument. 
 * The file can then be read into IDL with the commands provide under the idl/ directory, here in particular readgoftcube.pro.
 * Normally, one would use it to post-process the forward modelling results (when using with a RenderCube), but it is also very useful for debugging purposes when using it with a GoftCube before the rendering.
 * @param filename This parameter specifies which filename the data needs to be written to.
 */
void writegoftcube_binary(FoMo::GoftCube goftcube, const std::string filename)
{
	// write out goftcube to file "filename"
	int commrank;
	std::string space=" ";
	size_t chiantisize=goftcube.readchiantifile().size();
	size_t abundsize=goftcube.readabundfile().size();
	
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	if (commrank==0)
	{
		std::ofstream out(filename,std::ios::binary|std::ios::ate);
		if (out.is_open())
		{
			std::string versionstring("FoMo-");
			versionstring.append(FOMO_VER).append("#"); // this allows us to read the file into python, up to the marker string #
			// the longest possible string is 27 digits (e.g. FoMo-v3.3-79-gbeba0fc-dirty), 
			// but could grow longer if 3.3 changes to e.g. 3.10 or there are more than 99 commits since the last tag
			// in any case, it will be shorter than 40 digits
			out.write(versionstring.c_str(),versionstring.size());
			int nvars = goftcube.readnvars();
			int dim = goftcube.readdim();
			int ng = goftcube.readngrid();
			out.write(reinterpret_cast<const char*>(&dim),sizeof(dim));
			out.write(reinterpret_cast<const char*>(&ng),sizeof(ng));
			out.write(reinterpret_cast<const char*>(&nvars),sizeof(nvars));
			// units
			std::vector<std::string> unitvec=goftcube.readunit();
			for (unsigned int i=0; i<dim+nvars; i++)
			{
				size_t unitsize=unitvec.at(i).size();
				out.write(reinterpret_cast<const char*>(&unitsize),sizeof(unitsize));
				out.write(unitvec.at(i).c_str(),unitsize);
			}
			out.write(reinterpret_cast<const char*>(&chiantisize),sizeof(chiantisize));
			out.write(goftcube.readchiantifile().c_str(),chiantisize);
			out.write(reinterpret_cast<const char*>(&abundsize),sizeof(abundsize));
			out.write(goftcube.readabundfile().c_str(),abundsize);
			
			FoMo::tgrid grid=goftcube.readgrid();
			FoMo::tvars vars;
			for (int i=0; i<nvars; i++)
			{
				vars.push_back(goftcube.readvar(i));
			}
			for (int i=0; i<dim; i++)
			{ 
				out.write(reinterpret_cast<const char*>(grid[i].data()),ng*sizeof(grid[i][0]));
			}
			for (int i=0; i<nvars; i++)
			{
				out.write(reinterpret_cast<const char*>(vars[i].data()),ng*sizeof(vars[i][0]));
			}
		}
		else std::cerr << "Unable to write to " << filename << std::endl;
		out.close();
	}
}

const std::vector<std::string> extensions={".dat", ".txt"};

void zipfile(std::bitset<FoMo::noptions> woptions, std::string root)
{
	for (unsigned int i=0; i<woptions.size()-2; i++)
	{
		if (woptions[i])
		{
			// zipfile(root+extensions[i]);
			std::string filename=root+extensions[i];
			std::ifstream infile(filename.c_str());
			std::ofstream outfile(filename+".gz",std::ofstream::binary);
			boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
			in.push(boost::iostreams::gzip_compressor());
			in.push(infile);
			boost::iostreams::copy(in,outfile);
		}	
	}
}

/**
 * @brief This member returns the current writeoptions
 * 
 * The current writeoptions are returned as a std::bitset<FoMo::noptions>.
 * @return The current writeoptions.
 */
std::bitset<FoMo::noptions> FoMo::GoftCube::getwriteoptions()
{
	return writeoptions;
}

/**
 * @brief This members sets the output options
 * 
 * The valid options are given as a std::bitset<4>. The options are bits, which can be combined if needed, e.g. to zip a binary output file. 
 * options[0] says if a binary file should be written, options[1] states if a text file should be written, options[2] declares if the 
 * previous output files should be zipped. Options[4] says if the non-zipped files should be deleted. 
 * @param options A bitset<4> of options to be set. 
 */
void FoMo::GoftCube::setwriteoptions(std::bitset<FoMo::noptions> woptions)
{
	writeoptions=woptions;
}

/**
 * @brief This member forces the writeout of a binary output
 * 
 * If writebinary is true, the output will be written to a binary file (with extension .dat)
 * @param writebinary If true, a binary output will be written. If false, it will not be written. 
 */
void FoMo::GoftCube::setwriteoutbinary(const bool writebinary)
{
	writeoptions[0]=writebinary;
}

/**
 * @brief This member forces the writeout of a text output
 * 
 * If writetext is true (default), the output will be written to a text file (with extension .txt)
 * @param writetext If true, a text output will be written. If false, it will not be written.
 */
void FoMo::GoftCube::setwriteouttext(const bool writetext)
{
	writeoptions[1]=writetext;
}

/**
 * @brief This member forces the zipping of the output files
 * 
 * If writezip is true, the output files (in binary .dat or text .txt) will be zipped. 
 * @param writezip If true, the output files will be zipped. If false, they remain in their original form.
 */
void FoMo::GoftCube::setwriteoutzip(const bool writezip)
{
	writeoptions[2]=writezip;
}

/**
 * @brief This member forces the deletion of files that were previously zipped
 * 
 * If deletefiles is true, the previous output files (binary .dat or text .txt) will be removed after a .dat.zip or .txt.zip will have been created. 
 * @param deletefiles If true, previous files will be removed. If false, they will be kept. 
 */
void FoMo::GoftCube::setwriteoutdeletefiles(const bool deletefiles)
{
	writeoptions[3]=deletefiles;
}
 
/**
 * @brief This writes out the contents of the GoftCube.
 * 
 * This member function writes out the contents of the GoftCube to the hard disk, in the filename given by the argument (with the extension appropriately changed to .txt, .dat or .zip). 
 * The file can then be read into IDL with the commands provide under the idl/ directory, here in particular readgoftcube.pro.
 * Normally, one would use it to post-process the forward modelling results (when using with a RenderCube), but it is also very useful for debugging purposes when using it with a GoftCube before the rendering.
 * @param filename This parameter specifies which filename the data needs to be written to.
 */
void FoMo::GoftCube::writegoftcube(const std::string filename)
{
	// find root of filename
	// append correct extension
	size_t pos = filename.rfind(".");
	std::string root;
    if ((pos == std::string::npos) ||   //No extension.
		(pos == 0)) //. is at the front. Not an extension.
		{
			root=filename;
		}
	else
	{
		root=filename.substr(0, pos);
	}
	// if binary requested: writeoptions[0] is true
	if (writeoptions[0])
		writegoftcube_binary(*this,root+extensions[0]);
	// if text requested: writeoptions[1] is true
	if (writeoptions[1])
		writegoftcube_txt(*this,root+extensions[1]);
	// if zip requested: writeoptions[2] is true
	if (writeoptions[2] & (writeoptions[0] | writeoptions[1]))
	// zip filename
	{
		zipfile(writeoptions, root);
		// if delete files requested, remove uncompressed files
		if (writeoptions[3])
		{
			// delete files 
			for (unsigned int i=0; i<writeoptions.size()-2; i++)
			{
				std::string tmpfile=root+extensions[i];
				if (writeoptions[i]) remove(tmpfile.c_str());
			}
		}
	}
	else if (writeoptions[2])
		std::cerr << "Zipping of output only possible if output set with goftcube.setwriteouttext or goftcube.setwriteoutbinary." << std::endl << std::flush;
}

/**
 * @brief This is a routine for reading a previously saved GoftCube.
 * 
 * This member reads the information from file emissionsave, which was previously written by GoftCube::writegoftcube, and stores it in the GoftCube.
 * @param emissionsave The filename from which the GoftCube should be read.
 */
void FoMo::GoftCube::readgoftcube(const std::string emissionsave)
{
	// Read GoftCube datacube from file "emissionsave"
	// Otherwise read in data cube (simulation snapshots) from directory specified in main function

	int commrank;
#ifdef HAVEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
#else
	commrank = 0;
#endif
	if (commrank==0)
	{
		std::ifstream in(emissionsave,std::ios::binary);
		if (in.is_open())
		{
			in >> dim;
			in >> ng;
			in >> nvars;
			in >> chiantifile;
			in >> abundfile;
			
			grid.resize(dim);
			vars.resize(nvars);
			lambda0=readgoftfromchianti(chiantifile);
			
			for (unsigned int i=0; i<dim; i++)
			{
				grid[i].resize(ng);
			}
			for (unsigned int i=0; i<nvars; i++)
			{
				vars[i].resize(ng);
			}
			for (unsigned int j=0; j<ng; j++)
			{
				for (unsigned int i=0; i<dim; i++) in >> grid[i][j];
				for (unsigned int i=0; i<nvars; i++) in >> vars[i][j];
			}
		}
		else std::cerr << "Unable to read " << emissionsave << std::endl;
		in.close();
	}
}