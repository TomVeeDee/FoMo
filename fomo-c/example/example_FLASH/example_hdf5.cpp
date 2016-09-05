// N. Magyar - 20th August 2015
// Procedure to read 3D AMR FLASH data file (plotfile) into FoMo Object. 
// HDF5 installed is required (with C++ API).

#include <iostream>
#include <string>
#include <cmath>
#include <omp.h>
#include "H5Cpp.h"
#include "FoMo.h"

using namespace std; 

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#ifndef H5_NO_STD
    using std::cout;
    using std::cin;
    using std::endl;
#endif
#endif

//---------------------------------------------
//The name of the file being read - edit this
//---------------------------------------------
const H5std_string	FILE_NAME("example.hdf5");

const H5std_string	DENS_NAME("dens");
const H5std_string	TEMP_NAME("temp");
const H5std_string	VELX_NAME("velx");
const H5std_string	VELY_NAME("vely");
const H5std_string	VELZ_NAME("velz");
const H5std_string	BOUND_NAME("bounding box");
const H5std_string      REFINE_NAME("refine level");
const int       RANK = 3;


int main (void)
{
       
    int i, j, k, l, nr;
    int nxa,nxb,nxc,blocks;
    float dx,dy,dz;
    float lnorm, dnorm, tnorm, vnorm;

    try
    {
      double start = omp_get_wtime();
      // Open an existing file and datasets.
      H5File file(FILE_NAME, H5F_ACC_RDWR);
      DataSet densSet = file.openDataSet(DENS_NAME);
      DataSet tempSet = file.openDataSet(TEMP_NAME);
      DataSet velxSet = file.openDataSet(VELX_NAME);
      DataSet velySet = file.openDataSet(VELY_NAME);
      DataSet velzSet = file.openDataSet(VELZ_NAME);
      DataSet boundSet = file.openDataSet(BOUND_NAME);
      DataSet refineSet = file.openDataSet(REFINE_NAME);
       
      // Turn off the auto-printing when failure occurs so that we can
      // handle the errors appropriately
      Exception::dontPrint();
      
      //Find out extent of the dataset  
      H5T_class_t type_class = densSet.getTypeClass();
        
      //Find out endianness  
      FloatType inttype = densSet.getFloatType();
      H5std_string order_string;
      H5T_order_t order = inttype.getOrder(order_string); 
        
      cout << order_string << endl;
         
      //Print byte size (4 - single, 8 - double)
      size_t size = inttype.getSize();
      cout << "Data size is " << size << endl;
 
      DataSpace densSpace = densSet.getSpace();
      DataSpace tempSpace = tempSet.getSpace();
      DataSpace velxSpace = velxSet.getSpace();
      DataSpace velySpace = velySet.getSpace();
      DataSpace velzSpace = velzSet.getSpace();
      DataSpace boundSpace = boundSet.getSpace();
      DataSpace refineSpace = refineSet.getSpace();
 
      // rank - dimensions - 4 for 3D data
      int rank = densSpace.getSimpleExtentNdims();

      //Get dimension sizes and print
      hsize_t nd[rank];
      int ndims = densSpace.getSimpleExtentDims(nd, NULL);
      cout << "Rank: " << rank << ", dimensions: ";
      for (i = 0; i<rank; i++)
         {
              if (i < rank-1){
                 cout << (unsigned long)(nd[i]) << " x ";}
              else{
                 cout << (unsigned long)(nd[i]) << endl;}
         }
          
       hsize_t var_dims[3];
       hsize_t bound_dims[2];
       hsize_t block_dims[1];

       var_dims[0] = nd[1]; var_dims[1] = nd[2]; var_dims[2] = nd[3];
       bound_dims[0] = 3; bound_dims[1] = 2; 
       block_dims[0] = nd[0];

       nxa = nd[3]; nxb = nd[2]; nxc = nd[1]; blocks = nd[0];

       DataSpace space(RANK,var_dims);
       DataSpace bound_space(2,bound_dims);
       DataSpace refine_space(1,block_dims);

       //Arrays to read in variables
       float* dens = new float[nxa*nxb*nxc]();
       float* temp = new float[nxa*nxb*nxc]();
       float* velx = new float[nxa*nxb*nxc]();
       float* vely = new float[nxa*nxb*nxc]();
       float* velz = new float[nxa*nxb*nxc]();
       
       //Arrays for coordinates
       float* bound = new float[6]();
       float* coord = new float[nxa*nxb*nxc*3]();
       float* refine = new float[blocks]();
  
       //myfile.open("test.dat");
 
       //Create FoMoObject
       FoMo::FoMoObject Object;

       //
       cout << "No. of cells in one block: " << nxa*nxb*nxc << endl;

       //Read all the blocks, one by one, for coordinates and variables
       cout << "Reading in ... " << endl;

       //Read level of refinement of each block
       refineSet.read(refine,PredType::NATIVE_FLOAT,refine_space,refineSpace);
       
       for (i = 0; i < blocks; i++)
       {
          //Deal only with the highest refinement level blocks
          if (refine[i] >= refine[i+1]) 
          {
          //Read 1 block from the variables and put it in 1D array
          hsize_t offset[4] = {i,0,0,0};
          hsize_t count[4] = {1,nxc,nxb,nxa};
 
          densSpace.selectHyperslab(H5S_SELECT_SET, count, offset);
          densSet.read(dens, PredType::NATIVE_FLOAT, space, densSpace);

          tempSpace.selectHyperslab(H5S_SELECT_SET, count, offset);          
          tempSet.read(temp, PredType::NATIVE_FLOAT, space, tempSpace);

          velxSpace.selectHyperslab(H5S_SELECT_SET, count, offset);
          velxSet.read(velx, PredType::NATIVE_FLOAT, space, velxSpace);
          
          velySpace.selectHyperslab(H5S_SELECT_SET, count, offset);
          velySet.read(vely, PredType::NATIVE_FLOAT, space, velySpace);
          
          velzSpace.selectHyperslab(H5S_SELECT_SET, count, offset);
          velzSet.read(velz, PredType::NATIVE_FLOAT, space, velzSpace);

          //Read the coordinate and size of corresponding blocks
          hsize_t bound_offset[3] = {i,0,0}; 
          hsize_t bound_count[3] = {1,3,2};
          boundSpace.selectHyperslab(H5S_SELECT_SET, bound_count, bound_offset);
          boundSet.read(bound, PredType::NATIVE_FLOAT, bound_space, boundSpace); 
 
          //Construct cell coordinates
          dx = (bound[1] - bound[0])/nxa;
          dy = (bound[3] - bound[2])/nxb;
          dz = (bound[5] - bound[4])/nxc;
 
          nr = 0;
          for (j = 0; j < nxc; j++)
            for (k = 0; k < nxb; k++)
               for (l = 0; l < nxa; l++){ 
                 coord[nr] = bound[0] + dx/2 + dx*l; 
                 coord[nr+1] = bound[2] + dy/2 + dy*k;
                 coord[nr+2] = bound[4] + dz/2 + dz*j;
                 nr = nr + 3;
               }
          
          //Put the data into FoMoObject
          for (j = 0; j < nxa*nxb*nxc; j++){
 
            //Arrays for FoMoObject
            vector<double> coordinates;
            vector<double> variables;
            
            //here you can define normalization units, so that lengths are in km, densities in cm^-3,
            //temperature in K, and speeds in m/s.
             
            lnorm = 1.e5; //from 100 Mm to km
            dnorm = 1.e-12 * 1.204 * 1.e21; // from 10^-12 kg m^-3 to cm^-3 for hydrogen plasma
            tnorm = 1.e0; // K, unchanged
            vnorm = 1.e6; // from Mm/s to m/s

            coordinates.push_back(coord[3*j]*lnorm); 
            coordinates.push_back(coord[3*j+1]*lnorm); 
            coordinates.push_back(coord[3*j+2]*lnorm);

            variables.push_back(dens[j]*dnorm); 
            variables.push_back(temp[j]*tnorm); 
            variables.push_back(velx[j]*vnorm); 
            variables.push_back(vely[j]*vnorm); 
            variables.push_back(velz[j]*vnorm);
        
            /*myfile << coordinates[0] << " " <<
                      coordinates[1] << " " <<
                      coordinates[2] << " " <<
                      variables[0] << " " <<
                      variables[1] << " " <<
                      variables[2] << " " <<
                      variables[3] << " " <<
                      variables[4] << endl;*/
            
            Object.push_back_datapoint(coordinates,variables);
          }
          }  
        }

    //Deallocating
    
    delete [] dens;
    delete [] temp;
    delete [] velx;
    delete [] vely;
    delete [] velz;
    delete [] bound;
    delete [] coord;
    delete [] refine;
 
    //myfile.close();

    // data is in structure, now start the rendering
	
    /// [Set rendering options]
    Object.setchiantifile("/users/cpa/tomvd/data/idl/FoMo/chiantitables/goft_table_aia171_abco.dat"); // the default value is "../chiantitables/goft_table_fe_12_0194small_abco.dat"
    Object.setabundfile("/empty"); //use "/empty" or do not set it at all for the default sun_coronal_2012_schmelz.abund file
    Object.setrendermethod("CGAL"); // CGAL is the default rendermethod
    Object.setobservationtype(FoMo::Imaging);
    // adjust the resolution with these parameters
    int x_pixel=128;
    int y_pixel=200;
    int z_pixel=200;
    int lambda_pixel=1;
    double lambda_width=.13;
    Object.setresolution(x_pixel,y_pixel,z_pixel,lambda_pixel,lambda_width);
    // determine where the output will be written
    Object.setoutfile("fomo-hdf5-example-out.");
    /// [Set rendering options]
	
    /// [Render]
    //Object.render(2*atan(1.),2*atan(1.));
      Object.render(0.,0.);
    // alternatively, you could add different angles in radians as argument, e.g. Object.render(1.57/2.,1.57/2.) to obtain some nice Doppler shifts.
    /// [Render]
	
    /// [Details] 
    //    //Write out goftcube
    FoMo::GoftCube cube=Object.readgoftcube();
    cube.writegoftcube("goftcube.txt");
    
    // read the rendering
	
    FoMo::RenderCube rendercube=Object.readrendering();
    
    // read the resolution
    int nx,ny,nz,nlambda;
    double lambdawidth;
    rendercube.readresolution(nx,ny,nz,nlambda,lambdawidth);
    cout << "The rendering has resolution of x and y " << nx << " " << ny << ".\n";
    cout << "and " << nlambda << " pixels in the wavelength." << endl;
    cout << "It was done using the rendermethod " << rendercube.readrendermethod() << endl;
    cout << "with a resolution of " << nz << " along the line-of-sight." << endl;
    // This writes out the rendering results to the file fomo-output.txt.
    
    rendercube.writegoftcube("fomo-hdf5-output.txt");
    /// [Details]
    
    double end_time = omp_get_wtime() - start;
    
    cout << end_time << " s" << endl;
    }  // end of try block

    // catch failure caused by the H5File operations
    catch(FileIException error)
    {
	error.printError();
	return -1;
    }

    // catch failure caused by the DataSet operations
    catch(DataSetIException error)
    {
	error.printError();
	return -1;
    }
  
    return 0;  // successfully terminated
}
