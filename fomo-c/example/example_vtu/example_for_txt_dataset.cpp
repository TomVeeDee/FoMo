#include "FoMo.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>

/**
 * @file 
 * This file contains a basic example for using FoMo. The example program takes as argument a file name.
 * The file name in the argument should contain a list of 3D coordinates with associated variables (rho, T, vx, vy, vz).
 */
namespace fs = std::filesystem;

using namespace std;

void processTxtFile(const string& inputFilename, const string& outputFilename) {
    // Initialize the FoMo object
    FoMo::FoMoObject Object;

    // Open the file as fstream
    ifstream filetoread(inputFilename);

    if (filetoread.is_open()) {
        double tmpvar;
        while (!filetoread.eof()) {
            vector<double> coordinates;
            // The points should be three-dimensional
            for (unsigned int i = 0; i < 3; i++) {
                filetoread >> tmpvar;
                coordinates.push_back(tmpvar);
            }

            vector<double> variables;
            // There should be 5 variables (n, T, vx, vy, vz)
            for (unsigned int i = 0; i < 5; i++) {
                filetoread >> tmpvar;
                variables.push_back(tmpvar);
            }

            if (!filetoread.eof()) 
                Object.push_back_datapoint(coordinates, variables);
        }
        filetoread.close();
    } else {
        cerr << "Failed to open file: " << inputFilename << endl;
        return;
    }

    // Set rendering options
    Object.setchiantifile("../chiantitables/goft_table_fe_12_0194_abco.dat");
    Object.setabundfile("/empty");
    Object.setrendermethod("Thomson");
    Object.setobservationtype(FoMo::Imaging);

    // Adjust the resolution with these parameters
    int x_pixel = 2000;
    int y_pixel = 2000;
    int z_pixel = 2000;
    int lambda_pixel = 1000;
    double lambda_width = 200000;  // in m/s
    Object.setresolution(x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width);

    // Determine where the output will be written
    Object.setoutfile(outputFilename);

    // Render the data
    Object.render(0.0,0.0);
//1.57/3.0,1.57/2.0);

    // Read the rendering result from the object
    FoMo::RenderCube rendercube = Object.readrendering();

    // Set the output to be written in binary and zipped format
    rendercube.setwriteoutbinary();
    rendercube.setwriteoutzip();
    rendercube.setwriteoutdeletefiles();

    rendercube.writegoftcube(outputFilename + ".txt");
    cout << "Processed and output saved to: " << outputFilename << ".txt" << endl;
}

void processDirectory(const std::string& inputDirectory, const std::string& outputDirectory) {
    for (const auto& entry : fs::directory_iterator(inputDirectory)) {
        if (entry.path().extension() == ".txt") {
            std::string inputFilename = entry.path().string();
            std::string outputFilename = outputDirectory + "/" + entry.path().stem().string() + "_fomo_output";
            processTxtFile(inputFilename, outputFilename);
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_directory> <output_directory>" << endl;
        return -1;
    }

    std::string inputDirectory = argv[1];
    std::string outputDirectory = argv[2];

    // Create the output directory if it doesn't exist
    fs::create_directories(outputDirectory);

    processDirectory(inputDirectory, outputDirectory);

    return 0;
}
