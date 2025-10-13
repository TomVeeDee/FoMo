#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include "FoMo.h"

namespace fs = std::filesystem;
using namespace std;

// Function to process a single .vtu file and convert it to a FoMoObject
void processVtuFile(const std::string& inputFilename, const std::string& outputFilename) {
    try {
        // Create a VTK reader for .vtu files
        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        reader->SetFileName(inputFilename.c_str());
        reader->Update();

        // Extract grid data
        vtkSmartPointer<vtkUnstructuredGrid> grid = reader->GetOutput();
        if (!grid) {
            cerr << "Failed to read file: " << inputFilename << endl;
            return;
        }

        // Get the grid points
        vtkSmartPointer<vtkPoints> points = grid->GetPoints();
        if (!points) {
            cerr << "No points found in file: " << inputFilename << endl;
            return;
        }

        // Extract data fields
        vtkFloatArray* pArray = dynamic_cast<vtkFloatArray*>(grid->GetPointData()->GetArray("p"));
        vtkFloatArray* rhoArray = dynamic_cast<vtkFloatArray*>(grid->GetPointData()->GetArray("rho"));
        vtkFloatArray* vArray = dynamic_cast<vtkFloatArray*>(grid->GetPointData()->GetArray("v"));

        // Check for missing arrays
        if (!pArray || !rhoArray || !vArray) {
            cerr << "Missing data arrays in file: " << inputFilename << endl;
            return;
        }

        // Ensure the 'v' array has three components
        if (vArray->GetNumberOfComponents() != 3) {
            cerr << "'v' array does not have 3 components in file: " << inputFilename << endl;
            return;
        }

        // Create a FoMo object
        FoMo::FoMoObject Object;

        // Read and store data
        vtkIdType numPoints = points->GetNumberOfPoints();
        for (vtkIdType i = 0; i < numPoints; i++) {
            vector<double> coordinates;
            double p[3];
            points->GetPoint(i, p);
            
            // Scale coordinates
            coordinates.push_back(p[0] * 696.0);
            coordinates.push_back(p[1] * 696.0);
            coordinates.push_back(p[2] * 696.0);

            vector<double> variables;
            float pressure = pArray->GetValue(i);
            float rho = rhoArray->GetValue(i);
            float rho_cm = rho * 1e8;
            float T = (pressure / rho) * 1.7756e7;

            // Get velocity components
            double vComponents[3];
            vArray->GetTuple(i, vComponents);
            float vX = vComponents[0] * 4.8e5;
            float vY = vComponents[1] * 4.8e5;
            float vZ = vComponents[2] * 4.8e5;

            // Store variables (n = rho_cm, T, vx, vy, vz)
            variables.push_back(rho_cm);
            variables.push_back(T);
            variables.push_back(vX);
            variables.push_back(vY);
            variables.push_back(vZ);

            Object.push_back_datapoint(coordinates, variables);
        }

        // Set rendering options
        Object.setchiantifile("../chiantitables/goft_table_fe_12_0194_abco.dat");
        Object.setabundfile("/empty");
        Object.setrendermethod("Thomson");
        Object.setobservationtype(FoMo::Imaging);

        // Set the resolution parameters
        int x_pixel = 2000;
        int y_pixel = 2000;
        int z_pixel = 2000;
        int lambda_pixel = 1000;
        double lambda_width = 200000;  // in m/s
        Object.setresolution(x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width);

        // Set the output file
        Object.setoutfile(outputFilename);

        // Render the data
        Object.render(0.0, 0.0);

        // Read the rendering result
        FoMo::RenderCube rendercube = Object.readrendering();

        // Set output options
        rendercube.setwriteoutbinary();
        rendercube.setwriteoutzip();
        rendercube.setwriteoutdeletefiles();

        // Write output file
        rendercube.writegoftcube(outputFilename + ".txt");
        cout << "Processed and output saved to: " << outputFilename << ".txt" << endl;

    } catch (const std::exception& e) {
        cerr << "Error processing file " << inputFilename << ": " << e.what() << endl;
    }
}

// Function to process all .vtu files in a directory
void processDirectory(const std::string& inputDirectory, const std::string& outputDirectory) {
    for (const auto& entry : fs::directory_iterator(inputDirectory)) {
        if (entry.path().extension() == ".vtu") {
            std::string inputFilename = entry.path().string();
            std::string outputFilename = outputDirectory + "/" + entry.path().stem().string() + "_fomo_output";
            processVtuFile(inputFilename, outputFilename);
        }
    }
}

// Main function to process a directory of .vtu files
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