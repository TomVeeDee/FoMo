#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include "FoMo.h"

namespace fs = std::filesystem;
using namespace std;

void processVtuFile(const std::string& inputFilename, const std::string& outputFilename) {
    try {
        // Read the .vtu file
        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        reader->SetFileName(inputFilename.c_str());
        reader->Update();

        vtkSmartPointer<vtkUnstructuredGrid> grid = reader->GetOutput();

        // Retrieve CellData arrays
        vtkFloatArray* rhoArray = vtkFloatArray::SafeDownCast(grid->GetCellData()->GetArray("rho"));
        vtkFloatArray* v1Array = vtkFloatArray::SafeDownCast(grid->GetCellData()->GetArray("v1"));
        vtkFloatArray* v2Array = vtkFloatArray::SafeDownCast(grid->GetCellData()->GetArray("v2"));
        vtkFloatArray* v3Array = vtkFloatArray::SafeDownCast(grid->GetCellData()->GetArray("v3"));
        vtkFloatArray* pArray  = vtkFloatArray::SafeDownCast(grid->GetCellData()->GetArray("p"));

        if (!rhoArray || !v1Array || !v2Array || !v3Array || !pArray) {
            std::cerr << "Missing one or more required arrays (rho, v1, v2, v3, p) in file: " << inputFilename << std::endl;
            return;
        }

        vtkSmartPointer<vtkPoints> points = grid->GetPoints();
        if (!points) {
            std::cerr << "No points found in file: " << inputFilename << std::endl;
            return;
        }

        // Create a FoMo object
        FoMo::FoMoObject Object;
        vtkIdType numCells = grid->GetNumberOfCells();

        for (vtkIdType i = 0; i < numCells; ++i) {
            vtkCell* cell = grid->GetCell(i);

            // Compute centroid of the cell in Cartesian space
            double centroid[3] = {0.0, 0.0, 0.0};
            int numPointsInCell = cell->GetNumberOfPoints();

            for (int j = 0; j < numPointsInCell; ++j) {
                double point[3];
                points->GetPoint(cell->GetPointId(j), point);
                centroid[0] += point[0];
                centroid[1] += point[1];
                centroid[2] += point[2];
            }

            centroid[0] /= numPointsInCell;
            centroid[1] /= numPointsInCell;
            centroid[2] /= numPointsInCell;

            // Use coordinates directly (no spherical conversion)
            std::vector<double> coordinates = {
                centroid[0] * 696.0, // scale if needed
                centroid[1] * 696.0,
                centroid[2] * 696.0
            };

            // Get physical quantities
            float rho = rhoArray->GetValue(i);
            float p   = pArray->GetValue(i);
            float T   = p / rho;

            float vX = v1Array->GetValue(i);
            float vY = v2Array->GetValue(i);
            float vZ = v3Array->GetValue(i);

            // Store in FoMo object
            std::vector<double> variables = {rho, T, vX, vY, vZ};
            Object.push_back_datapoint(coordinates, variables);
        }

        // Rendering setup
        Object.setrendermethod("Thomson");
        Object.setobservationtype(FoMo::Imaging);

        int x_pixel = 500, y_pixel = 500, z_pixel = 500, lambda_pixel = 1000;
        double lambda_width = 100000;  // in m/s
        Object.setresolution(x_pixel, y_pixel, z_pixel, lambda_pixel, lambda_width);
        Object.setoutfile(outputFilename);
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

void processDirectory(const std::string& inputDirectory, const std::string& outputDirectory) {
    cout << "Scanning directory: " << inputDirectory << endl;
    for (const auto& entry : fs::directory_iterator(inputDirectory)) {
        if (entry.path().extension() == ".vtu") {
            string inputFilename = entry.path().string();
            string outputFilename = outputDirectory + "/" + entry.path().stem().string() + "_fomo_output";
            processVtuFile(inputFilename, outputFilename);
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_directory> <output_directory>" << endl;
        return -1;
    }

    string inputDirectory = argv[1];
    string outputDirectory = argv[2];

    fs::create_directories(outputDirectory);
    processDirectory(inputDirectory, outputDirectory);
    
    return 0;
}

