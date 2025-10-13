#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "FoMo.h"

namespace fs = std::filesystem;
using namespace std;

void processCfmeshFile(const std::string& inputFilename, const std::string& outputFilename) {
    try {
        std::ifstream file(inputFilename);
        if (!file) {
            std::cerr << "Failed to open file: " << inputFilename << std::endl;
            return;
        }

        std::string line;
        std::vector<std::vector<double>> nodeCoords;  // Stores node coordinates (indexed)
        std::vector<std::vector<int>> elements;       // Stores element-to-node mappings
        std::vector<int> stateIDs;                    // Stores state IDs

        std::vector<float> p;   // Pressure
        std::vector<float> rho; // Density
        std::vector<std::vector<float>> v; // Velocity components

        bool readingNodes = false;
        bool readingElems = false;
        bool readingStates = false;

        while (std::getline(file, line)) {
            //line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
            if (line.empty()) continue;

            if (line.find("!LIST_NODE") != std::string::npos) {
                readingNodes = true;
                readingElems = false;
                readingStates = false;
                continue;
            } else if (line.find("!LIST_ELEM") != std::string::npos) {
                readingNodes = false;
                readingElems = true;
                readingStates = false;
                continue;
            } else if (line.find("!LIST_STATE") != std::string::npos) {
                readingNodes = false;
                readingElems = false;
                readingStates = true;
                continue;
            }
                // If line starts with "!" but is not part of the list markers, stop reading sections
    		else if (line[0] == '!' && 
             		line.find("!LIST_NODE") == std::string::npos &&
             		line.find("!LIST_ELEM") == std::string::npos &&
             		line.find("!LIST_STATE") == std::string::npos) {
             	readingNodes = false;
        		readingElems = false;
        		readingStates = false;
        		continue;
    		}

            std::istringstream iss(line);

            // Reading node coordinates
            if (readingNodes) {
                double x, y, z;
                iss >> x >> y >> z;
                nodeCoords.push_back({x, y, z});
            }

            else if (readingElems) {
    			std::vector<int> elem(6);
    			int state;
    			iss >> elem[0] >> elem[1] >> elem[2] >> elem[3] >> elem[4] >> elem[5] >> state;
    			if (iss) { // Only add if reading was successful
        			elements.push_back(elem);
        			stateIDs.push_back(state);
    			} else {
        			std::cerr << "Error: Invalid element line, skipping: " << line << std::endl;
    			}
			}

            // Reading state data (rho, p, velocity)
            else if (readingStates) {
    // Read the first value as density and ignore the rest of the values for this line
    			float rho_val;
    			if (iss >> rho_val) {
        // Store the density value
        			rho.push_back(rho_val);

        // Skip the rest of the values on the line (assuming they are irrelevant)
        			float dummy;
        			for (int i = 0; i < 8; ++i) {  // Skip the next 8 values (pressure, velocity, etc.)
           				iss >> dummy;
        			}
    			}
			}
            
        }
		std::cout << "elements.size(): " << elements.size() << std::endl;
		std::cout << "stateIDs.size(): " << stateIDs.size() << std::endl;
		std::cout << "rho.size(): " << rho.size() << std::endl;


        // Ensure data consistency
        if (elements.size() != stateIDs.size() || elements.size() != rho.size() ) {
            std::cerr << "Error: Mismatch in element data sizes!" << std::endl;
            return;
        }

        // FoMo object creation
        FoMo::FoMoObject Object;

        // Process elements and compute centroids and variables
        for (size_t i = 0; i < elements.size(); ++i) {
            const auto& elem = elements[i];
            int stateIndex = stateIDs[i];
            
            // Debug print to check element and stateIndex
    		std::cout << "Processing element " << i << " with stateIndex: " << stateIndex << std::endl;
    		
            // Compute average node coordinate for the element
            std::vector<double> avgCoord(3, 0.0);
            for (int nid : elem) {
                for (int d = 0; d < 3; ++d)
                    avgCoord[d] += nodeCoords[nid][d];
            }
            for (int d = 0; d < 3; ++d)
                avgCoord[d] = (avgCoord[d] / 6.0) * 696.0;  // scale to Mm

            if (stateIndex < rho.size()) {
            //float pressure = p[stateIndex];
            	float density = rho[stateIndex];
            	float rho_cm = density * 1e8;
            	float T = 0.0f;  // Zero temperature
        		float vx = 0.0f; // Zero velocity
        		float vy = 0.0f;
        		float vz = 0.0f;
            
            //float T = (pressure / density) * 1.7756e7;

            //float vx = v[stateIndex][0] * 4.8e5;
            //float vy = v[stateIndex][1] * 4.8e5;
            //float vz = v[stateIndex][2] * 4.8e5;

            	std::vector<double> variables = {rho_cm, T, vx, vy, vz};
            	Object.push_back_datapoint(avgCoord, variables);
        	} else {
        		std::cerr << "Error: stateIndex " << stateIndex << " is out of bounds!" << std::endl;
			}
		}
		
        // Rendering and output (same as before)
        //Object.setchiantifile("../chiantitables/goft_table_fe_12_0194_abco.dat");
        Object.setabundfile("/empty");
        Object.setrendermethod("Thomson");
        Object.setobservationtype(FoMo::Imaging);

        Object.setresolution(1500, 1500, 1500, 1000, 200000);
        Object.setoutfile(outputFilename);
        Object.render(0.0, 1.57);

        FoMo::RenderCube rendercube = Object.readrendering();
        rendercube.setwriteoutbinary();
        rendercube.setwriteoutzip();
        rendercube.setwriteoutdeletefiles();
        rendercube.writegoftcube(outputFilename + ".txt");

        std::cout << "Processed and saved: " << outputFilename << ".txt" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

// Function to process all .vtu files in a directory
void processDirectory(const std::string& inputDirectory, const std::string& outputDirectory) {
    for (const auto& entry : fs::directory_iterator(inputDirectory)) {
        if (entry.path().extension() == ".CFmesh") {
            std::string inputFilename = entry.path().string();
            std::string outputFilename = outputDirectory + "/" + entry.path().stem().string() + "_fomo_output";
            processCfmeshFile(inputFilename, outputFilename);
        }
    }
}

// Main function to process a directory of files
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