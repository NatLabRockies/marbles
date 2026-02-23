#include <stdio.h>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <string>
#include <filesystem>
#include <chrono>
#include <cstdlib>

typedef float myReal;

int main(int argc, char* argv[])
{
    // Get geometry seed from command line argument or use default
    int geometrySeed = 12345; // default seed
    if (argc > 1) {
        geometrySeed = std::atoi(argv[1]);
    }

    std::cout << "Using geometry seed: " << geometrySeed << std::endl;

    const int nX(60), nY(40), nZ(30); // grid

    std::vector<myReal> diameterVector{7, 6, 5}; // diameter of tubes

    std::vector<myReal> volumeFractions{
        0.2, 0.3, 0.4}; // relativevolume fractions of corresponding tubes
    // Need not sum up to 1, will be automatically rescaled to 1

    myReal clusteringFactor =
        0.001; // reduce for more clustering, increase for less clustering

    myReal tubeVolumeFraction =
        0.15; // total volume fraction of all tubes, i.e. the empty space

    bool forMARBLES_LBM_code =
        true; // makes tubes 0 and rest of the domain of solid space is 1.
              // Otherwise, tubes are >0 and solid space is 0

    // Scale the volume fractions to a total of 1.0
    myReal sumVF = 0.0;
    for (const auto& vf : volumeFractions) sumVF += vf;
    for (auto& vf : volumeFractions)
        vf = vf / sumVF; // rescale so that the sum is 1.

    for (size_t i = 0; i < volumeFractions.size(); ++i)
        volumeFractions[i] *= tubeVolumeFraction;

    // no user input below this line
    std::vector<myReal> radiusVector(diameterVector.size());
    for (size_t i = 0; i < diameterVector.size(); ++i)
        radiusVector[i] = diameterVector[i] / 2.0; // convert diameter to radius

    const long maxNumberOfSpheres = static_cast<long>(nX) * nY * nZ;
    int maxRadius =
        std::ceil(*std::max_element(radiusVector.begin(), radiusVector.end()));
    const int pX(maxRadius), pY(maxRadius), pZ(maxRadius); // dummy points
    const int nsX(nX + 2 * pX), nsY(nY + 2 * pY), nsZ(nZ + 2 * pZ);

    // Use a more memory-efficient 3D structure instead of 4D vector
    std::vector<std::vector<std::vector<int>>> grid(
        nsZ, std::vector<std::vector<int>>(nsY, std::vector<int>(nsX, 0)));

    std::default_random_engine generator(geometrySeed);
    std::uniform_int_distribution<int> distributionX(pX, pX + nX - 1);
    std::uniform_int_distribution<int> distributionY(pY, pY + nY - 1);
    std::uniform_int_distribution<int> distributionZ(pZ, pZ + nZ - 1);

    std::uniform_int_distribution<int> distributionPathFluctuation(-1, 1);

    std::vector<int> tags(radiusVector.size());
    for (int i = 0; i < radiusVector.size(); i++)
        tags[i] = i + 1; // radiusVector[i]; //i+1;

    std::vector<long> count(radiusVector.size() + 1, 0);
    std::vector<long> countTarget(radiusVector.size() + 1, 0);
    count[0] = static_cast<long>(nX) * nY * nZ;

    myReal porosity =
        1.0 - tubeVolumeFraction; // porosity is 1 - tube volume fraction

    countTarget[0] = porosity * count[0];
    for (int i = 0; i < countTarget.size(); i++)
        countTarget[i + 1] =
            volumeFractions[i] * static_cast<long>(nX) * nY * nZ;

    // Free pass
    for (int entity = 0; entity < std::max(1, int(radiusVector.size() / 1));
         entity++) {
        int y, z; // define y and z outside the loop to make tubes along x
                  // direction

        long maxIterations = static_cast<long>(
            static_cast<double>(clusteringFactor) *
            static_cast<double>(maxNumberOfSpheres));

        for (long i = 0; i < maxIterations; i++) {
            if (count[entity + 1] < countTarget[entity + 1]) {

                int radius = std::ceil(radiusVector[entity]);

                int x = pX + (i % nX);

                if (x == pX) // new y and z for each sweep along x
                {
                    y = distributionY(generator);
                    z = distributionZ(generator);
                } else {
                    y += distributionPathFluctuation(generator) * (i % 2);
                    z += distributionPathFluctuation(generator) * ((i + 1) % 2);
                }

                // make y and z periodic with modulo in case they go out of
                // bounds
                y = (y - pY + nY) % nY + pY;
                z = (z - pZ + nZ) % nZ + pZ;

                // check if y and z are within bounds
                if ((y < pY) or (y >= pY + nY) or (z < pZ) or (z >= pZ + nZ)) {

                } else {

                    if (grid[z][y][x] == 0) {
                        count[0] -= 1;
                        grid[z][y][x] = tags[entity];
                        count[1 + entity] += 1;

                        for (int zz = z - radius; zz <= z + radius; zz++)
                            for (int yy = y - radius; yy <= y + radius; yy++)
                                for (int xx = x - radius; xx <= x + radius;
                                     xx++) {
                                    if (grid[zz][yy][xx] == 0) {
                                        if (std::sqrt(
                                                (xx - x) * (xx - x) +
                                                (yy - y) * (yy - y) +
                                                (zz - z) * (zz - z)) <=
                                            radiusVector[entity]) {
                                            count[0] -= 1;
                                            grid[zz][yy][xx] = tags[entity];
                                            count[1 + entity] += 1;
                                        }
                                    }
                                }
                    }
                }
            }
        }
        if (count[0] <= countTarget[0]) {
            std::cout << "Target porosity reached during structure generation. "
                         "Stopping and continuing to write to file"
                      << std::endl;
            break;
        }
    }

    // Pass with proximity check
    long maxIterations = static_cast<long>(
        static_cast<double>(1.0) * static_cast<double>(maxNumberOfSpheres));

    for (int entity = 0; entity < radiusVector.size(); entity++) {
        int y, z; // define y and z outside the loop to make tubes along x
                  // direction

        for (long i = 0; i < maxIterations; i++) {

            if (count[entity + 1] < countTarget[entity + 1]) {

                int radius = std::ceil(radiusVector[entity]);

                // x sweeps along the X direction to make multiple tubes in
                // passes Algorithm: modulo if i/nX
                int x = pX + (i % nX);

                if (i % nX == 0) // new y and z for each sweep along x
                {
                    y = distributionY(generator);
                    z = distributionZ(generator);
                } else {
                    y += distributionPathFluctuation(generator) * (i % 2);
                    z += distributionPathFluctuation(generator) * ((i + 1) % 2);
                }

                // make y and z periodic with modulo in case they go out of
                // bounds
                y = (y - pY + nY) % nY + pY;
                z = (z - pZ + nZ) % nZ + pZ;

                // check if y and z are within bounds
                if ((y < pY) or (y >= pY + nY) or (z < pZ) or (z >= pZ + nZ)) {

                } else {

                    if (grid[z][y][x] == 0) {
                        int connectivityChecked =
                            0; // to check if connectivity is checked for this
                               // point

                        for (int zz = z - radius; zz <= z + radius; zz++)
                            for (int yy = y - radius; yy <= y + radius; yy++)
                                for (int xx = x - radius; xx <= x + radius;
                                     xx++) {

                                    if (grid[zz][yy][xx] == 0) {

                                        if ((grid[z - (zz - z)][y - (yy - y)]
                                                 [x - (xx - x)] != 0) or
                                            (connectivityChecked ==
                                             1)) // check if the point is
                                                 // connected to any other point
                                        {

                                            connectivityChecked = 1;

                                            if (std::sqrt(
                                                    (xx - x) * (xx - x) +
                                                    (yy - y) * (yy - y) +
                                                    (zz - z) * (zz - z)) <=
                                                radiusVector[entity]) {
                                                count[0] -= 1;
                                                grid[zz][yy][xx] = tags[entity];
                                                count[1 + entity] += 1;

                                                if (grid[z][y][x] == 0) {
                                                    count[0] -= 1;
                                                    grid[z][y][x] =
                                                        tags[entity];
                                                    count[1 + entity] += 1;
                                                }
                                            }
                                        }
                                    }
                                }
                    }
                }
            }
        }
        if (count[0] <= countTarget[0]) {
            std::cout << "Target porosity reached during structure generation. "
                         "Stopping and continuing to write to file"
                      << std::endl;
            break;
        }
    }

    if (count[0] > countTarget[0]) {
        std::cout << "Target porosity not reached during structure generation. "
                     "Continuing to bridge spheres"
                  << std::endl;

        // Do bridging if target porosity is not reached
        for (int entity = radiusVector.size() - 1; entity > 0; entity--) {
            if (count[0] <= countTarget[0]) break;
            std::cout << "Bridging entity " << entity << " with radius "
                      << radiusVector[entity] << std::endl;

            for (int i = 0; i < int(0.01 * maxNumberOfSpheres); i++) {

                if (count[0] <= countTarget[0]) {
                    break;
                }

                int radius = std::ceil(radiusVector[entity]);

                int x = distributionX(generator);
                int y = distributionY(generator);
                int z = distributionZ(generator);

                if (grid[z][y][x] == 0) {
                    int connectivityChecked =
                        0; // to check if connectivity is checked for this point

                    for (int zz = z - radius; zz <= z + radius; zz++)
                        for (int yy = y - radius; yy <= y + radius; yy++)
                            for (int xx = x - radius; xx <= x + radius; xx++) {
                                if (grid[zz][yy][xx] == 0) {
                                    if ((grid[z - (zz - z)][y - (yy - y)]
                                             [x - (xx - x)] != 0) or
                                        (connectivityChecked ==
                                         1)) // check if the point is connected
                                             // to any other point
                                    {
                                        connectivityChecked = 1;

                                        if (std::sqrt(
                                                (xx - x) * (xx - x) +
                                                (yy - y) * (yy - y) +
                                                (zz - z) * (zz - z)) <=
                                            radiusVector[entity]) {
                                            count[0] -= 1;
                                            grid[zz][yy][xx] = tags[entity];
                                            count[1 + entity] += 1;

                                            if (grid[z][y][x] == 0) {
                                                count[0] -= 1;
                                                grid[z][y][x] = tags[entity];
                                                count[1 + entity] += 1;
                                            }
                                        }
                                    }
                                }
                            }
                }
            }
        }
    }

    if (count[0] > countTarget[0]) {
        std::cout << "WARNING!!! Target porosity could not be achieved"
                  << std::endl;
    } else {
        std::cout << "Target porosity reached during bridging. Stopping and "
                     "continuing to write to file"
                  << std::endl;
    }

    // marbles compatibility: make tubes 0 and solid space 1
    if (forMARBLES_LBM_code) {
        for (int k = 0; k < grid.size(); k++)
            for (int j = 0; j < grid[0].size(); j++)
                for (int i = 0; i < grid[0][0].size(); i++) {
                    // Make tag > 0 as 0 and 0 as 1
                    if (grid[k][j][i] > 0)
                        grid[k][j][i] = 0;
                    else
                        grid[k][j][i] = 1;
                }

        // Mark the first physical X plane as fluid
        for (int k = 0; k < grid.size(); k++)
            for (int j = 0; j < grid[0].size(); j++) grid[k][j][pX] = 0;

        // Mark the last physical X plane as fluid
        for (int k = 0; k < grid.size(); k++)
            for (int j = 0; j < grid[0].size(); j++)
                grid[k][j][grid[0][0].size() - 1 - pX] = 0;
    }

    std::filesystem::create_directories("results");

    std::ofstream ofile;
    std::string filename = "./results/microstructure_nX" + std::to_string(nX) +
                           "_nY" + std::to_string(nY) + "_nZ" +
                           std::to_string(nZ);

    filename += "_seed" + std::to_string(geometrySeed) + ".csv";

    ofile.open(filename);

    // Build the entire output as a string first
    std::string output;
    output.reserve(1000000000); // Reserve space for ~1GB string

    output += "X,Y,Z,tag\n";

    for (int k = pZ; k < grid.size() - pZ; k++) {
        std::cout << 100.0 * (k - pZ) / (grid.size() - 1) << " %" << std::endl;
        for (int j = pY; j < grid[0].size() - pY; j++)
            for (int i = pX; i < grid[0][0].size() - pX; i++) {
                output += std::to_string(i - pX) + "," +
                          std::to_string(j - pY) + "," +
                          std::to_string(k - pZ) + "," +
                          std::to_string(grid[k][j][i]) + "\n";
            }
    }

    // Write everything at once
    std::cout << "Writing to file: " << filename << std::endl;
    ofile << output;
    ofile.close();
    std::cout << "File written: " << filename << std::endl;

    // Also write binary data for faster TIFF conversion
    std::string binFilename = "./results/microstructure_nX" +
                              std::to_string(nX) + "_nY" + std::to_string(nY) +
                              "_nZ" + std::to_string(nZ) + "_seed" +
                              std::to_string(geometrySeed) + ".bin";
    std::ofstream binFile(binFilename, std::ios::binary);

    std::cout << "Writing binary file for TIFF conversion..." << std::endl;
    for (int k = pZ; k < grid.size() - pZ; k++) {
        if (k % 100 == 0)
            std::cout << "Binary write progress: "
                      << 100.0 * (k - pZ) / (grid.size() - 2 * pZ) << " %"
                      << std::endl;
        for (int j = pY; j < grid[0].size() - pY; j++)
            for (int i = pX; i < grid[0][0].size() - pX; i++) {
                uint16_t value = static_cast<uint16_t>(grid[k][j][i]);
                binFile.write(
                    reinterpret_cast<const char*>(&value), sizeof(uint16_t));
            }
    }
    binFile.close();
    std::cout << "Binary file written: " << binFilename << std::endl;

    // Write dimensions file for Python script
    std::ofstream dimFile("./results/dimensions.txt");
    dimFile << nX << " " << nY << " " << nZ << std::endl;
    dimFile.close();

    std::cout << "Open with Paraview with CSV Reader > Apply" << std::endl;
    std::cout << "option + space > Table to Structured Grid" << std::endl;
    std::cout << "Whole Extent: 0 to " << nX - 1 << std::endl;
    std::cout << "              0 to " << nY - 1 << std::endl;
    std::cout << "              0 to " << nZ - 1 << std::endl;
    std::cout << "X Column X" << std::endl;
    std::cout << "Y Column Y" << std::endl;
    std::cout << "Z Column Z" << std::endl;

    auto run_with_uvx = [](const std::string& script,
                           const std::string& uvx_packages) -> int {
        std::string uvx_command = "uvx " + uvx_packages + " python " + script;
        int status = std::system(uvx_command.c_str());
        if (status != 0) {
            std::cout << "uvx command failed (status " << status
                      << "). Falling back to `python " << script << "`."
                      << std::endl;
            status = std::system(("python " + script).c_str());
            if (status != 0) {
                std::cout << "`python` executable not available or script "
                             "failed (status "
                          << status << "). Attempting `python3` fallback."
                          << std::endl;
                status = std::system(("python3 " + script).c_str());
            }
        }
        return status;
    };

    const std::string binary_uvx_packages =
        "--with numpy --with tifffile --with imagecodecs";
    const std::string cad_uvx_packages =
        "--with numpy --with tifffile --with imagecodecs --with scikit-image "
        "--with trimesh";

    if (std::filesystem::exists("binaryToTiff.py")) {
        std::cout << "Running Python TIFF conversion (via uvx)..." << std::endl;
        int result = run_with_uvx("binaryToTiff.py", binary_uvx_packages);
        if (result == 0) {
            std::cout << "TIFF conversion completed successfully!" << std::endl;

            if (std::filesystem::exists("tiffToCAD.py")) {
                std::cout << "Running Python CAD conversion (via uvx)..."
                          << std::endl;
                int cadResult = run_with_uvx("tiffToCAD.py", cad_uvx_packages);
                if (cadResult == 0) {
                    std::cout << "CAD conversion completed successfully!"
                              << std::endl;
                    std::cout << "STL files saved in cad_exports/ directory"
                              << std::endl;
                } else {
                    std::cout << "CAD conversion failed. Run manually: uvx "
                              << cad_uvx_packages << " python tiffToCAD.py"
                              << std::endl;
                }
            } else {
                std::cout << "tiffToCAD.py not found. To convert to CAD "
                             "formats, run: uvx "
                          << cad_uvx_packages << " python tiffToCAD.py"
                          << std::endl;
            }

        } else {
            std::cout << "TIFF conversion failed. Run manually: uvx "
                      << binary_uvx_packages << " python binaryToTiff.py"
                      << std::endl;
        }
    } else {
        std::cout << "Python script not found. To convert to TIFF, run: uvx "
                  << binary_uvx_packages << " python binaryToTiff.py"
                  << std::endl;
    }

    return 0;
}
