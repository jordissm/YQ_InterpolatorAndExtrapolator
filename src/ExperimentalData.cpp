#include "ExperimentalData.h"
#include "Utils.h"
#include <fstream>
#include <filesystem>
#include <sstream>

ExperimentalData ReadExperimentalData(const std::string& filename) {
    ExperimentalData experimentalData;
    PrintUpdateText("Reading in data from file: " + filename, outputWidth);

    std::filesystem::path completeFilename = "data/" + std::filesystem::path(filename).string();
    if (!std::filesystem::exists(completeFilename)) {
        PrintUpdateText("Error: File " + filename + " does not exist.", outputWidth, "error");
        return experimentalData;
    }

    std::ifstream infile(completeFilename);

    if (!infile.is_open()) {
        PrintUpdateText("Error opening file: " + filename, outputWidth, "error");
        return experimentalData;
    }

    std::string line;
    // Skip the header line
    std::getline(infile, line);

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        double YQ, value, min, max;
        if (!(iss >> YQ >> value >> min >> max)) {
            PrintUpdateText("Error reading line: " + line, outputWidth, "error");
            continue;  // Skip lines that cannot be parsed
        }

        experimentalData.datapoints.push_back({YQ, {value, {min, max}}});
    }

    infile.close();
    PrintUpdateText("Data loaded successfully", outputWidth);
    return experimentalData;
};