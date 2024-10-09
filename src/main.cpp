#include "ExperimentalData.h"
#include "Utils.h"
#include "InputParameters.h"
#include "PolynomialCoefficientsEstimation.h"
#include "YAMLParser.h"
#include <iostream>
#include <string>
// #include <numeric>
#include <fstream>
#include <vector>
#include <tuple>
#include <unordered_map>
// #include <chrono>
#include <iomanip>
#include <algorithm>
#include <filesystem>

int main() {
    // Print useful information
    PrintProgramInfo();

    // Display settings
    displaySettings();

    // Reading input parameters
    const std::string inputParametersTitle = "Input parameters";
    PrintSectionTitle(inputParametersTitle);
    std::vector<Yield> yields;
    parseYAML("input.yml", yields);

    // Print the parsed data to verify
    for (const auto& yield : yields) {
        for (const auto& approx : yield.approximations) {
            InputParameters inputParameters;
            inputParameters.yieldName = yield.name;
            inputParameters.intervalMin = approx.interval[0];
            inputParameters.intervalMax = approx.interval[1];
            inputParameters.delta = approx.delta;
            inputParameters.order = approx.order;
            PrintInputParameters(inputParameters);
        }
    }

    // Reading experimental yield ratio data
    const std::string dataLoadingTitle = "Experimental data";
    PrintSectionTitle(dataLoadingTitle);
    for (auto& yield : yields) {
        yield.data = (yield.isTurnedOn) ? ReadExperimentalData(yield.filename) : ExperimentalData();
    }

    // Find Estimation Parameters
    const std::string estimationParametersTitle = "Obtaining valid parametrizations";
    PrintSectionTitle(estimationParametersTitle);
    std::unordered_map<std::string, std::unordered_map<int,std::vector<std::vector<double>>>> estimationParametersCollection;
    for (const auto& yield : yields) {
        if (yield.isTurnedOn && yield.data.datapoints.size() > 0){
            for (const auto& approx : yield.approximations) {
                PrintUpdateText("Estimating parameters for " + yield.name + " with order " + std::to_string(approx.order), outputWidth);
                std::vector<std::vector<double>> estimatedParameters = FindEstimationParameters(yield.data, approx.interval[0], approx.interval[1], approx.delta, approx.order);
                estimationParametersCollection[yield.name][approx.order] = estimatedParameters;
            }
        }
    }

    const std::string exportingDataTitle = "Exporting data";
    PrintSectionTitle(exportingDataTitle);

    for (const auto& yield : yields) {
        if (yield.isTurnedOn && yield.data.datapoints.size() > 0) {
            // Define output file path
            std::filesystem::path outputPath("output");
            std::filesystem::create_directory(outputPath);
            outputPath /= yield.filename;
            
            std::ofstream OUT(outputPath);
            if (!OUT.is_open()) {
                PrintUpdateText("Error writing to file: " + yield.filename, outputWidth, "error");
                continue;
            }
            double yqInterpolationLowerBound = yield.data.datapoints[0].first;
            double yqInterpolationUpperBound = yield.data.datapoints[yield.data.datapoints.size() - 1].first;
            OUT << std::setw(6)  << std::setfill(' ') << std::left << "YQ";
            OUT << std::setw(16) << std::setfill(' ') << std::right << "Mean_value ";
            OUT << std::setw(16) << std::setfill(' ') << std::right << "Min_value ";
            OUT << std::setw(16) << std::setfill(' ') << std::right << "Max_value ";
            OUT << std::setw(8)  << std::setfill(' ') << std::right << "Order ";
            OUT << std::setw(15) << std::setfill(' ') << std::right << "Approximation";
            for (const auto& approx : yield.approximations) {
                PrintUpdateText("Combining " + yield.name + " parametrizations of order " + std::to_string(approx.order) + " and exporting to file: " + "output/" + yield.filename, outputWidth);
                for (double yq = 0.300; yq <= 0.601; yq += 0.001) {    
                    std::string approximationType = (yq <= yqInterpolationLowerBound || yq >= yqInterpolationUpperBound) ? "extrapolation" : "interpolation";
                    auto result = EstimateExtrapolatedYieldRatioFromYQ(yq, estimationParametersCollection[yield.name][approx.order], approx.order);
                    double mean = result[0];
                    double delta_minus = result[1];
                    double delta_plus = result[2];
                    OUT << "\n";
                    OUT << std::setw(5)  << std::setfill(' ') << std::left  << std::setprecision(3) << std::fixed           << yq          << ' ';
                    OUT << std::setw(15) << std::setfill(' ') << std::right << std::scientific      << std::setprecision(8) << mean        << ' ';
                    OUT << std::setw(15) << std::setfill(' ') << std::right << std::scientific      << std::setprecision(8) << mean-delta_minus << ' ';
                    OUT << std::setw(15) << std::setfill(' ') << std::right << std::scientific      << std::setprecision(8) << mean+delta_plus  << ' ';
                    OUT << std::setw(7)  << std::setfill(' ') << std::right << approx.order         << ' ';
                    OUT << std::setw(15) << std::setfill(' ') << std::right << approximationType;
                }
            }
            OUT.close();
        }
    }
 
    return 0;
}