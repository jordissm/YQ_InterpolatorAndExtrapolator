#include "PolynomialCoefficientsEstimation.h"

#include "Utils.h"
#include <iostream>
#include <numeric>
#include <vector>
#include <tuple>


#ifdef OPENMP_AVAILABLE
    #define OPENMP_ENABLED 1
    #include <omp.h>
#else
    #define OPENMP_ENABLED 0
#endif

using namespace std::chrono;


std::mutex mtx;



std::vector<int> FindRepresentation(long long n, const int& base, const int& fixed_size) {
    std::vector<int> representation(fixed_size, 0); // Initialize with zeros
    int index = fixed_size - 1;

    while (n > 0 && index >= 0) {
        representation[index] = n % base;
        n /= base;
        --index;
    }

    return representation;
}

bool CheckConditions(const std::vector<double>& testParameters, const ExperimentalData& experimentalData, const int& order) {
    for (const auto& datapoint : experimentalData.datapoints) {
        double yq = datapoint.first;
        double value = 0.0;
        double yqPower = 1.0;  // Start with yq^0

        for (int i = 0; i <= order; ++i) {
            value += testParameters[i] * yqPower;
            yqPower *= yq;  // Increment the power of yq
        }

        double lowerBound = datapoint.second.second.first;
        double upperBound = datapoint.second.second.second;
        if (!(lowerBound <= value && value <= upperBound)) {
            return false;
        }
    }
    return true;
}

std::vector<std::vector<double>> FindEstimationParameters(const ExperimentalData& experimentalData, const double& intervalMin, const double& intervalMax, const double& delta, const int& order) {
    if (intervalMin > intervalMax) {
        throw std::invalid_argument("intervalMin must be less than or equal to intervalMax");
    }
    if (delta <= 0) {
        throw std::invalid_argument("delta must be positive");
    }

    auto start = high_resolution_clock::now();

    // Generate a vector with all the possible values for the polynomial coefficients
    std::vector<double> possibleCoefficients;
    possibleCoefficients.reserve(static_cast<size_t>((intervalMax - intervalMin) / delta) + 1);

    for (double i = intervalMin; i <= intervalMax; i += delta) {
        possibleCoefficients.push_back(i);
    }

    int numCoefficients = order + 1;
    int numPossibleValues = possibleCoefficients.size();

    std::vector<std::vector<double>> estimatedParameterSet;
    long long numCombinations = static_cast<long long>(std::exp(numCoefficients * std::log(numPossibleValues)));

    PrintUpdateText("Number of possible combinations: " + FormatNumberWithThousandsSeparator(numCombinations), outputWidth);

    #if OPENMP_ENABLED
        long long numItemsComputed = 0;
        int numThreads = omp_get_max_threads();
        #pragma omp parallel for schedule(guided) reduction(+:numItemsComputed) num_threads(numThreads)
    #endif
    for (long long i = 0; i < numCombinations; i++) {
        std::vector<int> parameterRepresentation = FindRepresentation(i, numPossibleValues, numCoefficients);
        std::vector<double> testParameters(numCoefficients);
        for (int j = 0; j < numCoefficients; ++j) {
            testParameters[j] = possibleCoefficients[parameterRepresentation[j]];
        }

        bool isValid = CheckConditions(testParameters, experimentalData, order) && (order <= 1 || testParameters[numCoefficients - 1] != 0);

        if (isValid) {
            std::lock_guard<std::mutex> guard(mtx);
            estimatedParameterSet.push_back(testParameters);
        }

        #if OPENMP_ENABLED
            #pragma omp atomic
                numItemsComputed++;
            if (i % (numCombinations / 100) == 0) {
                #pragma omp critical
                if (omp_get_thread_num() == 0) {
                    float progress = static_cast<float>(numThreads) * numItemsComputed / numCombinations;
                    display_progress(progress);
                }
            }
        #endif
    }
    #if OPENMP_ENABLED
        display_progress(1.0);
        std::cerr << std::endl;
    #endif

    long long numValidCombinations = estimatedParameterSet.size();
    PrintUpdateText("Number of valid combinations: " + FormatNumberWithThousandsSeparator(numValidCombinations) + " (" + FormatNumber(100.0 * numValidCombinations / numCombinations, 2) + "%)", outputWidth);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    PrintUpdateText("Time taken: " + FormatNumber(duration.count(), 0) + " ms", outputWidth);

    return estimatedParameterSet;
}

std::vector<double> EstimateExtrapolatedYieldRatioFromYQ(double yq, const std::vector<std::vector<double>>& estimatedParameterSet, int order) {
    std::vector<double> values;
    for (const auto& params : estimatedParameterSet) {
        double value = 0.0;
        for (int i = 0; i <= order; i++) {
            value += params[i] * pow(yq, i);
        }
        values.push_back(value);
    }

    double lowerBound = *std::min_element(std::begin(values), std::end(values));
    double upperBound = *std::max_element(std::begin(values), std::end(values));
    double sum = std::accumulate(std::begin(values), std::end(values), 0.0);
    double mean = sum / values.size();
    // double mean = (lowerBound + upperBound) / 2.0;

    std::vector<double> result{};
    result.push_back(mean);
    result.push_back(mean - lowerBound);
    result.push_back(upperBound - mean);
    // return {mean, mean-lowerBound, upperBound-mean};
    return result;
}

