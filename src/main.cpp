#include <iostream> // std::cout, std::endl
#include <string>   // std::string
#include <numeric>  // std::iota
#include <complex>  // std::complex
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <thread>
#include <stdexcept>
#include <omp.h>

using namespace std::chrono;

void display_progress(float progress) {
    int barWidth = 50;
    std::cerr << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cerr << "=";
        else if (i == pos) std::cerr << ">";
        else std::cerr << " ";
    }
    std::cerr << "] " << int(progress * 100.0) << " %\r";
    std::cerr.flush();
}

std::vector<int> FindRepresentation(long long n, int base, int fixed_size) {
    std::vector<int> representation(fixed_size, 0); // Initialize with zeros
    int index = 0;

    while (n > 0 && index < fixed_size) {
        representation[index] = n % base;
        n /= base;
        index++;
    }

    std::reverse(representation.begin(), representation.end());
    return representation;
}

std::string FormatNumberWithThousandsSeparator(long long n) {
    std::string result = std::to_string(n); // Convert integer 'n' to a string and store it in 'result'

    // Loop to insert commas for thousands separators in the string
    for (int i = result.size() - 3; i > 0; i -= 3) {
        result.insert(i, ","); // Insert a comma after every three digits from the end
    }

    return result; // Return the formatted string
}

std::string FormatNumber(double number, int precision) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(precision) << number;
    return stream.str();
}


std::vector<std::vector<double>> FindEstimationParameters(const std::vector<std::pair<double, std::pair<double, double>>>& data, double intervalMin, double intervalMax, double delta, int order = 1) {
    auto start = high_resolution_clock::now();
    // Generate a vector with all the possible values for the polynomial coefficients
    std::vector<double> possibleCoefficients;
    for (double i = intervalMin; i <= intervalMax; i += delta) {
        possibleCoefficients.push_back(i);
    }

    int n_params = order + 1;
    int base = possibleCoefficients.size();

    auto CheckConditions = [&](const std::vector<double>& testParameters) {
        for (const auto& datapoint : data) {
            double yq = datapoint.first;
            double value = 0.0;
            for (int i = 0; i <= order; ++i) {
                value += testParameters[i] * std::pow(yq, i);
            }
            double lowerBound = datapoint.second.first;
            double upperBound = datapoint.second.second;
            if (!(lowerBound <= value && value <= upperBound)) {
                return false;
            }
        }
        return true;
    };

    std::vector<std::vector<double>> estimatedParameterSet;
    long long n_params_combinations = std::pow(base, n_params);
    std::cout << " > Number of possible combinations: " << FormatNumberWithThousandsSeparator(n_params_combinations) << std::endl;

    long num_items_computed = 0;
    int n_threads = omp_get_max_threads();
    #pragma omp parallel for schedule(dynamic,100) reduction(+:num_items_computed) num_threads(n_threads)
    for (long long i = 0; i < n_params_combinations; i++) {
        std::vector<int> parameterRepresentation = FindRepresentation(i, base, n_params);
        std::vector<double> testParameters(n_params);
        for (int j = 0; j < n_params; ++j) {
            testParameters[j] = possibleCoefficients[parameterRepresentation[j]];
        }
        if (order > 1) {
            if (CheckConditions(testParameters) && (testParameters[n_params-1] != 0)) {
            estimatedParameterSet.push_back(testParameters);
            }
        } else {
            if (CheckConditions(testParameters)) {
                estimatedParameterSet.push_back(testParameters);
            }
        }

        #pragma omp atomic
            num_items_computed++;
        if (i % (n_params_combinations / 100) == 0) {
            #pragma omp critical
            if (omp_get_thread_num() == 0) {
                float progress = (float)n_threads * num_items_computed / n_params_combinations;
                display_progress(progress);
            }
        }
    }
    display_progress(1.0);
    std::cerr << std::endl;

    long long n_valid_combinations = estimatedParameterSet.size();
    std::cout << " > Number of valid combinations: " << FormatNumberWithThousandsSeparator(n_valid_combinations) << " (" << FormatNumber(100.0*n_valid_combinations/n_params_combinations,2) << "%)" << std::endl;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << " > Time taken: " << duration.count() << " ms" << std::endl;
    return estimatedParameterSet;
}


std::tuple<double, double, double> EstimateExtrapolatedYieldRatioFromYQ(double yq, const std::vector<std::vector<double>>& estimatedParameterSet, int order) {
    std::vector<double> values;
    for (const auto& params : estimatedParameterSet) {
        double value = 0.0;
        for (int i = 0; i <= order; i++) {
            value += params[i] * pow(yq, i);
        }
        values.push_back(value);
    }

    double lowerBound = *std::min_element(values.begin(), values.end());
    double upperBound = *std::max_element(values.begin(), values.end());
    double mean = (lowerBound + upperBound) / 2.0;
    double error = upperBound - mean;

    return {mean, error, order};
}


std::pair<double, std::pair<double, double>> CombinedPrediction(const std::vector<std::tuple<double, double, double>>& predictions) {
    double mean = 0.0;
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();

    for (const auto& prediction : predictions) {
        mean += std::get<0>(prediction);
        min = std::min(min, std::get<1>(prediction));
        max = std::max(max, std::get<2>(prediction));
    }

    mean /= predictions.size();
    
    return {mean, {mean - min, max - mean}};
}


void PrintUpdateText(const std::string& text, size_t lineWidth) {
    const std::string firstLinePrefix = " > ";
    const std::string subsequentLinePrefix = "   ";
    size_t effectiveWidth = lineWidth - firstLinePrefix.length();
    std::istringstream words(text);
    std::string word;
    std::vector<std::string> lines;
    std::string line;
    bool isFirstLine = true;

    while (words >> word) {
        if (line.length() + word.length() + 1 > effectiveWidth) {
            if (isFirstLine) {
                lines.push_back(firstLinePrefix + line);
                effectiveWidth = lineWidth - subsequentLinePrefix.length();
                isFirstLine = false;
            } else {
                lines.push_back(subsequentLinePrefix + line);
            }
            line = word;
        } else {
            if (!line.empty()) {
                line += " ";
            }
            line += word;
        }
    }
    if (!line.empty()) {
        if (isFirstLine) {
            lines.push_back(firstLinePrefix + line);
        } else {
            lines.push_back(subsequentLinePrefix + line);
        }
    }

    std::ostringstream wrappedText;
    for (const auto& l : lines) {
        wrappedText << l << "\n";
    }

    std::cout << wrappedText.str() << std::endl;
}

std::vector<std::string> SplitIntoLines(const std::string& text) {
    std::istringstream stream(text);
    std::string line;
    std::vector<std::string> lines;
    while (std::getline(stream, line)) {
        lines.push_back(line);
    }
    return lines;
}

// Function to center a string within a given width
void PrintCenteredText(const std::string& text, size_t width) {
    std::vector<std::string> lines = SplitIntoLines(text);
    
    for (const std::string& line : lines) {
        size_t padding = (width > line.length()) ? (width - line.length()) / 2 : 0;
        std::cout << std::string(padding, ' ') << line << "\n";
    }
}

// Function to print the program title and useful information
void PrintProgramInfo() {
    const std::string description = "This program provides a polynomial interpolation and extrapolation in $Y_Q$ for different particle yield ratios in heavy-ion collisions.";

    std::cout << "************************************************\n";
    std::cout << "**  Polynomial interpolator/extrapolator for  **\n";
    std::cout << "**  yield ratios in heavy-ion collisions.     **\n";
    std::cout << "**                                            **\n";
    std::cout << "**  Jordi Salinas San Martin,                 **\n";
    std::cout << "**  U. of I. Urbana-Champaign, 2024           **\n";
    std::cout << "************************************************\n";
    std::cout << "\n";
    PrintUpdateText(description, 50);
}

void PrintSectionTitle(const std::string& title) {
    std::cout << "\n";
    std::cout << "================================================\n";
    PrintCenteredText(title,50);
    std::cout << "================================================\n";
    std::cout << "\n";
}


// Function to print a nicely formatted table
void PrintSystemsToInterpolate(const std::vector<double>& values) {

std::string message = "The following systems will be interpolated/extrapolated with the provided data points:";
PrintUpdateText(message, 50);
// Print the header
int whitespace_size = 3;
std::string whitespace(whitespace_size, ' ');
std::string column_1_header = "YQ";
std::string column_2_header = "interpolation/extrapolation";
std::string header = "|" + whitespace + column_1_header + whitespace + "|" 
                         + whitespace + column_2_header + whitespace + "|";
int column_1_width = column_1_header.size() + 2 * whitespace_size;
int column_2_width = column_2_header.size() + 2 * whitespace_size;
std::string separator = "+" + std::string(column_1_width, '-') + "+" + std::string(column_2_width, '-') + "+";
PrintCenteredText(separator, 50);
PrintCenteredText(header, 50);
PrintCenteredText(separator, 50);


// Print each value in the vector
for (const double& value : values) {
    int precision = 4;
    std::string rounded = std::to_string(std::round(value * std::pow(10, precision)) / std::pow(10, precision)).substr(0, std::to_string(std::round(value * std::pow(10, precision)) / std::pow(10, precision)).find(".") + precision + 1);
    int column_1_leftspace = (column_1_width - rounded.size()) / 2;
    int column_1_rightspace = column_1_width - rounded.size() - column_1_leftspace;
    std::string column_1_text = "|" + std::string(column_1_leftspace, ' ') + rounded + std::string(column_1_rightspace, ' ') + "|";
    std::string text = "interpolated";
    int column_2_leftspace = (column_2_width - text.size()) / 2;
    int column_2_rightspace = column_2_width - text.size() - column_2_leftspace;
    std::string column_2_text = std::string(column_2_leftspace, ' ') + text + std::string(column_2_rightspace, ' ') + "|";

    PrintCenteredText(column_1_text + column_2_text, 50);
}
PrintCenteredText(separator, 50);
}

void PrintInputParameters(const std::vector<std::pair<std::string, double>>& parameters) {
    std::string message = "The following input parameters have been provided:";
    PrintUpdateText(message, 50);
    // Print the header
    int whitespace_size = 3;
    std::string whitespace(whitespace_size, ' ');
    std::string column_1_header = "Parameter name";
    std::string column_2_header = "Value";
    std::string header = "|" + whitespace + column_1_header + whitespace + "|" 
                         + whitespace + column_2_header + whitespace + "|";
    int column_1_width = column_1_header.size() + 2 * whitespace_size;
    int column_2_width = column_2_header.size() + 2 * whitespace_size;
    std::string separator = "+" + std::string(column_1_width, '-') + "+" + std::string(column_2_width, '-') + "+";
    PrintCenteredText(separator, 50);
    PrintCenteredText(header, 50);
    PrintCenteredText(separator, 50);

    // Print each value in the vector
    for (const auto& parameter : parameters) {
        int precision = 4;  
        std::string rounded = std::to_string(std::round(parameter.second * std::pow(10, precision)) / std::pow(10, precision)).substr(0, std::to_string(std::round(parameter.second * std::pow(10, precision)) / std::pow(10, precision)).find(".") + precision + 1);
        int column_1_leftspace = (column_1_width - parameter.first.size()) / 2;
        int column_1_rightspace = column_1_width - parameter.first.size() - column_1_leftspace;
        std::string column_1_text = "|" + std::string(column_1_leftspace, ' ') + parameter.first + std::string(column_1_rightspace, ' ') + "|";
        int column_2_leftspace = (column_2_width - rounded.size()) / 2;
        int column_2_rightspace = column_2_width - rounded.size() - column_2_leftspace;
        std::string column_2_text = std::string(column_2_leftspace, ' ') + rounded + std::string(column_2_rightspace, ' ') + "|";

        PrintCenteredText(column_1_text + column_2_text, 50);
    }
    PrintCenteredText(separator, 50);
}

void ReadExperimentalData(const std::string& filename, std::vector<std::pair<double, std::pair<double, double>>>& data) {
    PrintUpdateText("Reading experimental data from files:", 50);
};


int main() {
    // Print useful information
    PrintProgramInfo();

    // Reading input parameters
    const std::string inputParametersTitle = "Input parameters";
    PrintSectionTitle(inputParametersTitle);
    // Parameters
    double intervalMin = -1000;
    double intervalMax = 1000;
    double delta = 10;
    int order = 3;
    std::vector<std::pair<std::string, double>> parameters = {
        {"Interval minimum", intervalMin},
        {"Interval maximum", intervalMax},
        {"Delta", delta},
        {"Order", order}
    };
    PrintInputParameters(parameters);

    // Reading experimental yield ratio data
    const std::string dataLoadingTitle = "Experimental data";
    PrintSectionTitle(dataLoadingTitle);
    std::vector<std::pair<double, std::pair<double, double>>> data = {
        {0.387, {1.005-0.056, 1.005+0.056}},
        {0.399, {1.015-0.051, 1.015+0.051}},
        {0.460, {1.003-0.050, 1.003+0.050}}
    };
    ReadExperimentalData("data.txt", data);

    // Reading in systems to predict
    const std::string interpolateSystemsTitle = "Systems to interpolate/extrapolate";
    PrintSectionTitle(interpolateSystemsTitle);
    std::vector<double> systemsToInterpolate = {0.417, 0.458, 0.5};
    PrintSystemsToInterpolate(systemsToInterpolate);

    // Find Estimation Parameters
    const std::string estimationParametersTitle = "Estimation Parameters";
    PrintSectionTitle(estimationParametersTitle);
    std::vector<std::vector<double>> estimatedParameterSet = FindEstimationParameters(data, intervalMin, intervalMax, delta, order);


    // Estimate Extrapolated Yield Ratio from YQ
    const std::string extrapolatedYieldRatioTitle = "Extrapolated Yield Ratio";
    PrintSectionTitle(extrapolatedYieldRatioTitle);

    std::vector<std::tuple<double, double, double>> predictions;
    for (int i = 0; i < systemsToInterpolate.size(); i++) {
        predictions.push_back(EstimateExtrapolatedYieldRatioFromYQ(systemsToInterpolate[i], estimatedParameterSet, order));
    }
    for (int i = 0; i < predictions.size(); i++) {
        std::cout << " > YQ = " << FormatNumber(systemsToInterpolate[i],3) << " : " << std::get<0>(predictions[i]) << " +/- " << std::get<1>(predictions[i]) << " (order = " << std::get<2>(predictions[i]) << ")" << std::endl;
    }


    // Combined Prediction
    const std::string combinedPredictionTitle = "Combined Prediction";
    PrintSectionTitle(combinedPredictionTitle);


    // Data export
    const std::string dataExportTitle = "Data Export";
    PrintSectionTitle(dataExportTitle);

 
    return 0;
}