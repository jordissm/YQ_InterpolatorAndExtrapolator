#include <iostream> // std::cout, std::endl
#include <string>   // std::string
#include <numeric>  // std::iota
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <filesystem>
#include <mutex>
#include <yaml-cpp/yaml.h>

#ifdef OPENMP_AVAILABLE
    #define OPENMP_ENABLED 1
    #include <omp.h>
#else
    #define OPENMP_ENABLED 0
#endif

using namespace std::chrono;


constexpr size_t outputWidth = 55;
std::mutex mtx;

struct ExperimentalData {
    std::vector<std::pair<double, std::pair<double, std::pair<double,double>>>> datapoints;
};

struct InputParameters {
    std::string yieldName;
    double intervalMin;
    double intervalMax;
    double delta;
    int order;
};

void display_progress(float progress) {
    int barWidth = outputWidth-8;
    std::cerr << "[";
    int pos = static_cast<int>(barWidth * progress);
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cerr << "=";
        else if (i == pos) std::cerr << ">";
        else std::cerr << " ";
    }
    std::cerr << "] " << int(progress * 100.0) << " %\r";
    std::cerr.flush();
}

void PrintUpdateText(const std::string& text, const size_t& lineWidth, const std::string& type = "info") {
    std::string firstLinePrefix = (type == "info") ? " > " : " ! ";
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

    if (type == "info") {
        std::cout << wrappedText.str() << std::endl;;
    } else if (type == "error") {
        std::cerr << wrappedText.str() << std::endl;;
    }
}

std::vector<int> FindRepresentation(long long n, const int& base, const int& fixed_size) {
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

std::string FormatNumberWithThousandsSeparator(const long long& n) {
    std::string result = std::to_string(n); // Convert integer 'n' to a string and store it in 'result'

    // Loop to insert commas for thousands separators in the string
    for (int i = result.size() - 3; i > 0; i -= 3) {
        result.insert(i, ","); // Insert a comma after every three digits from the end
    }

    return result; // Return the formatted string
}

std::string FormatNumber(const double& number, const int& precision) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(precision) << number;
    return stream.str();
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

std::vector<std::vector<double>> FindEstimationParameters(const ExperimentalData& experimentalData, const double& intervalMin, const double& intervalMax, const double& delta, const int& order = 1) {
    auto start = high_resolution_clock::now();

    // Generate a vector with all the possible values for the polynomial coefficients
    std::vector<double> possibleCoefficients;
    for (double i = intervalMin; i <= intervalMax; i += delta) {
        possibleCoefficients.push_back(i);
    }

    int n_params = order + 1;
    int base = possibleCoefficients.size();

    std::vector<std::vector<double>> estimatedParameterSet;
    long long n_params_combinations = static_cast<long long>(std::exp(n_params * std::log(base)));
    PrintUpdateText("Number of possible combinations: " + FormatNumberWithThousandsSeparator(n_params_combinations), outputWidth);

    #if OPENMP_ENABLED
        long num_items_computed = 0;
        int n_threads = omp_get_max_threads();
        #pragma omp parallel for schedule(dynamic,100) reduction(+:num_items_computed) num_threads(n_threads)
    #endif
    for (long long i = 0; i < n_params_combinations; i++) {
        std::vector<int> parameterRepresentation = FindRepresentation(i, base, n_params);
        std::vector<double> testParameters(n_params);
        for (int j = 0; j < n_params; ++j) {
            testParameters[j] = possibleCoefficients[parameterRepresentation[j]];
        }

        bool isValid = CheckConditions(testParameters, experimentalData, order) && (order <= 1 || testParameters[n_params - 1] != 0);

        
        if (isValid) {
            std::lock_guard<std::mutex> guard(mtx);
            estimatedParameterSet.push_back(testParameters);
        }

        #if OPENMP_ENABLED
            #pragma omp atomic
                num_items_computed++;
            if (i % (n_params_combinations / 100) == 0) {
                #pragma omp critical
                if (omp_get_thread_num() == 0) {
                    float progress = static_cast<float>(n_threads) * num_items_computed / n_params_combinations;
                    display_progress(progress);
                }
            }
        #endif
    }
    #if OPENMP_ENABLED
        display_progress(1.0);
        std::cerr << std::endl;
    #endif

    long long n_valid_combinations = estimatedParameterSet.size();
    PrintUpdateText("Number of valid combinations: " + FormatNumberWithThousandsSeparator(n_valid_combinations) + " (" + FormatNumber(100.0 * n_valid_combinations / n_params_combinations, 2) + "%)", outputWidth);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    PrintUpdateText("Time taken: " + FormatNumber(duration.count(), 0) + " ms", outputWidth);

    return estimatedParameterSet;
}

std::tuple<double, double> EstimateExtrapolatedYieldRatioFromYQ(double yq, const std::vector<std::vector<double>>& estimatedParameterSet, int order) {
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

    return {mean, error};
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
void PrintCenteredText(const std::string& text, const size_t& width) {
    std::vector<std::string> lines = SplitIntoLines(text);
    
    for (const std::string& line : lines) {
        size_t padding = (width > line.length()) ? (width - line.length()) / 2 : 0;
        std::cout << std::string(padding, ' ') << line << "\n";
    }
}

// Function to print the program title and useful information
void PrintProgramInfo() {
    const std::string description = "This program provides a polynomial interpolation and extrapolation in $Y_Q$ for different particle yield ratios in heavy-ion collisions.";
    std::cout << "\n";
    std::cout << "**************************************************\n";
    std::cout << "**   Polynomial interpolator/extrapolator for   **\n";
    std::cout << "**   yield ratios in heavy-ion collisions.      **\n";
    std::cout << "**                                              **\n";
    std::cout << "**   Jordi Salinas San Martin,                  **\n";
    std::cout << "**   U. of I. Urbana-Champaign, 2024            **\n";
    std::cout << "**************************************************\n";
    std::cout << "\n";
    PrintUpdateText(description, outputWidth);
}

void PrintSectionTitle(const std::string& title) {
    const std::string separator(outputWidth, '=');
    std::cout << "\n";
    std::cout << separator << std::endl;
    PrintCenteredText(title, outputWidth);
    std::cout << separator << std::endl;
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

void PrintInputParameters(const InputParameters& inputParameters) {
    std::vector<std::pair<std::string, double>> parameters = {
        {"Interval minimum", inputParameters.intervalMin},
        {"Interval maximum", inputParameters.intervalMax},
        {"Delta", inputParameters.delta},
        {"Order", inputParameters.order}
    };
    std::string message = "The following input parameters have been provided:";
    PrintUpdateText(message, outputWidth);
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
    int title_size = inputParameters.yieldName.size();
    std::string title_rightwhitespace((header.size() - title_size - 5), ' ');
    std::string title = " * " + inputParameters.yieldName + " *" + title_rightwhitespace;
    PrintCenteredText(title, outputWidth);
    PrintCenteredText(separator, outputWidth);
    PrintCenteredText(header, outputWidth);
    PrintCenteredText(separator, outputWidth);

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

        PrintCenteredText(column_1_text + column_2_text, outputWidth);
    }
    PrintCenteredText(separator, outputWidth);
    std::cout << std::endl;
}

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

void displaySettings() {
    #if OPENMP_ENABLED
        PrintUpdateText("OpenMP enabled: True", 50);
    #else
        PrintUpdateText("OpenMP enabled: False", 50);
    #endif
}

struct Approximation {
    int order;
    std::vector<int> interval;
    double delta;
};

struct Yield {
    std::string name;
    bool isTurnedOn;
    std::string filename;
    std::vector<Approximation> approximations;
    ExperimentalData data;
};

void parseYAML(const std::string& filename, std::vector<Yield>& yields) {
    YAML::Node config = YAML::LoadFile(filename);

    for (const auto& yieldNode : config["yields"]) {
        Yield yield;
        yield.name = yieldNode["name"].as<std::string>();
        yield.isTurnedOn = yieldNode["turned_on"].as<bool>();
        yield.filename = yieldNode["filename"].as<std::string>();

        for (const auto& approxNode : yieldNode["approximations"]) {
            Approximation approx;
            approx.order = approxNode["order"].as<int>();
            approx.interval = approxNode["interval"].as<std::vector<int>>();
            approx.delta = approxNode["delta"].as<double>();
            yield.approximations.push_back(approx);
        }

        yield.isTurnedOn ? yields.push_back(yield) : void();
    }
}



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
                    auto [mean, error] = EstimateExtrapolatedYieldRatioFromYQ(yq, estimationParametersCollection[yield.name][approx.order], approx.order);
                    OUT << "\n";
                    OUT << std::setw(5)  << std::setfill(' ') << std::left  << std::setprecision(3) << std::fixed           << yq         << ' ';
                    OUT << std::setw(15) << std::setfill(' ') << std::right << std::scientific      << std::setprecision(8) << mean       << ' ';
                    OUT << std::setw(15) << std::setfill(' ') << std::right << std::scientific      << std::setprecision(8) << mean-error << ' ';
                    OUT << std::setw(15) << std::setfill(' ') << std::right << std::scientific      << std::setprecision(8) << mean+error << ' ';
                    OUT << std::setw(7)  << std::setfill(' ') << std::right << approx.order         << ' ';
                    OUT << std::setw(15) << std::setfill(' ') << std::right << approximationType;
                }
            }
            OUT.close();
        }
    }
 
    return 0;
}