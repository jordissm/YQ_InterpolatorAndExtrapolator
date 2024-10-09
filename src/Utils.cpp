#include "Utils.h"
#include <iostream>
#include <iomanip>
#include <sstream>

#ifdef OPENMP_AVAILABLE
    #define OPENMP_ENABLED 1
#else
    #define OPENMP_ENABLED 0
#endif

void PrintUpdateText(const std::string& text, const size_t& lineWidth, const std::string& type) {
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

void displaySettings() {
    #if OPENMP_ENABLED
        PrintUpdateText("OpenMP enabled: True", 50);
    #else
        PrintUpdateText("OpenMP enabled: False", 50);
    #endif
}
