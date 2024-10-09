#include "InputParameters.h"
#include "Utils.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>


void PrintInputParameters(const InputParameters& inputParameters) {
    std::vector<std::pair<std::string, double>> parameters;
    parameters.push_back({"Interval minimum", inputParameters.intervalMin});
    parameters.push_back({"Interval maximum", inputParameters.intervalMax});
    parameters.push_back({"Delta", inputParameters.delta});
    parameters.push_back({"Order", static_cast<double>(inputParameters.order)});
    
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