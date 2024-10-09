#ifndef UTILS_H
#define UTILS_H

#include "ExperimentalData.h"
#include <string>
#include <vector>

// Global constant for the output width, used across multiple files for consistent formatting
constexpr size_t outputWidth = 55;

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

void PrintUpdateText(const std::string& text, const size_t& lineWidth, const std::string& type = "info");
void display_progress(float progress);
std::string FormatNumberWithThousandsSeparator(const long long& n);
std::string FormatNumber(const double& number, const int& precision);
std::vector<std::string> SplitIntoLines(const std::string& text);
void PrintCenteredText(const std::string& text, const size_t& width);
void PrintProgramInfo();
void PrintSectionTitle(const std::string& title);
void PrintSystemsToInterpolate(const std::vector<double>& values);
void displaySettings();

#endif // UTILS_H