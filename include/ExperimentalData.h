#ifndef EXPERIMENTAL_DATA_H
#define EXPERIMENTAL_DATA_H

#include <vector>
#include <string>

struct ExperimentalData {
    std::vector<std::pair<double, std::pair<double, std::pair<double,double>>>> datapoints;
};


ExperimentalData ReadExperimentalData(const std::string& filename);

#endif // EXPERIMENTAL_DATA_H