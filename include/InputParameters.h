#ifndef INPUT_PARAMETERS_H
#define INPUT_PARAMETERS_H

#include <string>

struct InputParameters {
    std::string yieldName;
    double intervalMin;
    double intervalMax;
    double delta;
    int order;
};

void PrintInputParameters(const InputParameters& inputParameters);

#endif // INPUT_PARAMETERS_H