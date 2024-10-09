#ifndef POLYNOMIALCOEFFICIENTSESTIMATION_H
#define POLYNOMIALCOEFFICIENTSESTIMATION_H

#include "ExperimentalData.h"
#include <mutex>
#include <vector>
#include <tuple>


std::vector<int> FindRepresentation(long long n, const int& base, const int& fixed_size);
bool CheckConditions(const std::vector<double>& testParameters, const ExperimentalData& experimentalData, const int& order);
std::vector<std::vector<double>> FindEstimationParameters(const ExperimentalData& experimentalData, const double& intervalMin, const double& intervalMax, const double& delta, const int& order);
std::vector<double> EstimateExtrapolatedYieldRatioFromYQ(double yq, const std::vector<std::vector<double>>& estimatedParameterSet, int order);

#endif // POLYNOMIALCOEFFICIENTSESTIMATION_H