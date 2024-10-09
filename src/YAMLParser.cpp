#include "YAMLParser.h"
#include "Utils.h"
#include <yaml-cpp/yaml.h>
#include <vector>
#include <string>

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