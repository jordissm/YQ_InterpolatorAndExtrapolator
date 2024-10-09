#ifndef YAML_PARSER_H
#define YAML_PARSER_H

#include <string>
#include <vector>

struct Yield;
void parseYAML(const std::string& filename, std::vector<Yield>& yields);

#endif // YAML_PARSER_H