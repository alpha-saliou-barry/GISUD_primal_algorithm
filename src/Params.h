#pragma once
#include <map>
#include <string>
#include <fstream>

class Params {
private:
	// Hash map between param name and param value
	std::map<std::string, std::string> params;

public:
	// Constructor of Params
	// Load params from path
	Params(std::string path);

	// Get the value associated to parameter name "parameter_name" as a double
	double getDoubleParam(std::string parameter_name);
};