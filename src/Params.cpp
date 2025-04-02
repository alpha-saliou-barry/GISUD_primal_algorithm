#include "Params.h"

// Constructor of Params
// Load params from path
Params::Params(std::string path) {
	std::ifstream file;
	file.open(path);

	while (file.good()) {
		std::string line;
		std::getline(file, line);

		std::string param_name;
		std::string param_value;

		bool isValue = false;
		for (int i = 0; i < line.size(); i++) {
			if (line[i] != ' ' && line[i] != '\n') {
				if (line[i] == ':') {
					isValue = true;
				}
				else {
					if (!isValue) {
						param_name += line[i];
					}
					else {
						param_value += line[i];
					}
				}
			}
		}

		params[param_name] = param_value;
	}
}


// Get the value associated to parameter name "parameter_name" as a double
double Params::getDoubleParam(std::string parameter_name) {
	return std::stod(params[parameter_name]);
}