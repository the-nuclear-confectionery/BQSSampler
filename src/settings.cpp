#include "settings.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

Settings::Settings(const std::string& filepath) {
    parse_file(filepath);
}

void Settings::parse_file(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open())
        throw std::runtime_error("Unable to open settings file: " + filepath);

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key, value;
        if (iss >> key >> value) {
            values[key] = value;
        }
    }
}

std::string Settings::get_string(const std::string& key) const {
    auto it = values.find(key);
    if (it == values.end())
        throw std::runtime_error("Missing key in settings: " + key);
    return it->second;
}

int Settings::get_int(const std::string& key) const {
    return std::stoi(get_string(key));
}

double Settings::get_double(const std::string& key) const {
    return std::stod(get_string(key));
}

bool Settings::get_bool(const std::string& key) const {
    std::string val = get_string(key);
    std::transform(val.begin(), val.end(), val.begin(), ::tolower);
    return val == "true" || val == "1";
}

bool Settings::has_key(const std::string& key) const {
    return values.find(key) != values.end();
}
