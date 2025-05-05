#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <unordered_map>

class Settings {
public:
    explicit Settings(const std::string& filepath);

    std::string get_string(const std::string& key) const;
    int get_int(const std::string& key) const;
    double get_double(const std::string& key) const;
    bool get_bool(const std::string& key) const;

    // Utility for checking optional keys
    bool has_key(const std::string& key) const;

private:
    std::unordered_map<std::string, std::string> values;
    void parse_file(const std::string& filepath);
};

#endif // SETTINGS_H
