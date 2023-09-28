#ifndef CSV_READER_H
#define CSV_READER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <string>
using namespace std;

class Csv_Reader {
private:
    map<double, double> data;

    map<double, double> initData(ifstream& s_file) {
        map <double, double> params;
        string line;
        while (getline(s_file, line)) {
            std::string delimiter = ",";
            std::string key;
            key = line.substr(0, line.find(delimiter));
            line.erase(0, line.find(delimiter) + 1);
            //istringstream is_line(line);
            //string key, value;
            //getline(is_line, key, ',');
            //getline(is_line, value);
            // size_t offset;
            // cout << stod(line, &offset) << "," << stod(line, &offset + 1) << endl;
            params.insert({string_to_double(key), string_to_double(line)});
        }
        return params;
    }

    double string_to_double(string& s) {
        stringstream ss;
        double d;
        ss << s;
        ss >> d;
        return d;
    }

public:
    Csv_Reader(const char* path) {
        read_new_file(path);
    }

    map<double, double> getData() {
        return data;
    }

    void read_new_file(const char* path) {
        ifstream file(path);
        data = initData(file);
    }

    void print_data() {
        for (auto itr = data.begin(); itr != data.end(); ++itr)
            cout << itr->first << " : " << itr->second << std::endl;
    }
};

#endif // CSV_READER_H
