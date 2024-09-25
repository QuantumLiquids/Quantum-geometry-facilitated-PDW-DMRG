#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include "json.hpp"

using json = nlohmann::json;
using namespace std;

const int Ly = 4;
const int Lx = 48;
const int N = 2 * Lx * Ly;
const double ts = 1;
const double td = -1;
const double tsd_xy = 1;
const double tsd_nn = 0;
const double Uss = 8;
const double Udd = 8;
const double Usd = 8;
const int Hole = 192;
const vector<int> D_values = {21000};

const double doping = 1.0 / 2.0;
const double nbar = 1 - doping;

vector<double> kz_set = {0, M_PI};
vector<double> ky_set = {0, M_PI / 2, M_PI, 3 * M_PI / 2};
vector<double> kx_set = []() {
    vector<double> kx;
    for (double k = -M_PI; k <= M_PI; k += 0.1) kx.push_back(k);
    kx.push_back(0.0);
    sort(kx.begin(), kx.end());
    return kx;
}();

void SiteIdx2Coor3D(int Ly, int site_idx, int &x, int &y, int &z) {
    z = site_idx % 2;
    int site_idx_2d = site_idx / 2;
    y = site_idx_2d % Ly;
    x = site_idx_2d / Ly;
}

void calculate_nk(int D, vector<vector<vector<double>>> &nk) {
    string FileNamePostfix = to_string(Ly) + "x" + to_string(Lx) + "ts" + to_string(int(ts)) + "td" + to_string(int(td)) +
                             "tsd_xy" + to_string(int(tsd_xy)) + "tsd_nn" + to_string(int(tsd_nn)) + "Uss" + to_string(int(Uss)) +
                             "Udd" + to_string(int(Udd)) + "Usd" + to_string(int(Usd)) + "Hole" + to_string(Hole) + "D" + to_string(D) + ".json";

    cout << "Processing files with postfix: " << FileNamePostfix << endl;

    ifstream file0("../../data/single_particle_correlation0" + FileNamePostfix);
    ifstream file1("../../data/single_particle_correlation1" + FileNamePostfix);
    ifstream file2("../../data/single_particle_correlation2" + FileNamePostfix);
    ifstream file3("../../data/single_particle_correlation3" + FileNamePostfix);

    if (!file0.is_open()) {
        cerr << "Error: Unable to open file0: " << "../../data/single_particle_correlation0" + FileNamePostfix << endl;
        return;
    }
    if (!file1.is_open()) {
        cerr << "Error: Unable to open file1: " << "../../data/single_particle_correlation1" + FileNamePostfix << endl;
        return;
    }
    if (!file2.is_open()) {
        cerr << "Error: Unable to open file2: " << "../../data/single_particle_correlation2" + FileNamePostfix << endl;
        return;
    }
    if (!file3.is_open()) {
        cerr << "Error: Unable to open file3: " << "../../data/single_particle_correlation3" + FileNamePostfix << endl;
        return;
    }

    json data0, data1, data2, data3;
    file0 >> data0;
    file1 >> data1;
    file2 >> data2;
    file3 >> data3;

    int data_size = data0.size();

    vector<double> SingleParticleCorrelation(data_size);
    vector<array<int, 3>> distance(data_size);

    for (int i = 0; i < data_size; ++i) {
        int site1 = data0[i][0][0];
        int site2 = data0[i][0][1];

        SingleParticleCorrelation[i] = data0[i][1].get<double>() +
                                       data1[i][1].get<double>() -
                                       data2[i][1].get<double>() -
                                       data3[i][1].get<double>();

        int x1, y1, z1, x2, y2, z2;
        SiteIdx2Coor3D(Ly, site1, x1, y1, z1);
        SiteIdx2Coor3D(Ly, site2, x2, y2, z2);
        distance[i] = {x1 - x2, y1 - y2, z1 - z2};

        for (size_t kx_idx = 0; kx_idx < kx_set.size(); ++kx_idx) {
            double kx = kx_set[kx_idx];
            for (size_t ky_idx = 0; ky_idx < ky_set.size(); ++ky_idx) {
                double ky = ky_set[ky_idx];
                for (size_t kz_idx = 0; kz_idx < kz_set.size(); ++kz_idx) {
                    double kz = kz_set[kz_idx];
                    double phase = kx * distance[i][0] + ky * distance[i][1] + kz * distance[i][2];
                    nk[kx_idx][ky_idx][kz_idx] += cos(phase) * SingleParticleCorrelation[i];
                }
            }
        }
    }
}

void normalize_nk(vector<vector<vector<double>>> &nk) {
    for (auto &row : nk) {
        for (auto &col : row) {
            for (auto &val : col) {
                val = val / N + 1 - doping;
            }
        }
    }
}

int main() {
    for (int D : D_values) {
        vector<vector<vector<double>>> nk(kx_set.size(), vector<vector<double>>(Ly, vector<double>(2, 0)));
        calculate_nk(D, nk);
        normalize_nk(nk);

        // Dump nk data to a file
        string output_filename = "nkD" + to_string(D) + ".json";
        ofstream output(output_filename);
        if (!output.is_open()) {
            cerr << "Error: Unable to open output file: " << output_filename << endl;
            return 1;
        }
        output << json(nk);
        output.close();

        cout << "Data successfully written to " << output_filename << endl;
    }

    return 0;
}

