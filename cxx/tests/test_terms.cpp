// Compare the C++ terms with reference values written by scripts/make_golden.py.
// Usage: test_terms <golden.csv>
// Each non-comment line: r,q,beta,logA,B,alpha,c,C6,d,elst,exch,ind,disp

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "ionff_terms.hpp"

namespace {

bool close(double a, double b) { return std::fabs(a - b) <= 1e-9 * std::fmax(1.0, std::fabs(b)); }

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s golden.csv\n", argv[0]);
        return 2;
    }
    std::ifstream in(argv[1]);
    if (!in) {
        std::fprintf(stderr, "cannot open %s\n", argv[1]);
        return 2;
    }
    std::string line;
    int rows = 0, bad = 0;
    std::vector<double> distances;
    std::vector<ionff::SiteParams> sites;
    ionff::Energy expected_sum;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string cell;
        double v[13];
        int n = 0;
        while (std::getline(ss, cell, ',') && n < 13) v[n++] = std::stod(cell);
        if (n != 13) {
            std::fprintf(stderr, "bad line: %s\n", line.c_str());
            return 2;
        }
        const ionff::SiteParams p{v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]};
        const ionff::Energy e = ionff::site_energy(v[0], p);
        const double got[4] = {e.elst, e.exch, e.ind, e.disp};
        for (int k = 0; k < 4; ++k) {
            if (!close(got[k], v[9 + k])) {
                std::fprintf(stderr, "row %d term %d: got %.12g expected %.12g\n", rows, k, got[k], v[9 + k]);
                ++bad;
            }
        }
        distances.push_back(v[0]);
        sites.push_back(p);
        expected_sum += ionff::Energy{v[9], v[10], v[11], v[12]};
        ++rows;
    }
    const ionff::Energy sum = ionff::total_energy(distances.data(), sites.data(), sites.size());
    if (!close(sum.total(), expected_sum.total())) {
        std::fprintf(stderr, "summed total: got %.12g expected %.12g\n", sum.total(), expected_sum.total());
        ++bad;
    }
    if (rows == 0) {
        std::fprintf(stderr, "no rows read\n");
        return 2;
    }
    std::printf("%d rows checked, %d mismatches\n", rows, bad);
    return bad == 0 ? 0 : 1;
}
