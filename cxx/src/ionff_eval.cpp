// Evaluate the Na+ / site terms for distances read from stdin.
// Usage: ionff_eval q beta logA B alpha c C6 d < distances.txt
// One distance (Angstrom) per line in; prints r, elst, exch, ind, disp, total (kJ/mol).

#include <cstdio>
#include <cstdlib>

#include "ionff_terms.hpp"

int main(int argc, char** argv) {
    if (argc != 9) {
        std::fprintf(stderr, "usage: %s q beta logA B alpha c C6 d < distances.txt\n", argv[0]);
        return 2;
    }
    double v[8];
    for (int i = 0; i < 8; ++i) v[i] = std::atof(argv[i + 1]);
    const ionff::SiteParams p{v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]};
    std::printf("# r_A elst exch ind disp total (kJ/mol)\n");
    double r;
    while (std::scanf("%lf", &r) == 1) {
        const ionff::Energy e = ionff::site_energy(r, p);
        std::printf("%.4f %.6f %.6f %.6f %.6f %.6f\n", r, e.elst, e.exch, e.ind, e.disp, e.total());
    }
    return 0;
}
