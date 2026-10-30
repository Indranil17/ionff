// Analytic Na+ / solvent-site interaction terms, header only.
//
// Units: distances in Angstrom, energies in kJ/mol, charges in e, polarisabilities
// in Angstrom^3. The functional forms and constants match src/ionff/terms.py, and
// cxx/tests/test_terms.cpp checks them against values written by the Python code.
#pragma once

#include <cmath>
#include <cstddef>

namespace ionff {

constexpr double K_COUL = 1389.35457520;  // kJ mol^-1 Angstrom e^-2
constexpr double Q_NA = 1.0;

// Tang-Toennies damping, 1 - exp(-x) * sum_{k<=n} x^k / k!
inline double tt_damp(double x, int n) {
    double term = 1.0, total = 1.0;
    for (int k = 1; k <= n; ++k) {
        term *= x / k;
        total += term;
    }
    return 1.0 - std::exp(-x) * total;
}

inline double e_elst(double r, double q, double beta) {
    return K_COUL * Q_NA * q / r * (1.0 - std::exp(-beta * r) * (1.0 + beta * r));
}

inline double e_exch(double r, double logA, double B) { return std::exp(logA - B * r); }

inline double e_ind(double r, double alpha, double c) {
    const double r2 = r * r;
    return -0.5 * K_COUL * Q_NA * Q_NA * alpha / (r2 * r2) * tt_damp(c * r, 4);
}

inline double e_disp(double r, double C6, double d) {
    const double r2 = r * r;
    return -C6 / (r2 * r2 * r2) * tt_damp(d * r, 6);
}

// Parameters of one solvent site type.
struct SiteParams {
    double q, beta, logA, B, alpha, c, C6, d;
};

// The four terms and their sum for one Na+ to site distance.
struct Energy {
    double elst = 0.0, exch = 0.0, ind = 0.0, disp = 0.0;
    double total() const { return elst + exch + ind + disp; }
    Energy& operator+=(const Energy& o) {
        elst += o.elst;
        exch += o.exch;
        ind += o.ind;
        disp += o.disp;
        return *this;
    }
};

inline Energy site_energy(double r, const SiteParams& p) {
    return {e_elst(r, p.q, p.beta), e_exch(r, p.logA, p.B), e_ind(r, p.alpha, p.c), e_disp(r, p.C6, p.d)};
}

// Pair-additive energy of Na+ with n sites at distances r[i], using parameters p[i].
inline Energy total_energy(const double* r, const SiteParams* p, std::size_t n) {
    Energy sum;
    for (std::size_t i = 0; i < n; ++i) sum += site_energy(r[i], p[i]);
    return sum;
}

}  // namespace ionff
