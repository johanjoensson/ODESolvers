#include <numerov.h>
#include <cmath>
#include <gtest/gtest.h>
#include <vector>
#include <fstream>

const double TOL = 1e-12;

TEST(Numerov, QuantumHarmonicOscillator)
{
	int n = 2;
    auto exact_energy = n + 0.5;

    size_t len = 3001;
	double l = 15;
    double dl = l/static_cast<double>(len - 1);
	double j = -l/2;
	std::vector<double> v(len), g(len), x(len), zero(len, 0);
	for(auto vi = v.begin(), xi = x.begin(), gi = g.begin(); vi != v.end(); vi++, xi++, gi++, j += dl){
        *xi = j;
		*vi = std::pow(j, 2);
        *gi = 2*exact_energy - *vi;
	}
    std::vector<double> res(len);
    // Load initial conditions
    res[0] = std::exp(-std::pow(l/2, 2)/2);
    res[1] = std::exp(-std::pow(-l/2+dl, 2)/2);

    res[len - 2] = std::exp(-std::pow(-l/2+dl, 2)/2);
    res[len - 1] = std::exp(-std::pow(l/2, 2)/2);

    Numerov sol;
    auto i_meeting = sol.solve(res.begin() + 2, res.begin() + len/2 + 1, g.begin() + 2, zero.begin() + 2, dl);
    double dPsi = *(i_meeting - 1);
    auto ri_meeting = sol.solve(res.rbegin() + 2, res.rbegin() + len/2 + 1, g.rbegin() + 2, zero.rbegin() + 2, -dl);
    dPsi -= *(ri_meeting - 1);

    std::ofstream outfile;
    outfile.open("harmonic_oscillator_groundstate.dat", std::ios::out);
    outfile << "# x  y  v\n";
    for(auto xi = x.begin(), yi = res.begin(), vi = v.begin(); 
             xi != x.end(); xi++, yi++, vi++){
        outfile << *xi << "  " << *yi << "  " << *vi << "\n";
    }

    ASSERT_LT(std::abs(dPsi), TOL);
}
