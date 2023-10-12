#include <adams_moulton.h>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <limits>

const double TOL = 1e-12;

constexpr double eps = std::numeric_limits<double>::epsilon();

TEST(AdamsMoulton, Exponential)
{
    size_t len = 1001;
    double dx = 0.001;
    std::vector<double> y(len), g(len, 1);
    y[0] = 1;
    y[1] = std::exp(dx);
    y[2] = std::exp(2*dx);

    LinearAdamsMoulton4 sol;
    sol.solve(y.begin() + 3, y.end(), g.begin() + 3, dx);


    std::ofstream outfile;
    outfile.open("adams_moulton_exponential.dat", std::ios::out);
    outfile << "# x  y \n";
    double x = 0;
    for(auto yi = y.begin(); yi != y.end(); yi++, x += dx){
        outfile << x << "  " << *yi << "\n";
    }

    ASSERT_LT(std::abs(y.back() - std::exp(1)), TOL);
}

TEST(AdamsMoulton, TwoExponential)
{
    size_t len = 1001;
    double dx = 0.001;
    std::vector<double> y(2*len), g(4*len, 0);
    std::transform(g.begin(), g.begin() + len, g.begin(),
            [](const double){return 1;});
    std::transform(g.begin() + 3*len, g.begin() + 4*len, g.begin() + 3*len,
            [](const double){return 1;});
    y[0] = 1;
    y[1] = std::exp(dx);
    y[2] = std::exp(2*dx);
    y[len] = 1;
    y[len + 1] = std::exp(dx);
    y[len + 2] = std::exp(2*dx);
    LinearAdamsMoulton4 sol;
    sol.solve(std::tuple{y.begin() + 3, y.begin() + 3 + len}, y.begin() + len, 
              std::tuple{g.begin() + 3        , g.begin() + 3 +   len, 
                         g.begin() + 3 + 2*len, g.begin() + 3 + 3*len}, dx);


    double x = 0;
    std::ofstream outfile;
    outfile.open("adams_moulton_two_exponential.dat", std::ios::out);
    outfile << "# x  y  y'  v \n";
    for(auto yi = y.begin(), gi = g.begin(); yi != y.begin() + len;
              yi++, gi++, x += dx){
        outfile << x << "  " << *yi << "  " << *(yi + len) << "  " << *gi << "  " << *(gi + len) << "  " << *(gi + 2*len) << "  " << *(gi + 3*len) << "\n";
    }

    ASSERT_LT(std::abs(y[len - 1] - std::exp(1)), TOL);
    ASSERT_LT(std::abs(y.back() - std::exp(1)), TOL);
}

TEST(AdamsMoulton, HArmonicOScillator)
{
    size_t n = 0;
    double exact_energy = n + 0.5;

    size_t len = 1001;
    double l = std::floor(std::sqrt(-2*std::log(eps)));
    double dx = l/(len - 1), j = -l/2;
    std::vector<double> y(2*len), g(4*len, 0), x(len), v(len);
    for(auto xi = x.begin(), vi = v.begin(); xi != x.end(); xi++, vi++, j += dx){
        *xi = j;
        *vi = std::pow(j, 2);
    }
    std::transform(g.begin() + len, g.begin() + 2*len, g.begin() + len,
            [](const double){return 1;});
    std::transform(v.begin(), v.end(), g.begin() + 2*len,
            [exact_energy](const double v){return v - 2*exact_energy;});

    y[0]       = std::pow(1/M_PI, 0.25)*std::exp(-std::pow(-l/2       , 2)/2);
    y[1]       = std::pow(1/M_PI, 0.25)*std::exp(-std::pow(-l/2 +   dx, 2)/2);
    y[2]       = std::pow(1/M_PI, 0.25)*std::exp(-std::pow(-l/2 + 2*dx, 2)/2);
    y[len - 3] = std::pow(1/M_PI, 0.25)*std::exp(-std::pow(-l/2 + 2*dx, 2)/2);
    y[len - 2] = std::pow(1/M_PI, 0.25)*std::exp(-std::pow(-l/2 +   dx, 2)/2);
    y[len - 1] = std::pow(1/M_PI, 0.25)*std::exp(-std::pow(-l/2       , 2)/2);
    
    y[len]       =    l/2        *y[0];
    y[len + 1]   = -(-l/2 +   dx)*y[1];
    y[len + 2]   = -(-l/2 + 2*dx)*y[2];
    y[2*len - 3] = -( l/2 - 2*dx)*y[len - 3];
    y[2*len - 2] = -( l/2 -   dx)*y[len - 2];
    y[2*len - 1] = -  l/2        *y[len - 1];

    LinearAdamsMoulton4 sol;
    auto i_f = sol.solve(std::tuple{y.begin() + 3, y.begin() + 3 + len}, y.begin() + len/2 + 1, 
              std::tuple{g.begin() + 3        , g.begin() + 3 +   len, 
                         g.begin() + 3 + 2*len, g.begin() + 3 + 3*len}, dx);
    double ydiff = *(std::get<0>(i_f) - 1);
    double ypdiff = *(std::get<1>(i_f) - 1);
    auto i_r = sol.solve(std::tuple{y.rbegin() + 3 + len, y.rbegin() + 3}, y.rbegin() + len + len/2 +1, 
              std::tuple{g.rbegin() + 3 + 3*len, g.rbegin() + 3 + 2*len, 
                         g.rbegin() + 3 +   len, g.rbegin() + 3         }, -dx);
    ydiff = std::abs(ydiff - *(std::get<0>(i_r) - 1));
    ypdiff = std::abs(ypdiff - *(std::get<1>(i_r) - 1));


    std::ofstream outfile;
    outfile.open("adams_moulton_quantum_harmonic_oscillator.dat", std::ios::out);
    outfile << "# x  y  y'  v  g11  g12  g21  g22\n";
    for(auto yi = y.begin(), gi = g.begin(), xi = x.begin(), vi = v.begin(); 
             yi != y.begin() + len;
              yi++, gi++, xi++, vi++){
        outfile << *xi << "  " << *yi << "  " << *(yi + len) << "  " << *vi << "  " << *gi << "  " << *(gi + len) << "  " << *(gi + 2*len) << "  " << *(gi + 3*len) << "\n";
    }

    ASSERT_LT(ydiff, TOL);
    ASSERT_LT(ypdiff, std::sqrt(TOL));
}
