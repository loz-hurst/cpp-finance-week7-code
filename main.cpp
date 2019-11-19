/*
 * Code for week 7 exercises of C++ for Finance.
 *
 * Copyright 2019 Laurence Alexander Hurst
 *
 * This file is part of C++ for Finance.
 *
 *     C++ for Finance is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     C++ for Finance is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 * See the file LICENCE in the original source code repository for the
 * full licence.
 */

#include <ctime>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include "MonteCarlo.hpp"

namespace MC = MonteCarlo;

// Print out the MC simulation result - including value, error and calculation of time taken from (end-start)/(ticks/s)
void print_result(const double value, const double error, const time_t start, const time_t end) {
    std::cout << std::setprecision(8) << std::fixed
              << "Calculated estimate value of " << value
              << " with error " << error
              << " in " << ((end-start)/CLOCKS_PER_SEC) << "s"
              << std::endl;
}

// Convenience function for IS benchmark - sets the rate to mu and sigma to sigma in both options passed by value
void set_mu_sigma_pair(MC::Data & opt1, MC::Data & opt2, const double mu, const double sigma) {
    opt1.rate = mu; opt1.sigma = sigma;
    opt2.rate = mu; opt2.sigma = sigma;
}

int main() {
    // S_0, rate, sigma, maturity, strike, paths, steps
    const MC::Data option {100, 0.05, 0.20, 1, 100, 1000000, 100};

    std::cout << "Begin plain Monte Carlo:" << std::endl;
    const time_t plain_mc_start {std::clock()};
    const std::unique_ptr<MC::Result> plain_result {MC::Plain(option)};
    const time_t plain_mc_end {std::clock()};
    print_result(plain_result->value, plain_result->error, plain_mc_start, plain_mc_end);

    std::cout << "Begin ln_S Monte Carlo:" << std::endl;
    const time_t ln_S_mc_start {std::clock()};
    const std::unique_ptr<MC::Result> ln_S_result {MC::Ln_S(option)};
    const time_t ln_S_mc_end {std::clock()};
    print_result(ln_S_result->value, ln_S_result->error, ln_S_mc_start, ln_S_mc_end);

    std::cout << "Begin CV Monte Carlo:" << std::endl;
    const time_t cv_mc_start {std::clock()};
    double correlation;
    const std::unique_ptr<MC::Result> cv_result {MC::Cv(option, correlation)};
    const time_t cv_mc_end {std::clock()};
    print_result(cv_result->value, cv_result->error, cv_mc_start, cv_mc_end);
    std::cout << "Correlation was: " << correlation << std::endl;

    // Everything but paths and steps is irrelevant for the plain IS benchmark
    MC::Data is_bench_10k {0, 0, 0, 0, 0, 10000, 100};
    MC::Data is_bench_100k {0, 0, 0, 0, 0, 100000, 100};

    std::cout << "Begin IS benchmark Monte Carlo:" << std::endl;
    const time_t is_b_p_10k_mc_start {std::clock()};
    const std::unique_ptr<MC::Result> is_b_p_10k_result {MC::Is_Benchmark_Plain(is_bench_10k)};
    const time_t is_b_p_10k_mc_end {std::clock()};
    const time_t is_b_p_100k_mc_start {std::clock()};
    const std::unique_ptr<MC::Result> is_b_p_100k_result {MC::Is_Benchmark_Plain(is_bench_100k)};
    const time_t is_b_p_100k_mc_end {std::clock()};
    std::cout << "plain results:" << std::endl;
    std::cout << " 10k samples: "; print_result(is_b_p_10k_result->value, is_b_p_10k_result->error, is_b_p_10k_mc_start, is_b_p_10k_mc_end);
    std::cout << "100k samples: "; print_result(is_b_p_100k_result->value, is_b_p_100k_result->error, is_b_p_100k_mc_start, is_b_p_100k_mc_end);
    // Create a vector of pairs of mu, sigma values to try
    std::vector<std::pair<double,double>> is_b_trials {
        std::make_pair(1, 2), std::make_pair(2, 2), std::make_pair(2, 1), std::make_pair(2.3, 1), std::make_pair(2.3, 0.5)
    };
    std::cout << "IS results:" << std::endl;
    // Do each trial
    for (const auto item : is_b_trials) {
        set_mu_sigma_pair(is_bench_10k, is_bench_100k, item.first, item.second);
        const time_t trial_10k_start {std::clock()};
        std::unique_ptr<MC::Result> trail_10k_result {MC::Is_Benchmark_Is(is_bench_10k)};
        const time_t trial_10k_end {std::clock()};
        const time_t trial_100k_start {std::clock()};
        std::unique_ptr<MC::Result> trail_100k_result {MC::Is_Benchmark_Is(is_bench_100k)};
        const time_t trial_100k_end {std::clock()};
        std::cout << std::setprecision(1) << " 10k samples mu: " << item.first << " sigma: " << item.second << ": ";
        print_result(trail_10k_result->value, trail_10k_result->error, trial_10k_start, trial_10k_end);
        std::cout << std::setprecision(1) << "100k samples mu: " << item.first << " sigma: " << item.second << ": ";
        print_result(trail_100k_result->value, trail_100k_result->error, trial_100k_start, trial_100k_end);
    }



    return 0;
}