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
#include <iostream>
#include "MonteCarlo.hpp"

namespace MC = MonteCarlo;

// Print out the MC simulation result - including value, error and calculation of time taken from (end-start)/(ticks/s)
void print_result(const double value, const double error, const time_t start, const time_t end) {
    std::cout << "Calculated estimate value of " << value
              << " with error " << error
              << " in " << ((end-start)/CLOCKS_PER_SEC) << "s"
              << std::endl;
}

int main() {
    // S_0, rate, sigma, maturity, strike, paths, steps
    MC::Data option {100, 0.05, 0.20, 1, 100, 1000000, 100};

    std::cout << "Begin plain Monte Carlo:" << std::endl;
    time_t plain_mc_start {std::clock()};
    const std::unique_ptr<MC::Result> plain_result {MC::Plain(option)};
    time_t plain_mc_end {std::clock()};
    print_result(plain_result->value, plain_result->error, plain_mc_start, plain_mc_end);

    std::cout << "Begin ln_S Monte Carlo:" << std::endl;
    time_t ln_S_mc_start {std::clock()};
    const std::unique_ptr<MC::Result> ln_S_result {MC::Ln_S(option)};
    time_t ln_S_mc_end {std::clock()};
    print_result(ln_S_result->value, ln_S_result->error, ln_S_mc_start, ln_S_mc_end);

    std::cout << "Begin CV Monte Carlo:" << std::endl;
    time_t cv_mc_start {std::clock()};
    double correlation;
    const std::unique_ptr<MC::Result> cv_result {MC::Cv(option, correlation)};
    time_t cv_mc_end {std::clock()};
    print_result(cv_result->value, cv_result->error, cv_mc_start, cv_mc_end);
    std::cout << "Correlation was: " << correlation << std::endl;



    return 0;
}