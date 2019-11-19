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

#ifndef CPP_FINANCE_WEEK7_CODE_MONTECARLO_HPP
#define CPP_FINANCE_WEEK7_CODE_MONTECARLO_HPP

#include <memory>

namespace MonteCarlo {
    struct Data {
        double S_0;
        double rate;
        double sigma;
        double maturity;
        double strike;

        long paths;
        long steps;
    };

    struct Result {
        double value;
        double error;
    };

    /* Calculates an estimate for the value using the exact solution to the GBM
     * equation (dS_t = (mu*S_t)dt + (sigma*S_t)dz_t):
     * S_t+delta_t = S_t*exp( (mu-0.5sigma^2)delta_t + sigma*sqrt(delta_t)*w_i
     * using plain Monte Carlo.
     * Arguments:
     *   Data containing required simulation parameters
     * Results:
     *   Result contain the estimate and a measure of its error
     */
    std::unique_ptr<Result> Plain(const Data &);

    // Exactly the same as Plain but finds R_t=ln(S_t) internally before recovering S_t with a single exponent call.
    std::unique_ptr<Result> Ln_S(const Data &);

    /* Plain Monte Carlo modified to use a control variate (CV) d = S_t - e^(-rT)*S_0
     * Puts the calculated correlation in the double passed by reference as the 2nd argument.
     */
    std::unique_ptr<Result> Cv(const Data &, double & correlation);

    // Importance Sampling benchmark using plain Monte Carlo
    std::unique_ptr<Result> Is_Benchmark_Plain(const Data &);

    // Importance Sampling benchmark using Importance Sampling
    std::unique_ptr<Result> Is_Benchmark_Is(const Data &);
}

#endif //CPP_FINANCE_WEEK7_CODE_MONTECARLO_HPP
