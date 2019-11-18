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

#include <algorithm>
#include <cmath>
#include <memory>
#include "MonteCarlo.hpp"
#include "Random.hpp"
#include "Utility.hpp"

namespace MonteCarlo {
    std::unique_ptr<Result> Plain(const Data& data) {
        std::unique_ptr<Result> result {std::make_unique<Result>(Result{0, 0})};

        /* Firstly we want to pre-calculate everything we can outside of any loop.
         * (Only calculating values once, where possible, is clearly more efficient.)
         */

        // Size of each step
        const double delta_t {data.maturity/data.steps};

        /* Use (mu-0.5sigma^2)delta_t as the drift - remembering mu is set to
         * rate under risk neutrality.
         */
        const double drift {(data.rate-0.5*data.sigma*data.sigma)*delta_t};

        /* We need sigma*sqrt(delta_t) at each step - but it's a constant value
         * for this simulation.
         */
        const double sigma_sqrt_delta_t {data.sigma*std::sqrt(delta_t)};

        // Discount: e^(-rT)
        const double discount {std::exp(-data.rate * data.maturity)};

        // Accumulate the values found (i.e. max(S_t-X, 0))
        double accumulator_values {0};
        /* Accumulate the squares so we can calculate the error (using standard
         *  deviation) via sqrt(sum(c_t^2)/M-c^2)
         */
        double accumulator_squares {0};

        // Simulate each of the paths
        for(long i {0}; data.paths > i; ++i){

            // Print the progress (so there's some output)
            Utility::print_progress(i+1, 50000);
            double S {data.S_0}; // Start each path from S_0
            // For each path, simulate each of the steps
            for(long j {0}; data.steps > j; ++j) {
                // Get our Brownian value
                double w {Random::GetNormalValue()};
                /* S_t+delta_t = S_t * exp(drift + sigma_sqrt_delta_t*w)
                 * Where drift = (mu-0.5sigma^2)delta_t
                 *       sigma_sqrt_delta_t = sigma * sqrt(delta_t)
                 */
                S = S * std::exp(drift + sigma_sqrt_delta_t*w);
            }
            // Done all the time steps.

            // Work out the payoff
            const double payoff {std::max(S-data.strike, 0.0)};

            // Add it to our running sums
            accumulator_values += payoff;
            accumulator_squares += payoff*payoff;
        }
        Utility::print_clear();

        /* est(c) = 1/paths * discount * accumulator_values
         * Where: paths = number of simulations
         *        discount = e^(-rT)
         *        accumulator_values = sum(max(S_t - X), 0)
         */
        result->value = discount*accumulator_values/data.paths;
        result->error = std::sqrt(accumulator_squares/data.paths - result->value*result->value);

        return result;

    }
}