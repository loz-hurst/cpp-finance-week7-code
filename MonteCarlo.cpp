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
#include <utility>
#include <vector>
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
        result->error = discount*std::sqrt(accumulator_squares/data.paths - result->value*result->value)/data.paths;

        return result;

    }

    std::unique_ptr<Result> Ln_S(const Data& data) {
        // This follows the same process as Plain - so see the comments in there for the most part
        std::unique_ptr<Result> result {std::make_unique<Result>(Result{0, 0})};

        const double delta_t {data.maturity/data.steps};
        const double drift {(data.rate-0.5*data.sigma*data.sigma)*delta_t};
        const double sigma_sqrt_delta_t {data.sigma*std::sqrt(delta_t)};
        const double discount {std::exp(-data.rate * data.maturity)};

        // Find log of S_0 - we will need that for each path
        const double ln_S_0 {std::log(data.S_0)};

        double accumulator_values {0};
        double accumulator_squares {0};

        for(long i {0}; data.paths > i; ++i){
            Utility::print_progress(i+1, 50000);

            // Starting from ln(S_0) instead of S_0
            double ln_S {ln_S_0};

            for(long j {0}; data.steps > j; ++j) {
                double w {Random::GetNormalValue()};
                // Calculating ln_S now, so addition instead of product and no exp to find
                ln_S += drift + sigma_sqrt_delta_t*w;
            }

            const double S {std::exp(ln_S)}; // This is the only call to exp for each path

            const double payoff {std::max(S-data.strike, 0.0)};

            accumulator_values += payoff;
            accumulator_squares += payoff*payoff;
        }
        Utility::print_clear();

        result->value = discount*accumulator_values/data.paths;
        result->error = discount*std::sqrt(accumulator_squares/data.paths - result->value*result->value)/data.paths;

        return result;

    }

    std::unique_ptr<Result> Cv(const Data& data, double & correlation) {
        // This follows the same process as Plain - so see the comments in there for the most part
        std::unique_ptr<Result> result {std::make_unique<Result>(Result{0, 0})};

        const double delta_t {data.maturity/data.steps};
        const double drift {(data.rate-0.5*data.sigma*data.sigma)*delta_t};
        const double sigma_sqrt_delta_t {data.sigma*std::sqrt(delta_t)};
        const double discount {std::exp(-data.rate * data.maturity)};

        // expected future stock price (e^(rT)S_0)
        const double discount_S_0 {std::exp(data.rate*data.maturity) * data.S_0};

        // This time we need to store each of the values as well as some accumulators

        // Since we know the size (number of paths), using an array is going to be more memory efficient than a vector
        std::vector<double> payoffs (data.paths); // Payoffs
        std::vector<double> cvs (data.paths); // Control variants

        // This time we will also need to know the sum of the control variants
        double accumulator_cv_values {0};
        // ...and the sum of the control variants squared
        double accumulator_cv_squares {0};
        // ...and the sum of the product of payoff and control variants
        double accumulator_p_cv {0};
        // ... as well as the sum of the payoffs and the payoffs squared (as before, but new variable names)
        double accumulator_p_values {0};
        double accumulator_p_squares {0};

        for(long i {0}; data.paths > i; ++i) {
            Utility::print_progress(i + 1, 50000);

            double S{data.S_0};

            for (long j{0}; data.steps > j; ++j) {
                double w{Random::GetNormalValue()};
                S = S * std::exp(drift + sigma_sqrt_delta_t * w);
            }

            // Store the payoff
            payoffs[i] = std::max(S - data.strike, 0.0);
            cvs[i] = S - discount_S_0;

            accumulator_cv_values += cvs[i];
            accumulator_cv_squares += cvs[i] * cvs[i];
            accumulator_p_cv += payoffs[i] * cvs[i];
            accumulator_p_values += payoffs[i];
            accumulator_p_squares += payoffs[i]*payoffs[i];
        }
        Utility::print_clear();

        // Crudely calculate beta = cov(c,d)/var(d)
        const double beta {(accumulator_p_cv-accumulator_p_values*accumulator_cv_values/data.paths)/
            (accumulator_cv_squares-accumulator_cv_values*accumulator_cv_values/data.paths)};

        correlation = beta*std::sqrt(
            (data.paths*accumulator_cv_squares - accumulator_cv_values*accumulator_cv_values)/
            (data.paths*accumulator_p_squares - accumulator_p_values*accumulator_p_values)
        );

        // Now calculate the corrected payoff using cp = payoff_i - beta*cv_i
        // We will need the accumulated sum and squares to calculate the final value and error, respectively
        double accumulator_values {0};
        double accumulator_squares {0};
        for (int i {0}; data.paths > i; ++i) {
            const double correct_p {payoffs[i] - beta*cvs[i]}; // corrected payoff
            accumulator_values += correct_p;
            accumulator_squares += correct_p*correct_p;
        }

        result->value = discount*accumulator_values/data.paths;
        result->error = discount*std::sqrt(accumulator_squares/data.paths - result->value*result->value)/data.paths;

        return result;

    }

    std::unique_ptr<Result> Is_Benchmark_Plain(const Data& data) {
        std::unique_ptr<Result> result {std::make_unique<Result>(Result{0, 0})};

        const double quartile {0.99}; // Nintey-ninth quartile
        // Ideally would calculate from quartile but non-trival and easily looked up.
        const double inv_quartile {2.3263};

        const double quartile_over {1-quartile}; // proportion over the quartile

        // Store the square of the "accuracy" (the difference^2 of samples were over the quartile vs proportion expected)
        double accumulator_acc_square {0};
        // Store the values so we can return an average
        double accumulator_value {0};

        // Slightly fudging what we use steps and paths for - will do "steps" repeats of "paths" samples
        for (int i {0}; data.steps > i; ++i) {
            double accumulator_p {0};

            for (int j {0}; data.paths > j; ++j) {
                const double x {Random::GetNormalValue()};
                double p {(x >= inv_quartile) ? 1.0 : 0.0};

                accumulator_p += p;
            }

            double value {accumulator_p/data.paths}; // should be our proportion over the quartile
            accumulator_acc_square += (value-quartile_over)*(value-quartile_over);
            accumulator_value += value;
        }


        result->value = accumulator_value/data.steps;
        result->error = std::sqrt(accumulator_acc_square)/data.steps;

        return result;
    }

    std::unique_ptr<Result> Is_Benchmark_Is(const Data& data) {
        // Mostly the same as Is_Benchmark_Plain so see those comments where required
        std::unique_ptr<Result> result {std::make_unique<Result>(Result{0, 0})};

        const double quartile {0.99};
        const double inv_quartile {2.3263};

        const double quartile_over {1-quartile};

        double accumulator_acc_square {0};
        double accumulator_value {0};

        for (int i {0}; data.steps > i; ++i) {
            double accumulator_p {0};

            for (int j {0}; data.paths > j; ++j) {
                // Random value from our non-standard normal distribution
                const double x {data.sigma*Random::GetNormalValue()+data.rate};
                double p {0};
                if (x >= inv_quartile) {
                    // This needs to be weighted to fit our new distribution back into the old one
                    double x_minus_mu_div_sig {(x-data.rate)/data.sigma};
                    p = data.sigma*std::exp(-0.5*(x*x - x_minus_mu_div_sig*x_minus_mu_div_sig));
                }

                accumulator_p += p;
            }

            double value {accumulator_p/data.paths};
            accumulator_acc_square += (value-quartile_over)*(value-quartile_over);
            accumulator_value += value;
        }


        result->value = accumulator_value/data.steps;
        result->error = std::sqrt(accumulator_acc_square)/data.steps;

        return result;

    }

}