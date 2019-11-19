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

#include <chrono>
#include <random>
#include "Random.hpp"


// Anonymous namespace - initialise the random number generator as early as possible
namespace {
    std::mt19937 engine ((unsigned int)std::chrono::system_clock::now().time_since_epoch().count());
    //std::default_random_engine engine;
    std::normal_distribution<double> standard_normal_distribution(0, 1);
    std::uniform_real_distribution<double> uniform_0_1(0, 1);
}

namespace Random {
    double GetRandom() {
        return ::uniform_0_1(::engine);
    }

    double GetNormalValue() {
        return ::standard_normal_distribution(::engine);
    }
}