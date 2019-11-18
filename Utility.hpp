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

#ifndef CPP_FINANCE_WEEK7_CODE_UTILITY_HPP
#define CPP_FINANCE_WEEK7_CODE_UTILITY_HPP

namespace Utility {
    /* Prints a progress message on a single line (does not clear the line)
     * Arguments:
     *     done - how many units are complete
     *     update_freq - only print after each (this many) units
     */
    void print_progress(long done, long update_freq=1);

    // Convenience function to move to the next output line after using print_progress
    void print_clear();
}

#endif //CPP_FINANCE_WEEK7_CODE_UTILITY_HPP
