/*!
 * \file
 * \brief linspace functions test program
 * \author Bogdan Cristea
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include <itpp/itbase.h>

using namespace std;
using namespace itpp;

#define TOL 1e-12
int main(void)
{
    cout << "================================" << endl;
    cout << "    Test of linspace() function " << endl;
    cout << "================================" << endl;

    double from = 0;
    double to = 2;
    int points = 1;
    vec v = linspace(from, to, points);
    cout << "Generating a vector of " << points << " points, from " <<
    from << " to " << to << std::endl;
    if (points != v.length())
    {
        cout << "Wrong vector length: " << v.length() << endl;
        cout << "Vector length should be " << points << endl;
    } else
    {
        cout << "Vector length is " << v.length() << endl;
    }
    if (to != v(0))
    {
        cout << "Wrong vector element: " << v << endl;
        cout << "Vector element should be " << to << endl;
    } else
    {
        cout << "Vector element is " << v(0) << endl;
    }

    from = 0;
    to = 10;
    points = 30;
    v = linspace(from, to, points);
    cout << "Generating a vector of " << points << " points, from " <<
    from << " to " << to << endl;
    if (points != v.length())
    {
        cout << "Wrong vector length: " << v.length() << endl;
        cout << "Vector length should be " << points << endl;
    } else
    {
        cout << "Vector length is " << v.length() << endl;
    }
    vec actual_v(v.length());
    int n;
    for (n = 0; n < v.length()-1; ++n)
    {
        actual_v(n) = n*(to-from)/(points-1);
    }
    actual_v(n) = to;
    if (any(abs(actual_v - v) > TOL) || from != v(0) || to != v(v.length()-1))
    {
        cout << "Wrong vector elements: " << v << endl;
        cout << "Vector elements should be: " << actual_v << endl;
    } else
    {
        cout << "Vector elements are: " << v << endl;
    }

    cout << "==========================================" << endl;
    cout << "    Test of linspace_fixed_step() function " << endl;
    cout << "==========================================" << endl;

    from = 0;
    to = 10;
    double step = 1.0;
    v = linspace_fixed_step(from, to, step);
    cout << "Generating a vector with a step of " << step <<
    ", from " << from << " to " << to << endl;
    int actual_len = floor_i((to-from)/step)+1;
    if (actual_len != v.length())
    {
        cout << "Wrong vector length: " << v.length() << endl;
        cout << "Vector length should be: " << actual_len << endl;
    } else
    {
        cout << "Vector length is: " << v.length() << endl;
    }
    actual_v.set_size(actual_len);
    for (n = 0; n < actual_v.length(); ++n)
    {
        actual_v(n) = from+n*step;
    }
    if (any(abs(actual_v - v) > TOL) || from != v(0) || to < v(v.length()-1))
    {
        cout << "Wrong vector elements: " << v << endl;
        cout << "Vector elements should be: " << actual_v << endl;
    } else
    {
        cout << "Vector elements are: " << v << endl;
    }

    from = 0;
    to = 10;
    step = 1.1;
    v = linspace_fixed_step(from, to, step);
    cout << "Generating a vector with a step of " << step <<
    ", from " << from << " to " << to << endl;
    actual_len = floor_i((to-from)/step)+1;
    if (actual_len != v.length())
    {
        cout << "Wrong vector length: " << v.length() << endl;
        cout << "Vector length should be: " << actual_len << endl;
    } else
    {
        cout << "Vector length is: " << v.length() << endl;
    }
    actual_v.set_size(actual_len);
    for (n = 0; n < actual_v.length(); ++n)
    {
        actual_v(n) = from+n*step;
    }
    if (any(abs(actual_v - v) > TOL) || from != v(0) || to < v(v.length()-1))
    {
        cout << "Wrong vector elements: " << v << endl;
        cout << "Vector elements should be: " << actual_v << endl;
    } else
    {
        cout << "Vector elements are: " << v << endl;
    }

    from = 0;
    to = 1;
    step = 2;
    v = linspace_fixed_step(from, to, step);
    cout << "Generating a vector with a step of " << step <<
    ", from " << from << " to " << to << endl;
    actual_len = 1;
    if (actual_len != v.length())
    {
        cout << "Wrong vector length: " << v.length() << endl;
        cout << "Vector length should be: " << actual_len << endl;
    } else
    {
        cout << "Vector length is: " << v.length() << endl;
    }
    actual_v.set_size(actual_len);
    for (n = 0; n < actual_v.length(); ++n)
    {
        actual_v(n) = from+n*step;
    }
    if (from != v(0))
    {
        cout << "Wrong vector element: " << v << endl;
        cout << "Vector element should be: " << actual_v << endl;
    } else
    {
        cout << "Vector element is: " << v << endl;
    }

    from = 0;
    to = -1;
    step = -2;
    v = linspace_fixed_step(from, to, step);
    cout << "Generating a vector with a step of " << step <<
    ", from " << from << " to " << to << endl;
    actual_len = 1;
    if (actual_len != v.length())
    {
        cout << "Wrong vector length: " << v.length() << endl;
        cout << "Vector length should be: " << actual_len << endl;
    } else
    {
        cout << "Vector length is: " << v.length() << endl;
    }
    actual_v.set_size(actual_len);
    for (n = 0; n < actual_v.length(); ++n)
    {
        actual_v(n) = from+n*step;
    }
    if (from != v(0))
    {
        cout << "Wrong vector element: " << v << endl;
        cout << "Vector element should be: " << actual_v << endl;
    } else
    {
        cout << "Vector element is: " << v << endl;
    }

    from = 0;
    to = -1;
    step = 2;
    v = linspace_fixed_step(from, to, step);
    cout << "Generating a vector with a step of " << step <<
    ", from " << from << " to " << to << endl;
    actual_len = 0;
    if (actual_len != v.length())
    {
        cout << "Wrong vector length: " << v.length() << endl;
        cout << "Vector length should be: " << actual_len << endl;
    } else
    {
        cout << "Vector length is: " << v.length() << endl;
    }
    actual_v.set_size(actual_len);
    if (v != actual_v)
    {
        cout << "Wrong vector element: " << v << endl;
        cout << "Vector element should be: " << actual_v << endl;
    } else
    {
        cout << "Vector element is: " << v << endl;
    }

    from = -1;
    to = 0;
    step = 2;
    v = linspace_fixed_step(from, to, step);
    cout << "Generating a vector with a step of " << step <<
    ", from " << from << " to " << to << endl;
    actual_len = 1;
    if (actual_len != v.length())
    {
        cout << "Wrong vector length: " << v.length() << endl;
        cout << "Vector length should be: " << actual_len << endl;
    } else
    {
        cout << "Vector length is: " << v.length() << endl;
    }
    actual_v.set_size(actual_len);
    actual_v(0) = from;
    if (v != actual_v)
    {
        cout << "Wrong vector element: " << v << endl;
        cout << "Vector element should be: " << actual_v << endl;
    } else
    {
        cout << "Vector element is: " << v << endl;
    }

    from = -1;
    to = 0;
    step = -2;
    v = linspace_fixed_step(from, to, step);
    cout << "Generating a vector with a step of " << step <<
    ", from " << from << " to " << to << endl;
    actual_len = 0;
    if (actual_len != v.length())
    {
        cout << "Wrong vector length: " << v.length() << endl;
        cout << "Vector length should be: " << actual_len << endl;
    } else
    {
        cout << "Vector length is: " << v.length() << endl;
    }
    actual_v.set_size(actual_len);
    if (v != actual_v)
    {
        cout << "Wrong vector element: " << v << endl;
        cout << "Vector element should be: " << actual_v << endl;
    } else
    {
        cout << "Vector element is: " << v << endl;
    }

    from = -5;
    to = 5;
    step = 2;
    v = linspace_fixed_step(from, to, step);
    cout << "Generating a vector with a step of " << step <<
    ", from " << from << " to " << to << endl;
    actual_len = floor_i((to-from)/step)+1;
    if (actual_len != v.length())
    {
        cout << "Wrong vector length: " << v.length() << endl;
        cout << "Vector length should be: " << actual_len << endl;
    } else
    {
        cout << "Vector length is: " << v.length() << endl;
    }
    actual_v.set_size(actual_len);
    for (n = 0; n < actual_v.length(); ++n)
    {
        actual_v(n) = from+n*step;
    }
    if (v != actual_v)
    {
        cout << "Wrong vector element: " << v << endl;
        cout << "Vector element should be: " << actual_v << endl;
    } else
    {
        cout << "Vector element is: " << v << endl;
    }

    from = 5;
    to = -5;
    step = -2;
    v = linspace_fixed_step(from, to, step);
    cout << "Generating a vector with a step of " << step <<
    ", from " << from << " to " << to << endl;
    actual_len = floor_i((to-from)/step)+1;
    if (actual_len != v.length())
    {
        cout << "Wrong vector length: " << v.length() << endl;
        cout << "Vector length should be: " << actual_len << endl;
    } else
    {
        cout << "Vector length is: " << v.length() << endl;
    }
    actual_v.set_size(actual_len);
    for (n = 0; n < actual_len; ++n)
    {
        actual_v(n) = from+n*step;
    }
    if (v != actual_v)
    {
        cout << "Wrong vector element: " << v << endl;
        cout << "Vector element should be: " << actual_v << endl;
    } else
    {
        cout << "Vector element is: " << v << endl;
    }

    int ifrom = 0;
    int ito = 10;
    ivec iv = linspace_fixed_step(ifrom, ito);
    cout << "Generating an integer vector with a step of 1" <<
    ", from " << ifrom << " to " << ito << endl;
    actual_len = floor_i(ito-ifrom)+1;
    if (actual_len != iv.length())
    {
        cout << "Wrong vector length: " << iv.length() << endl;
        cout << "Vector length should be: " << actual_len << endl;
    } else
    {
        cout << "Vector length is: " << iv.length() << endl;
    }
    ivec actual_iv(actual_len);
    for (n = 0; n < actual_iv.length(); ++n)
    {
        actual_iv(n) = ifrom+n;
    }
    if (iv != actual_iv)
    {
        cout << "Wrong vector element: " << iv << endl;
        cout << "Vector element should be: " << actual_iv << endl;
    } else
    {
        cout << "Vector element is: " << iv << endl;
    }

    short int sfrom = 0;
    short int sto = 10;
    svec sv = linspace_fixed_step(sfrom, sto);
    cout << "Generating a short integer vector with a step of 1" <<
    ", from " << sfrom << " to " << sto << endl;
    actual_len = floor_i(sto-sfrom)+1;
    if (actual_len != sv.length())
    {
        cout << "Wrong vector length: " << sv.length() << endl;
        cout << "Vector length should be: " << actual_len << endl;
    } else
    {
        cout << "Vector length is: " << sv.length() << endl;
    }
    svec actual_sv(actual_len);
    for (n = 0; n < actual_sv.length(); ++n)
    {
        actual_sv(n) = sfrom+n;
    }
    if (sv != actual_sv)
    {
        cout << "Wrong vector element: " << sv << endl;
        cout << "Vector element should be: " << actual_sv << endl;
    } else
    {
        cout << "Vector element is: " << sv << endl;
    }
}
