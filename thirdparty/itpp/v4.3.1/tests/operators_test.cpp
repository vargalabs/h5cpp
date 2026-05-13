/*!
 * \file
 * \brief Operators test program
 * \author Bogdan Cristea
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2011  (see AUTHORS file for a list of contributors)
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

#include "itpp/itbase.h"

using namespace itpp;
using std::cout;
using std::endl;
using std::complex;

int main()
{
	//complex plus/minus scalar
	{
		//complex scalar
		complex<double> x(1.0, 1.0);
		int y = 2;
		complex<double> z = x+y;
		if (complex<double>(x.real()+y, x.imag()) != z)
		{
			cout << "complex scalar plus int scalar operator is wrong" << endl;
			return EXIT_FAILURE;
		}
		z = x-y;
		if (complex<double>(x.real()-y, x.imag()) != z)
		{
			cout << "complex scalar minus int scalar operator is wrong" << endl;
			return EXIT_FAILURE;
		}

		//complex vector
		vec x1_re("1 3");
		vec x1_im("2 4");
		cvec x1 = to_cvec(x1_re, x1_im);
		cvec z1 = x1+y;
		if (to_cvec(x1_re+y, x1_im) != z1)
		{
			cout << "complex vector plus int scalar operator is wrong" << endl;
			return EXIT_FAILURE;
		}
		z1 = x1-y;
		if (to_cvec(x1_re-y, x1_im) != z1)
		{
			cout << "complex vector minus int scalar operator is wrong" << endl;
			return EXIT_FAILURE;
		}

		//complex scalar
		float y1 = 2.7;
		z = x+y1;
		if (complex<double>(x.real()+y1, x.imag()) != z)
		{
			cout << "complex plus float scalar operator is wrong" << endl;
			return EXIT_FAILURE;
		}
		z = x-y1;
		if (complex<double>(x.real()-y1, x.imag()) != z)
		{
			cout << "complex minus float scalar operator is wrong" << endl;
			return EXIT_FAILURE;
		}

		//complex vector
		z1 = x1+y1;
		if (to_cvec(x1_re+y1, x1_im) != z1)
		{
			cout << "complex vector plus float scalar operator is wrong" << endl;
			return EXIT_FAILURE;
		}
		z1 = x1-y1;
		if (to_cvec(x1_re-y1, x1_im) != z1)
		{
			cout << "complex vector minus float scalar operator is wrong" << endl;
			return EXIT_FAILURE;
		}

		//complex matrix
		mat x2_re("1 3; 5 7");
		mat x2_im("2 4; 6 8");
		cmat x2 = to_cmat(x2_re, x2_im);
		cmat z2 = x2+complex<double>(y1, 0);
		if (to_cmat(x2_re+double(y1), x2_im) != z2)
		{
			cout << "complex matrix plus complex scalar operator is wrong" << endl;
			return EXIT_FAILURE;
		}
		z2 = x2-complex<double>(y1, 0);
		if (to_cmat(x2_re-double(y1), x2_im) != z2)
		{
			cout << "complex matrix minus complex scalar operator is wrong" << endl;
			return EXIT_FAILURE;
		}
	}
	//scalar plus/minus complex
	{
		//complex scalar
		int x = 2;
		complex<double> y(1.0, 1.0);
		complex<double> z = x+y;
		if (complex<double>(x+y.real(), y.imag()) != z)
		{
			cout << "int scalar plus complex scalar operator is wrong" << endl;
			return EXIT_FAILURE;
		}
		z = x-y;
		if (complex<double>(x-y.real(), -y.imag()) != z)
		{
			cout << "int scalar minus complex scalar operator is wrong" << endl;
			return EXIT_FAILURE;
		}

		//complex vector
		vec y1_re("1 3");
		vec y1_im("2 4");
		cvec y1 = to_cvec(y1_re, y1_im);
		cvec z1 = x+y1;
		if (to_cvec(x+y1_re, y1_im) != z1)
		{
			cout << "int scalar plus complex vector operator is wrong" << endl;
			return EXIT_FAILURE;
		}
		z1 = x-y1;
		if (to_cvec(x-y1_re, -y1_im) != z1)
		{
			cout << "int scalar minus complex vector operator is wrong" << endl;
			return EXIT_FAILURE;
		}

		//complex scalar
		float x1 = 2.7;
		z = x1+y;
		if (complex<double>(x1+y.real(), y.imag()) != z)
		{
			cout << "float scalar plus complex operator is wrong" << endl;
			return EXIT_FAILURE;
		}
		z = x1-y;
		if (complex<double>(x1-y.real(), -y.imag()) != z)
		{
			cout << "float scalar minus complex operator is wrong" << endl;
			return EXIT_FAILURE;
		}

		//complex vector
		z1 = x1+y1;
		if (to_cvec(x1+y1_re, y1_im) != z1)
		{
			cout << "float scalar plus complex vector operator is wrong" << endl;
			return EXIT_FAILURE;
		}
		z1 = x1-y1;
		if (to_cvec(x1-y1_re, -y1_im) != z1)
		{
			cout << "float scalar minus complex vector operator is wrong" << endl;
			return EXIT_FAILURE;
		}

		//complex matrix
		mat y2_re("1 3; 5 7");
		mat y2_im("2 4; 6 8");
		cmat y2 = to_cmat(y2_re, y2_im);
		cmat z2 = complex<double>(x1, 0)+y2;
		if (to_cmat(x1+y2_re, y2_im) != z2)
		{
			cout << "complex scalar plus complex matrix operator is wrong" << endl;
			return EXIT_FAILURE;
		}
		z2 = complex<double>(x1, 0)-y2;
		if (to_cmat(x1-y2_re, -y2_im) != z2)
		{
			cout << "complex scalar minus complex matrix operator is wrong" << endl;
			return EXIT_FAILURE;
		}
	}
	cout << "operators test successful" << endl;
	return EXIT_SUCCESS;
}
