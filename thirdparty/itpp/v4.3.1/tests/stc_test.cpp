/*!
 * \file
 * \brief STC class test program
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

#include "itpp/itcomm.h"

//returns true if the absolute difference of all elements of matrices is below err
bool near(const itpp::cmat &M1, const itpp::cmat &M2, double err)
{
	if (M1.rows() != M2.rows() || M1.cols() != M2.cols())
	{
		return false;
	}

	register int n, k;
	for (n = 0; n < M1.rows(); ++n)
	{
		for (k = 0; k < M1.cols(); ++k)
		{
			if (std::abs(M1(n,k)-M2(n,k)) >= err)
			{
				return false;
			}
		}
	}
	return true;
}

int main(void)
{
	int em_antenna = 2;
	int channel_uses = 2;
	const int const_size = 4;//QAM
	const double abs_err = 1e-12;

	std::cout << "Space Time Code class test started" << std::endl;

	std::cout << "Checking V-BLAST" << std::endl;
	std::string code_name("V-BLAST_MxN");

	//should succeed
	itpp::STC stc(code_name, const_size, em_antenna, channel_uses);

	//check parameters
	if (2 != stc.get_nb_emission_antenna())
	{
		std::cout << "Number of emission antenna is wrong" << std::endl;
	}
	if (2 != stc.get_channel_uses())
	{
		std::cout << "Channel uses is wrong" << std::endl;
	}
	if (4 != stc.get_nb_symbols_per_block())
	{
		std::cout << "Number of symbols per block is wrong" << std::endl;
	}

	//check generated matrices
	itpp::cmat ref_vblast("1 0; 0 0; 0 1; 0 0; 0 0; 1 0; 0 0; 0 1");
	if (!near(ref_vblast, stc.get_1st_gen_matrix(), abs_err) ||
		!near(ref_vblast, stc.get_2nd_gen_matrix(), abs_err))
	{
		std::cout << "Generator matrices are different" << std::endl;
		std::cout << "ref_vblast = " << ref_vblast << std::endl;
		std::cout << "A = " << stc.get_1st_gen_matrix() << std::endl;
	}

	std::cout << "Checking improved V-BLAST" << std::endl;
	code_name = "imp_V-BLAST_MxN";

	//should succeed
	stc.setup(code_name, const_size, em_antenna, channel_uses);

	//check parameters
	if (2 != stc.get_nb_emission_antenna())
	{
		std::cout << "Number of emission antenna is wrong" << std::endl;
	}
	if (2 != stc.get_channel_uses())
	{
		std::cout << "Channel uses is wrong" << std::endl;
	}
	if (4 != stc.get_nb_symbols_per_block())
	{
		std::cout << "Number of symbols per block is wrong" << std::endl;
	}

	//check generated matrices
	itpp::cmat ref_imp_vblast("1 0; 0 1; 0 1; 1 0; 1 0; 0 -1; 0 1; -1 0");
	ref_imp_vblast /= std::sqrt(2.0);
	if (!near(ref_imp_vblast, stc.get_1st_gen_matrix(), abs_err) ||
		!near(ref_imp_vblast, stc.get_2nd_gen_matrix(), abs_err))
	{
		std::cout << "Generator matrices are different" << std::endl;
		std::cout << "ref_imp_vblast = " << ref_imp_vblast << std::endl;
		std::cout << "A = " << stc.get_1st_gen_matrix() << std::endl;
	}

	std::cout << "Checking Alamouti code" << std::endl;
	code_name = "Alamouti_2xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	if (2 != stc.get_nb_emission_antenna())
	{
		std::cout << "Number of emission antenna is wrong" << std::endl;
	}
	if (2 != stc.get_channel_uses())
	{
		std::cout << "Channel uses is wrong" << std::endl;
	}
	if (2 != stc.get_nb_symbols_per_block())
	{
		std::cout << "Number of symbols per block is wrong" << std::endl;
	}

	//check generated matrices
	itpp::cmat ref_A("1 0; 0 1; 0 1; -1 0");
	itpp::cmat ref_B("1 0; 0 -1; 0 1; 1 0");
	if (!near(ref_A, stc.get_1st_gen_matrix(), abs_err) ||
		!near(ref_B, stc.get_2nd_gen_matrix(), abs_err))
	{
		std::cout << "Generator matrices are different" << std::endl;
		std::cout << "ref_A = " << ref_A << std::endl;
		std::cout << "A = " << stc.get_1st_gen_matrix() << std::endl;
		std::cout << "ref_B = " << ref_B << std::endl;
		std::cout << "B = " << stc.get_2nd_gen_matrix() << std::endl;
	}

	std::cout << "Checking switched Alamouti code" << std::endl;
	code_name = "Switched_Alamouti_4xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	if (4 != stc.get_nb_emission_antenna())
	{
		std::cout << "Number of emission antenna is wrong" << std::endl;
	}
	if (4 != stc.get_channel_uses())
	{
		std::cout << "Channel uses is wrong" << std::endl;
	}
	if (4 != stc.get_nb_symbols_per_block())
	{
		std::cout << "Number of symbols per block is wrong" << std::endl;
	}

	std::cout << "Checking double Alamouti code" << std::endl;
	code_name = "Double_Alamouti_4xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	if (4 != stc.get_nb_emission_antenna())
	{
		std::cout << "Number of emission antenna is wrong" << std::endl;
	}
	if (2 != stc.get_channel_uses())
	{
		std::cout << "Channel uses is wrong" << std::endl;
	}
	if (4 != stc.get_nb_symbols_per_block())
	{
		std::cout << "Number of symbols per block is wrong" << std::endl;
	}

	std::cout << "Checking Jafarkhani code" << std::endl;
	code_name = "Jafarkhani_4xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	if (4 != stc.get_nb_emission_antenna())
	{
		std::cout << "Number of emission antenna is wrong" << std::endl;
	}
	if (4 != stc.get_channel_uses())
	{
		std::cout << "Channel uses is wrong" << std::endl;
	}
	if (4 != stc.get_nb_symbols_per_block())
	{
		std::cout << "Number of symbols per block is wrong" << std::endl;
	}

	std::cout << "Checking Golden code" << std::endl;
	code_name = "Golden_2x2";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	if (2 != stc.get_nb_emission_antenna())
	{
		std::cout << "Number of emission antenna is wrong" << std::endl;
	}
	if (2 != stc.get_channel_uses())
	{
		std::cout << "Channel uses is wrong" << std::endl;
	}
	if (4 != stc.get_nb_symbols_per_block())
	{
		std::cout << "Number of symbols per block is wrong" << std::endl;
	}

	std::cout << "Checking Damen code" << std::endl;
	code_name = "Damen_2x2";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	if (2 != stc.get_nb_emission_antenna())
	{
		std::cout << "Number of emission antenna is wrong" << std::endl;
	}
	if (2 != stc.get_channel_uses())
	{
		std::cout << "Channel uses is wrong" << std::endl;
	}
	if (4 != stc.get_nb_symbols_per_block())
	{
		std::cout << "Number of symbols per block is wrong" << std::endl;
	}

	std::cout << "Checking orthogonal code 3xN" << std::endl;
	code_name = "34ortho_3xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	if (3 != stc.get_nb_emission_antenna())
	{
		std::cout << "Number of emission antenna is wrong" << std::endl;
	}
	if (4 != stc.get_channel_uses())
	{
		std::cout << "Channel uses is wrong" << std::endl;
	}
	if (3 != stc.get_nb_symbols_per_block())
	{
		std::cout << "Number of symbols per block is wrong" << std::endl;
	}

	std::cout << "Checking 36 LD code 3xN" << std::endl;
	code_name = "36LD_3xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	if (3 != stc.get_nb_emission_antenna())
	{
		std::cout << "Number of emission antenna is wrong" << std::endl;
	}
	if (4 != stc.get_channel_uses())
	{
		std::cout << "Channel uses is wrong" << std::endl;
	}
	if (4 != stc.get_nb_symbols_per_block())
	{
		std::cout << "Number of symbols per block is wrong" << std::endl;
	}

	std::cout << "Checking 37 LD code 3xN" << std::endl;
	code_name = "37LD_3xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	if (3 != stc.get_nb_emission_antenna())
	{
		std::cout << "Number of emission antenna is wrong" << std::endl;
	}
	if (6 != stc.get_channel_uses())
	{
		std::cout << "Channel uses is wrong" << std::endl;
	}
	if (6 != stc.get_nb_symbols_per_block())
	{
		std::cout << "Number of symbols per block is wrong" << std::endl;
	}

	std::cout << "Checking 39 LD code 3xN" << std::endl;
	code_name = "39LD_3xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	if (3 != stc.get_nb_emission_antenna())
	{
		std::cout << "Number of emission antenna is wrong" << std::endl;
	}
	if (6 != stc.get_channel_uses())
	{
		std::cout << "Channel uses is wrong" << std::endl;
	}
	if (6 != stc.get_nb_symbols_per_block())
	{
		std::cout << "Number of symbols per block is wrong" << std::endl;
	}

	std::cout << "Space Time Code class test finished" << std::endl;
	return EXIT_SUCCESS;
}
