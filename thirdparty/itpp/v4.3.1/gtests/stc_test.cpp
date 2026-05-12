/*!
 * \file
 * \brief STC class test program
 * \author Bogdan Cristea
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2013  (see AUTHORS file for a list of contributors)
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
#include "gtest/gtest.h"

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

TEST(STC, All)
{
	int em_antenna = 2;
	int channel_uses = 2;
	const int const_size = 4;//QAM
	const double abs_err = 1e-12;

	// Checking V-BLAST
	std::string code_name("V-BLAST_MxN");

	//should succeed
	itpp::STC stc(code_name, const_size, em_antenna, channel_uses);

	//check parameters
	ASSERT_EQ(2, stc.get_nb_emission_antenna());
	ASSERT_EQ(2, stc.get_channel_uses());
	ASSERT_EQ(4, stc.get_nb_symbols_per_block());

	//check generated matrices
	itpp::cmat ref_vblast("1 0; 0 0; 0 1; 0 0; 0 0; 1 0; 0 0; 0 1");
	ASSERT_TRUE(near(ref_vblast, stc.get_1st_gen_matrix(), abs_err));
	ASSERT_TRUE(near(ref_vblast, stc.get_2nd_gen_matrix(), abs_err));

	// Checking improved V-BLAST
	code_name = "imp_V-BLAST_MxN";

	//should succeed
	stc.setup(code_name, const_size, em_antenna, channel_uses);

	//check parameters
	ASSERT_EQ(2, stc.get_nb_emission_antenna());
	ASSERT_EQ(2, stc.get_channel_uses());
	ASSERT_EQ(4, stc.get_nb_symbols_per_block());

	//check generated matrices
	itpp::cmat ref_imp_vblast("1 0; 0 1; 0 1; 1 0; 1 0; 0 -1; 0 1; -1 0");
	ref_imp_vblast /= std::sqrt(2.0);
	ASSERT_TRUE(near(ref_imp_vblast, stc.get_1st_gen_matrix(), abs_err));
	ASSERT_TRUE(near(ref_imp_vblast, stc.get_2nd_gen_matrix(), abs_err));

	// Checking Alamouti code
	code_name = "Alamouti_2xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	ASSERT_EQ(2, stc.get_nb_emission_antenna());
	ASSERT_EQ(2, stc.get_channel_uses());
	ASSERT_EQ(2, stc.get_nb_symbols_per_block());

	//check generated matrices
	itpp::cmat ref_A("1 0; 0 1; 0 1; -1 0");
	itpp::cmat ref_B("1 0; 0 -1; 0 1; 1 0");
	ASSERT_TRUE(near(ref_A, stc.get_1st_gen_matrix(), abs_err));
	ASSERT_TRUE(near(ref_B, stc.get_2nd_gen_matrix(), abs_err));

	// Checking switched Alamouti code
	code_name = "Switched_Alamouti_4xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	ASSERT_EQ(4, stc.get_nb_emission_antenna());
	ASSERT_EQ(4, stc.get_channel_uses());
	ASSERT_EQ(4, stc.get_nb_symbols_per_block());

	// Checking double Alamouti code
	code_name = "Double_Alamouti_4xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	ASSERT_EQ(4, stc.get_nb_emission_antenna());
	ASSERT_EQ(2, stc.get_channel_uses());
	ASSERT_EQ(4, stc.get_nb_symbols_per_block());

	// Checking Jafarkhani code
	code_name = "Jafarkhani_4xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	ASSERT_EQ(4, stc.get_nb_emission_antenna());
	ASSERT_EQ(4, stc.get_channel_uses());
	ASSERT_EQ(4, stc.get_nb_symbols_per_block());

	// Checking Golden code
	code_name = "Golden_2x2";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	ASSERT_EQ(2, stc.get_nb_emission_antenna());
	ASSERT_EQ(2, stc.get_channel_uses());
	ASSERT_EQ(4, stc.get_nb_symbols_per_block());

	// Checking Damen code
	code_name = "Damen_2x2";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	ASSERT_EQ(2, stc.get_nb_emission_antenna());
	ASSERT_EQ(2, stc.get_channel_uses());
	ASSERT_EQ(4, stc.get_nb_symbols_per_block());

	// Checking orthogonal code 3xN
	code_name = "34ortho_3xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	ASSERT_EQ(3, stc.get_nb_emission_antenna());
	ASSERT_EQ(4, stc.get_channel_uses());
	ASSERT_EQ(3, stc.get_nb_symbols_per_block());

	// Checking 36 LD code 3xN
	code_name = "36LD_3xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	ASSERT_EQ(3, stc.get_nb_emission_antenna());
	ASSERT_EQ(4, stc.get_channel_uses());
	ASSERT_EQ(4, stc.get_nb_symbols_per_block());

	// Checking 37 LD code 3xN
	code_name = "37LD_3xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	ASSERT_EQ(3, stc.get_nb_emission_antenna());
	ASSERT_EQ(6, stc.get_channel_uses());
	ASSERT_EQ(6, stc.get_nb_symbols_per_block());

	// Checking 39 LD code 3xN
	code_name = "39LD_3xN";

	//should succeed
	stc.setup(code_name, const_size);

	//check parameters
	ASSERT_EQ(3, stc.get_nb_emission_antenna());
	ASSERT_EQ(6, stc.get_channel_uses());
	ASSERT_EQ(6, stc.get_nb_symbols_per_block());
}
