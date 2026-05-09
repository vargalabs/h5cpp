/*!
 * \file
 * \brief Include file for the IT++ communications module
 * \author Tony Ottosson
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

#ifndef ITCOMM_H
#define ITCOMM_H

/*!
 * \defgroup comm Communications Module
 * @{
 */

//! \defgroup channels Channel Modeling
//! \defgroup modulators Digital Modulation
//! \defgroup fec Forward Error Correcting Codes
//! \defgroup interl Interleavers
//! \defgroup misccommfunc Miscellaneous Communications Functions
//! \defgroup sequence Sequences

/*!
 * @}
 */

#include <itpp/itsignal.h>
#include <itpp/comm/bch.h>
#include <itpp/comm/channel.h>
#include <itpp/comm/channel_code.h>
#include <itpp/comm/commfunc.h>
#include <itpp/comm/convcode.h>
#include <itpp/comm/crc.h>
#include <itpp/comm/egolay.h>
#include <itpp/comm/error_counters.h>
#include <itpp/comm/galois.h>
#include <itpp/comm/hammcode.h>
#include <itpp/comm/interleave.h>
#include <itpp/comm/ldpc.h>
#include <itpp/comm/llr.h>
#include <itpp/comm/modulator.h>
#include <itpp/comm/modulator_nd.h>
#include <itpp/comm/ofdm.h>
#include <itpp/comm/pulse_shape.h>
#include <itpp/comm/punct_convcode.h>
#include <itpp/comm/rec_syst_conv_code.h>
#include <itpp/comm/reedsolomon.h>
#include <itpp/comm/sequence.h>
#include <itpp/comm/spread.h>
#include <itpp/comm/turbo.h>
#include <itpp/comm/siso.h>
#include <itpp/comm/exit.h>
#include <itpp/comm/stc.h>
#include <itpp/comm/multilateration.h>

#endif // #ifndef ITCOMM_H
