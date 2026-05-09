/*!
 * \file
 * \brief Definition of a Packet channel class
 * \author Anders Persson and Tony Ottosson
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

#ifndef PACKET_CHANNEL_H
#define PACKET_CHANNEL_H

#include <itpp/base/vec.h>

#if (defined(_MSC_VER) && defined(ITPP_SHARED_LIB) && !(defined(itpp_EXPORTS) || defined(itpp_debug_EXPORTS)))

#ifndef ITPP_PROTOCOL_EXCLUDED
#define ITPP_PROTOCOL_EXCLUDED
#pragma message( "PROTOCOL definitions are not available for MSVC shared builds" )
#endif

#else

#include <itpp/protocol/packet.h>

namespace itpp
{

//! \addtogroup protocol
//@{

//! ADD DOCUMENTATION HERE
class Packet_Channel
{
public:
  //! ADD DOCUMENTATION HERE
  Packet_Channel();
  //! ADD DOCUMENTATION HERE
  Packet_Channel(const double Pr, const Ttype Delay, const double Block_rate, const int Max_slots = 0);

  //! ADD DOCUMENTATION HERE
  ~Packet_Channel();

  // -- Slots -- //
  //! ADD DOCUMENTATION HERE
  Slot<Packet_Channel, bool> start;
  //! ADD DOCUMENTATION HERE
  Slot<Packet_Channel, Link_Packet*> input;
  //! ADD DOCUMENTATION HERE
  Slot<Packet_Channel, int> nof_inputs;

  // -- Signals -- //
  //! ADD DOCUMENTATION HERE
  Signal<Link_Packet*> output;
  //! ADD DOCUMENTATION HERE
  Signal<int> input_request;
  //! ADD DOCUMENTATION HERE
  Signal<void*> get_nof_inputs;

  //! ADD DOCUMENTATION HERE
  void set_parameters(const double Pr, const Ttype Delay, const double Block_rate, const int Max_slots);

  //! ADD DOCUMENTATION HERE
  void set_errors(const ivec &Lost);

private:
  void block_rate_loop();
  void handle_input(Link_Packet* M);
  void handle_start(const bool start);
  void handle_nof_inputs(const int N);

  bool keep_running;
  bool parameters_ok;
  bool explicit_errors;
  bool lose;
  double pr;
  Ttype delay;
  double block_time;
  int max_slots;
  ivec lost;
  int k, K, L;
};


//! ADD DOCUMENTATION HERE
class ACK_Channel
{
public:
  //! ADD DOCUMENTATION HERE
  ACK_Channel();

  //! ADD DOCUMENTATION HERE
  ACK_Channel(const double Pr, const Ttype Delay);

  //! ADD DOCUMENTATION HERE
  ~ACK_Channel();

  // -- Slots -- //
  //! ADD DOCUMENTATION HERE
  Slot<ACK_Channel, ACK*> input;

  // -- Signals -- //
  //! ADD DOCUMENTATION HERE
  Signal<ACK*> output;

  //! ADD DOCUMENTATION HERE
  void set_parameters(const double Pr, const Ttype Delay);
  //! ADD DOCUMENTATION HERE
  void set_errors(const ivec& Lost);

private:
  void handle_input(ACK* M);

  bool parameters_ok;
  bool explicit_errors;
  bool lose;
  double pr;
  Ttype delay;
  ivec lost;
  int k, K, L;
};

//@}

} // namespace itpp

#endif

#endif // #ifndef PACKET_CHANNEL_H

