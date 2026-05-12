/*!
 * \file
 * \brief Definition of a Packet class
 * \author Krister Norlund, Anders Persson and Tony Ottosson
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

#ifndef PACKET_H
#define PACKET_H

#include <itpp/itexports.h>

#if (defined(_MSC_VER) && defined(ITPP_SHARED_LIB) && !(defined(itpp_EXPORTS) || defined(itpp_debug_EXPORTS)))

#ifndef ITPP_PROTOCOL_EXCLUDED
#define ITPP_PROTOCOL_EXCLUDED
#pragma message( "PROTOCOL definitions are not available for MSVC shared builds" )
#endif

#else

#include <itpp/protocol/signals_slots.h>


namespace itpp
{

//! \addtogroup protocol
//@{

/*! ADD DOCUMENTATION HERE

 */
class Packet
{
public:
  //! ADD DOCUMENTATION HERE
  Packet(const int packet_size = 0) { set_bit_size(packet_size); }
  //! ADD DOCUMENTATION HERE
  virtual ~Packet() {}

  //! set size of packet in bits
  void set_bit_size(int packet_size) { it_assert(packet_size >= 0, "Packet size must be positive"); size_bits = packet_size; }

  //! get size of packet in bits
  int bit_size() { return size_bits; }

private:
  int size_bits; // size of packet in bits
};


/*! ADD DOCUMENTATION HERE

 */
class L3_Packet_Info
{
public:
  //! ADD DOCUMENTATION HERE
  L3_Packet_Info(Packet *packet) { timestamp = 0; pkt_pointer = packet; }

  //! ADD DOCUMENTATION HERE
  ~L3_Packet_Info() {}

  //! ADD DOCUMENTATION HERE
  Ttype timestamp;

  //! ADD DOCUMENTATION HERE
  Packet *pkt_pointer;
};


/*! ADD DOCUMENTATION HERE

 */
class Link_Packet : public Packet
{
public:
  //! ADD DOCUMENTATION HERE
  Link_Packet(const int Seq_no, const unsigned long int Link_packet_id, L3_Packet_Info *Cp) { seq_no = Seq_no; link_packet_id = Link_packet_id; l3_pkt_info_p = Cp; }

  //! ADD DOCUMENTATION HERE
  ~Link_Packet() {}

  //! ADD DOCUMENTATION HERE
  unsigned long int link_packet_id;
  //! ADD DOCUMENTATION HERE
  int seq_no;
  //! ADD DOCUMENTATION HERE
  L3_Packet_Info *l3_pkt_info_p;
};

/*! ADD DOCUMENTATION HERE

 */
class ACK : public Packet
{
public:
  //! ADD DOCUMENTATION HERE
  ACK(const int Seq_no = -1, const int Id = 0) { seq_no = Seq_no; id = Id; }

  //! ADD DOCUMENTATION HERE
  ~ACK() {}

  //! ADD DOCUMENTATION HERE
  int id;
  //! ADD DOCUMENTATION HERE
  int seq_no;
};

//@}

} // namespace itpp

#endif

#endif // #ifndef PACKET_H

