/*!
 * \file
 * \brief Definition of a Packet generator class
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

#ifndef PACKET_GENERATOR_H
#define PACKET_GENERATOR_H

#include <itpp/base/random.h>

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

/*! ADD DOCUMENTATION HERE

 */
class Packet_Generator
{
public:
  //! ADD DOCUMENTATION HERE
  Packet_Generator(const int Packet_size = 150, const unsigned long int Max_packets = 0);
  //! ADD DOCUMENTATION HERE
  virtual ~Packet_Generator();
  //! ADD DOCUMENTATION HERE
  Signal<Packet*> output;
  //! ADD DOCUMENTATION HERE
  Slot<Packet_Generator, bool> start;
  //! ADD DOCUMENTATION HERE
  void set_parameters(const int Packet_size, const unsigned long int Max_packets);
  //! ADD DOCUMENTATION HERE
  int get_packet_size();
  //! ADD DOCUMENTATION HERE
  int get_max_packets();
protected:
  //! ADD DOCUMENTATION HERE
  virtual Ttype delta_t() = 0;
private:
  Slot<Packet_Generator, Packet*> next;
  void handle_next(Packet*);
  void handle_start(const bool run);
  bool keep_running;
  unsigned long int id;
  int packet_size;
  unsigned long int max_packets;
};


/*! ADD DOCUMENTATION HERE

 */
class Poisson_Packet_Generator : public Packet_Generator
{
public:
  //! ADD DOCUMENTATION HERE
  Poisson_Packet_Generator(const double Avg_bit_rate = 1.0, const int Packet_size = 150, const unsigned long int Max_packets = 0);
  //! ADD DOCUMENTATION HERE
  virtual ~Poisson_Packet_Generator();
  //! ADD DOCUMENTATION HERE
  void set_parameters(const double Avg_bit_rate, const int Packet_size, const unsigned long int Max_packets);
  //! ADD DOCUMENTATION HERE
  double get_avg_bit_rate();
protected:
  //! ADD DOCUMENTATION HERE
  virtual Ttype delta_t();
  //! ADD DOCUMENTATION HERE
  double avg_delta_t;
  //! ADD DOCUMENTATION HERE
  double avg_bit_rate;
  //! ADD DOCUMENTATION HERE
  Exponential_RNG ee;
};


/*! ADD DOCUMENTATION HERE

 */
class Constant_Rate_Packet_Generator : public Poisson_Packet_Generator
{
public:
  //! ADD DOCUMENTATION HERE
  Constant_Rate_Packet_Generator(const double Avg_bit_rate = 1.0, const int Packet_size = 150, const unsigned long int Max_packets = 0);
  //! ADD DOCUMENTATION HERE
  virtual ~Constant_Rate_Packet_Generator();
protected:
  //! ADD DOCUMENTATION HERE
  virtual Ttype delta_t();
};

/*! ADD DOCUMENTATION HERE

 */
class Burst_WWW_Packet_Generator : public Poisson_Packet_Generator
{
public:
  //! ADD DOCUMENTATION HERE
  Burst_WWW_Packet_Generator(const double Avg_bit_rate = 1.0, const int Packet_size = 150, const int Max_packets = 0);
  //! ADD DOCUMENTATION HERE
  virtual ~Burst_WWW_Packet_Generator();
protected:
  //! ADD DOCUMENTATION HERE
  virtual Ttype delta_t();
  //! ADD DOCUMENTATION HERE
  int N;
  //! ADD DOCUMENTATION HERE
  int Navg;
  //! ADD DOCUMENTATION HERE
  double Ti;
  //! ADD DOCUMENTATION HERE
  double Tr;
};


/*! ADD DOCUMENTATION HERE

 */
class Sink
{
public:
  //! ADD DOCUMENTATION HERE
  Sink(const unsigned long int Max_packets = 1000);
  //! ADD DOCUMENTATION HERE
  ~Sink();
  // -- Slots -- //
  //! ADD DOCUMENTATION HERE
  Slot<Sink, Packet*> packet_input;
private:
  void handle_packet_input(Packet* packet);
  unsigned long int Ncp;
  unsigned long int Nbytes;
  unsigned long int max_packets;
  Ttype start_time;
};

//@}

} // namespace itpp

#endif

#endif // #ifndef PACKET_GENERATOR_H

