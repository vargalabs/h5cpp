/*!
 * \file
 * \brief Implementation of a Packet generator class
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

#include <itpp/protocol/packet_generator.h>


namespace itpp
{

Packet_Generator::Packet_Generator(const int Packet_size, const unsigned long int Max_packets)
{
  keep_running = false;
  start.forward(this, &Packet_Generator::handle_start);
  next.forward(this, &Packet_Generator::handle_next);
  output.connect(&next);
  set_parameters(Packet_size, Max_packets);
}

Packet_Generator::~Packet_Generator() { }

void Packet_Generator::set_parameters(const int Packet_size, const unsigned long int Max_packets)
{
  it_assert(Packet_size > 0, "Packet_Generator::set_parameters(): ");
  packet_size = Packet_size;
  max_packets = Max_packets;
  id = 0;
}

int Packet_Generator::get_packet_size()
{
  return packet_size;
}

int Packet_Generator::get_max_packets()
{
  return max_packets;
}

void Packet_Generator::handle_next(Packet*)
{
  if (keep_running) {
    output(new Packet(8*packet_size), delta_t());
    id++;
    if (max_packets && id >= max_packets)
      start(false);
  }
}

void Packet_Generator::handle_start(const bool run)
{
  if (run && !keep_running) {
    keep_running = run;
    handle_next(NULL);
  }
  keep_running = run;
}


// ---------------------------- Poisson_Packet_Generator -------------------------------------------------

Poisson_Packet_Generator::Poisson_Packet_Generator(const double Avg_bit_rate,
    const int Packet_size,
    const unsigned long int Max_packets): Packet_Generator(Packet_size, Max_packets)
{
  set_parameters(Avg_bit_rate, Packet_size, Max_packets);
}

Poisson_Packet_Generator::~Poisson_Packet_Generator() {}

void Poisson_Packet_Generator::set_parameters(const double Avg_bit_rate,
    const int Packet_size,
    const unsigned long int Max_packets)
{
  Packet_Generator::set_parameters(Packet_size, Max_packets);
  it_assert(Avg_bit_rate > 0.0, "Packet_Generator::set_parameters(): ");
  avg_bit_rate = Avg_bit_rate;
  avg_delta_t = 8.0 * get_packet_size() / avg_bit_rate;
  ee.setup(1.0);
}

double Poisson_Packet_Generator::get_avg_bit_rate()
{
  return avg_bit_rate;
}


Ttype Poisson_Packet_Generator::delta_t()
{
  return ee()*avg_delta_t;
}


// ---------------------------- Constant_Rate_Packet_Generator -------------------------------------------------

Constant_Rate_Packet_Generator::Constant_Rate_Packet_Generator(const double Avg_rate, const int Packet_size, const unsigned long int Max_packets): Poisson_Packet_Generator(Avg_rate, Packet_size, Max_packets) {}

Constant_Rate_Packet_Generator::~Constant_Rate_Packet_Generator() {}

Ttype Constant_Rate_Packet_Generator::delta_t()
{
  return avg_delta_t;
}


// ---------------------------- Burst_WWW_Packet_Generator -------------------------------------------------


Burst_WWW_Packet_Generator::Burst_WWW_Packet_Generator(const double Avg_bit_rate, const int Packet_size, const int Max_packets): Poisson_Packet_Generator(Avg_bit_rate, Packet_size, Max_packets)
{
  Navg = 50; // Average number of packets per burst [packets].
  Ti = 1.1960e-4; // Average inter-arrival time between packets in burst [s].
  Tr = Navg * Packet_size * 8.0 / Avg_bit_rate - Ti * (Navg - 1); // Average time between bursts.
  N = 0;
}

Burst_WWW_Packet_Generator::~Burst_WWW_Packet_Generator()
{

}

Ttype Burst_WWW_Packet_Generator::delta_t()
{
  if (N == 0) { // Start of a new burst.
    N = Navg;
    N--; // First packet is triggered at ...
    return ee()*Tr; // ... start time of next burst.
  }
  else { // Within a burst.
    N--; // One packet less in the burst ...
    return ee()*Ti; // ... arrival time for next packet within the burst.
  }
}


// ----------------------------Sink -------------------------------------------------

Sink::Sink(const unsigned long int Max_packets)
{
  it_assert(Max_packets > 0, "Sink::Sink(): ");
  max_packets = Max_packets;
  Ncp = 0;
  Nbytes = 0;
  packet_input.forward(this, &Sink::handle_packet_input);
  start_time = Event_Queue::now();
}

Sink::~Sink()
{
  std::cout << "Time = " << Event_Queue::now() << ", Sink : " << std::endl;
  std::cout << "Received " << Ncp << " packets in sequence." << std::endl;
  std::cout << "Receive average bit rate = " << Nbytes*8.0 / (Event_Queue::now() - start_time) << " [bits/second]." << std::endl;
}


void Sink::handle_packet_input(Packet *P)
{
  it_assert(P != NULL, "Sink::handle_packet_input(): ");
  Ncp++;
  Nbytes += (P->bit_size() / 8);
  delete P;
  if (Ncp >= max_packets) {
    std::cout << "Time = " << Event_Queue::now() << ", Sink : " << std::endl;
    std::cout << "Simulation stopped because : Ncp > max_packets" << std::endl;
    Event_Queue::stop();
  }
}


} // namespace itpp
