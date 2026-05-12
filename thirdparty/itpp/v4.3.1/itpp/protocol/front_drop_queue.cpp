/*!
 * \file
 * \brief Implementation of a Front Drop Queue class
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

#include <itpp/protocol/front_drop_queue.h>


namespace itpp
{

void Front_Drop_Queue::push(Packet *packet)
{
  if (debug) {
    std::cout << "Front_Drop_Queue::push_packet"
              //                << " byte_size=" << packet->bit_size()/8
              << " ptr=" << packet
              << " time=" << Event_Queue::now() << std::endl;
  }

  Packet *hol_packet;
  while ((!std::queue<Packet*>::empty()) &&
         ((8*bytes_in_queue + packet->bit_size()) >  8*max_bytes_in_queue)) {
    hol_packet = std::queue<Packet*>::front();
    Front_Drop_Queue::pop();
    delete hol_packet;

    //      TTCPPacket *tcp_packet = (TTCPPacket *) hol_packet;
    //      delete tcp_packet;

    if (debug) {
      std::cout << "Link_With_Input_Q::received_packet, "
                << "Packet Dropped, buffer overflow."
                << std::endl;
    }
  }

  bytes_in_queue += packet->bit_size() / 8;
  std::queue<Packet*>::push(packet);

}

void Front_Drop_Queue::pop()
{
  Packet *hol_packet;
  hol_packet = std::queue<Packet*>::front();
  bytes_in_queue -= (hol_packet->bit_size() / 8);
  if (debug) {
    std::cout << "Front_Drop_Queue::pop_packet"
              << " ptr=" << hol_packet
              << " time=" << Event_Queue::now() << std::endl;
  }
  std::queue<Packet*>::pop();

}

} // namespace itpp
