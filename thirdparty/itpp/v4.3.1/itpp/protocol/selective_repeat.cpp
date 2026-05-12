/*!
 * \file
 * \brief Implementation of Selective Repeat ARQ classes
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

#include <itpp/protocol/selective_repeat.h>
#include <cstdlib>

//! \cond

namespace itpp
{

bool in_sequence(const int a, const int b, const int L)
{
  it_assert(a >= 0 && a < L, "in_sequence(): ");
  it_assert(b >= 0 && b < L, "in_sequence(): ");
  return ((b - a + L) % L) < L / 2;
}

Selective_Repeat_ARQ_Sender::Selective_Repeat_ARQ_Sender()
{
  parameters_ok = false;
  packet_input.set_name("Selective_Repeat_ARQ_Sender packet_input Slot");
  packet_input.forward(this, &Selective_Repeat_ARQ_Sender::handle_packet_input);
  ack_input.set_name("Selective_Repeat_ARQ_Sender ack_input Slot");
  ack_input.forward(this, &Selective_Repeat_ARQ_Sender::handle_ack_input);
  query_nof_ready_packets.set_name("Selective_Repeat_ARQ_Sender query_nof_ready_packets Slot");
  query_nof_ready_packets.forward(this, &Selective_Repeat_ARQ_Sender::handle_query_nof_ready_packets);
  packet_output_request.set_name("Selective_Repeat_ARQ_Sender packet_output_request Slot");
  packet_output_request.forward(this, &Selective_Repeat_ARQ_Sender::handle_packet_output_request);

}

Selective_Repeat_ARQ_Sender::Selective_Repeat_ARQ_Sender(const int Seq_no_size, const int Buffer_size_factor, const int Link_packet_size, const Ttype Time_out)
{
  set_parameters(Seq_no_size, Buffer_size_factor, Link_packet_size, Time_out);
  packet_input.set_name("Selective_Repeat_ARQ_Sender packet_input Slot");
  packet_input.forward(this, &Selective_Repeat_ARQ_Sender::handle_packet_input);
  ack_input.set_name("Selective_Repeat_ARQ_Sender ack_input Slot");
  ack_input.forward(this, &Selective_Repeat_ARQ_Sender::handle_ack_input);
  query_nof_ready_packets.set_name("Selective_Repeat_ARQ_Sender query_nof_ready_packets Slot");
  query_nof_ready_packets.forward(this, &Selective_Repeat_ARQ_Sender::handle_query_nof_ready_packets);
  packet_output_request.set_name("Selective_Repeat_ARQ_Sender packet_output_request Slot");
  packet_output_request.forward(this, &Selective_Repeat_ARQ_Sender::handle_packet_output_request);
}

Selective_Repeat_ARQ_Sender::~Selective_Repeat_ARQ_Sender()
{
  std::cout << "no_retransmit = " << no_retransmit << std::endl;
}

void Selective_Repeat_ARQ_Sender::set_parameters(const int Seq_no_size,
    const int Buffer_size_factor,
    const int Link_packet_size,
    const Ttype Time_out)
{
  it_assert((0 < Seq_no_size) && (Seq_no_size <= 30),
            "Selective_Repeat_ARQ_Sender::set_parameters(): ");
  it_assert((0 < Buffer_size_factor) && (Buffer_size_factor <= 10),
            "Selective_Repeat_ARQ_Sender::set_parameters(): ");
  it_assert(Link_packet_size > 0, "Selective_Repeat_ARQ_Sender::set_parameters(): ");
  it_assert(Time_out > 0, "Selective_Repeat_ARQ_Sender::set_parameters(): ");
  seq_no_size = Seq_no_size;
  link_packet_size = Link_packet_size;
  seq_no_max = 1 << Seq_no_size;
  input_buffer_size = seq_no_max * Buffer_size_factor;
  input_buffer.set_size(input_buffer_size);
  for (int l = 0; l < input_buffer_size; input_buffer(l++) = NULL);
  input_free_space = input_buffer_size;
  input_next = 0;
  tx_next = 0;
  tx_last = 0;
  time_out = Time_out;
  timer.set_size(seq_no_max);
  for (int l = 0; l < seq_no_max; timer(l++).forward(this, &Selective_Repeat_ARQ_Sender::retransmit));
  outstanding = 0;
  seq_no = 0;
  output_indexes.set_size(seq_no_max);
  output_indexes.ones();
  output_indexes *= -1;
  retransmission_indexes.set_size(seq_no_max);
  retransmission_indexes.ones();
  retransmission_indexes *= -1;
  rd_pos = 0;
  rt_pos = 0;
  scheduled_total = 0;
  scheduled_retransmissions = 0;
  no_retransmit = 0;
  parameters_ok = true;
  ip_pkt_queue.set_max_byte_size(1500*32);
  id = 0;
}

void Selective_Repeat_ARQ_Sender::handle_ack_input(Array<Packet*> packet_array)
{
  Packet *packet = packet_array(0);
  ACK *A = (ACK *) packet;

  it_assert(parameters_ok, "Selective_Repeat_ARQ_Sender::handle_ack_input(): ");
  it_assert(A, "Selective_Repeat_ARQ_Sender::handle_ack_input(): ");
  it_assert(A->seq_no >= 0 && A->seq_no < seq_no_max, "Selective_Repeat_ARQ_Sender::handle_ack_input(): ");
  if (outstanding) {
    if (in_sequence(tx_last % seq_no_max, A->seq_no, seq_no_max))
      remove(A->seq_no);
    while (!input_buffer(tx_last) && outstanding) {
      outstanding--;
      input_free_space++;
      tx_last = (tx_last + 1) % input_buffer_size;
    }
  }
  delete A;
  fill_output();
}

void Selective_Repeat_ARQ_Sender::handle_packet_input(Packet *packet)
{
  it_assert(parameters_ok, "Selective_Repeat_ARQ_Sender::handle_packet_input(): ");
  it_assert(packet, "Selective_Repeat_ARQ_Sender::handle_packet_input(): ");
  ip_pkt_queue.push(packet);

}

// The number of blocks in the ip_pkt_queue that can be scheduled to be
// transmitted (in the tx buffer)
int Selective_Repeat_ARQ_Sender::feasable_blocks()
{
  div_t q = div(ip_pkt_queue.byte_size(), link_packet_size);
  int blocks_in_ip_queue = (q.rem) ? q.quot + 1 : q.quot;
  return std::min(free_sequence_numbers(),
                  buffered_non_outstanding() +
                  std::min(blocks_in_ip_queue, input_free_space));
}


void Selective_Repeat_ARQ_Sender::handle_query_nof_ready_packets(void*)
{
  it_assert(parameters_ok, "Selective_Repeat_ARQ_Sender::handle_query_nof_ready_packets(): ");
  nof_ready_packets(scheduled_total + feasable_blocks());
}

void Selective_Repeat_ARQ_Sender::handle_packet_output_request(const int nbr_blocks_requested)
{
  int nbr_blocks_to_tx;
  int feasable_blks = feasable_blocks();
  if (nbr_blocks_requested <= scheduled_total + feasable_blks) {
    nbr_blocks_to_tx = nbr_blocks_requested;
  }
  else {
    it_warning("Number of requested blocks is more than what is possible to transmitt");
    nbr_blocks_to_tx = scheduled_total + feasable_blks;
  }

  //int nbr_ip_pkts_in_q = ip_pkt_queue.size();
  while (nbr_blocks_to_tx > scheduled_total) {
    it_assert(!ip_pkt_queue.empty(), "Selective_Repeat_ARQ_Sender::handle_packet_output_request(): ");
    Packet *packet = ip_pkt_queue.front();
    ip_pkt_queue.pop();
    push_packet_on_tx_buffer(packet);
  }

  Array<Packet*> tmp;
  get_link_packets(nbr_blocks_requested, tmp);
  packet_output(tmp);
}

void Selective_Repeat_ARQ_Sender::push_packet_on_tx_buffer(Packet *packet)
{
  L3_Packet_Info *pkt_info = new L3_Packet_Info(packet);
  int packet_byte_size = pkt_info->pkt_pointer->bit_size() / 8;
  int nbr_blocks = packet_byte_size / link_packet_size;
  if (nbr_blocks*link_packet_size != packet_byte_size)
    nbr_blocks++;
  if (input_free_space >= nbr_blocks) {
    pkt_info->timestamp = Event_Queue::now();
    for (int n = nbr_blocks - 1; n >= 0; n--) {
      input_buffer(input_next) = new Link_Packet(-1, n, pkt_info);
      input_free_space--;
      input_next = (input_next + 1) % input_buffer_size;
    }
  }
  else {
    buffer_overflow(0);
    it_error("Selective_Repeat_ARQ_Sender::push_packet_on_tx_buffer(): "
             "Stopped due to buffer overflow");
  }
  fill_output();

}

void Selective_Repeat_ARQ_Sender::fill_output()
{
  int packets_2_output = std::min(free_sequence_numbers(), buffered_non_outstanding());
  while (packets_2_output) {
    input_buffer(tx_next)->seq_no = seq_no;
    outstanding++;
    schedule_output(tx_next, seq_no, false);
    seq_no = (seq_no + 1) % seq_no_max;
    tx_next = (tx_next + 1) % input_buffer_size;
    packets_2_output--;
  }
}

void Selective_Repeat_ARQ_Sender::schedule_output(const int Buffer_index, const int Sequence_number, const bool Retransmission)
{
  it_assert(input_buffer(Buffer_index) != NULL, "Selective_Repeat_ARQ_Sender::schedule_output(): ");
  if (output_indexes(Sequence_number) == -1)
    scheduled_total++;
  output_indexes(Sequence_number) = Buffer_index;
  if (Retransmission) {
    if (retransmission_indexes(Sequence_number) != 1) // This is a new retransmission.
      scheduled_retransmissions++;
    retransmission_indexes(Sequence_number) = 1; // Mark packet (index) for retransmission.
  }
  else // Mark packet (index) for first time transmission.
    retransmission_indexes(Sequence_number) = 0;
}

void Selective_Repeat_ARQ_Sender::get_link_packets(const int K, Array<Packet*> &pa)
{
  int packets_2_retransmit = std::min(K, scheduled_retransmissions);
  int new_packets_2_transmit = std::min(K, scheduled_total) - packets_2_retransmit;
  scheduled_retransmissions -= packets_2_retransmit;
  scheduled_total -= packets_2_retransmit + new_packets_2_transmit;
  pa.set_size(packets_2_retransmit + new_packets_2_transmit);
  int l = 0;
  while (packets_2_retransmit) { // Retransmissions have priority over ...
    if (retransmission_indexes(rt_pos) == 1) {
      timer(rt_pos).set(rt_pos, time_out);
      pa(l++) = (Packet *) new Link_Packet(*input_buffer(output_indexes(rt_pos)));
      output_indexes(rt_pos) = -1;
      retransmission_indexes(rt_pos) = -1;
      packets_2_retransmit--;
    }
    rt_pos = (rt_pos + 1) % seq_no_max;
  }
  while (new_packets_2_transmit) { // new packets.
    if (output_indexes(rd_pos) != -1) {
      timer(rd_pos).set(rd_pos, time_out);
      pa(l++) = (Packet *) new Link_Packet(*input_buffer(output_indexes(rd_pos)));
      output_indexes(rd_pos) = -1;
      new_packets_2_transmit--;
    }
    rd_pos = (rd_pos + 1) % seq_no_max;
  }
}

void Selective_Repeat_ARQ_Sender::remove(const int Sequence_number)
{
  if (output_indexes(Sequence_number) != -1) {
    output_indexes(Sequence_number) = -1;
    scheduled_total--;
    if (retransmission_indexes(Sequence_number) == 1)
      scheduled_retransmissions--;
    retransmission_indexes(Sequence_number) = -1;
  }
  const int i = sequence_number_2_buffer_index(Sequence_number);
  if (input_buffer(i)) {
    timer(Sequence_number).cancel(); // Cancel the retransmission timer.
    it_assert(input_buffer(i)->seq_no == Sequence_number, "Selective_Repeat_ARQ_Sender::remove(): ");
    delete input_buffer(i);
    input_buffer(i) = NULL;
  }
}

void Selective_Repeat_ARQ_Sender::retransmit(const int Sequence_number)
{
  no_retransmit++;
  const int buffer_index = sequence_number_2_buffer_index(Sequence_number);
  schedule_output(buffer_index, Sequence_number, true);
}

int Selective_Repeat_ARQ_Sender::buffered_non_outstanding()
{
  return input_buffer_size - input_free_space - outstanding;
}

int Selective_Repeat_ARQ_Sender::free_sequence_numbers()
{
  return seq_no_max / 2 - outstanding;
}

int Selective_Repeat_ARQ_Sender::sequence_number_2_buffer_index(const int Sequence_number)
{
  it_assert(input_buffer(tx_last), "Selective_Repeat_ARQ_Sender::sequence_number_2_buffer_index(): ");
  it_assert(input_buffer(tx_last)->seq_no != -1, "Selective_Repeat_ARQ_Sender::sequence_number_2_buffer_index(): ");
  return (tx_last + (Sequence_number - input_buffer(tx_last)->seq_no + seq_no_max) % seq_no_max) % input_buffer_size;
}

int Selective_Repeat_ARQ_Sender::link_packets_buffered()
{
  it_assert(parameters_ok, "Selective_Repeat_ARQ_Sender::link_packets_buffered(): ");
  return input_buffer_size - input_free_space;
}

int Selective_Repeat_ARQ_Sender::nof_ready_link_packets()
{
  it_assert(parameters_ok, "Selective_Repeat_ARQ_Sender::nof_ready_link_packets(): ");
  return scheduled_total + feasable_blocks();
}

int Selective_Repeat_ARQ_Sender::link_packets_queued_waiting_for_transmission()
{
  it_assert(parameters_ok, "Selective_Repeat_ARQ_Sender::link_packets_queued_waiting_for_transmission(): ");
  div_t q = div(ip_pkt_queue.byte_size(), link_packet_size);
  int blocks_in_ip_queue = (q.rem) ? q.quot + 1 : q.quot;
  return buffered_non_outstanding() + scheduled_total + blocks_in_ip_queue;
}

// int Selective_Repeat_ARQ_Sender::time_stamp_HOL_packet(){
//   assert(parameters_ok);
//   return buffered_non_outstanding()+feasable_blocks();
// }

int Selective_Repeat_ARQ_Sender::buffer_size()
{
  it_assert(parameters_ok, "Selective_Repeat_ARQ_Sender::buffer_size(): ");
  return input_buffer_size;
}

Ttype Selective_Repeat_ARQ_Sender::link_packets_max_queuing_time()
{
  it_assert(parameters_ok, "Selective_Repeat_ARQ_Sender::link_packets_max_queuing_time(): ");
  it_assert(input_buffer(tx_last), "Selective_Repeat_ARQ_Sender::link_packets_max_queuing_time(): ");
  return Event_Queue::now() - input_buffer(tx_last)->l3_pkt_info_p->timestamp;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Selective_Repeat_ARQ_Receiver::Selective_Repeat_ARQ_Receiver()
{
  parameters_ok = false;
  packet_input.forward(this, &Selective_Repeat_ARQ_Receiver::handle_packet_input);
  packet_input.set_name("Selective_Repeat_ARQ_Receiver packet_input Slot");
}

Selective_Repeat_ARQ_Receiver::Selective_Repeat_ARQ_Receiver(const int Seq_no_size)
{
  set_parameters(Seq_no_size);
  packet_input.forward(this, &Selective_Repeat_ARQ_Receiver::handle_packet_input);
  packet_input.set_name("Selective_Repeat_ARQ_Receiver packet_input Slot");
}

Selective_Repeat_ARQ_Receiver::~Selective_Repeat_ARQ_Receiver() {}

void Selective_Repeat_ARQ_Receiver::set_parameters(const int Seq_no_size)
{
  seq_no_size = Seq_no_size;
  seq_no_max = 1 << seq_no_size;
  rx_buffer.set_size(seq_no_max);
  for (int l = 0; l < seq_no_max; rx_buffer(l++) = NULL);
  Rnext = 0;
  id = 0;
  parameters_ok = true;
}

void Selective_Repeat_ARQ_Receiver::handle_packet_input(Array<Packet*> packet_array)
{
  it_assert(parameters_ok, "Selective_Repeat_ARQ_Receiver::handle_packet_input(): ");

  int nbr_pkts = packet_array.length();
  Link_Packet *packet;
  for (int i = 0;i < nbr_pkts;i++) {
    packet = (Link_Packet *) packet_array(i);
    it_assert(packet, "Selective_Repeat_ARQ_Receiver::handle_packet_input(): ");
    it_assert(packet->seq_no >= 0 && packet->seq_no < seq_no_max, "Selective_Repeat_ARQ_Receiver::handle_packet_input(): ");
    Array<Packet*> ack_pkt;
    ack_pkt.set_size(1);
    ack_pkt(0) = (Packet *) new ACK(packet->seq_no, id++);
    ack_output(ack_pkt); // Acknowledge the receipt of this packet.
    if (in_sequence(Rnext, packet->seq_no, seq_no_max) && !rx_buffer(packet->seq_no)) // Is this a new packet in-sequence packet?
      rx_buffer(packet->seq_no) = packet; // This is a new in-sequence packet.
    else // This either is a duplicate packet or an out-of-sequence packet.
      delete packet;
    while (rx_buffer(Rnext)) { // Is there an unbroken sequence of packets that we can output?

      if (rx_buffer(Rnext)->link_packet_id == 0) {
        packet_output(rx_buffer(Rnext)->l3_pkt_info_p->pkt_pointer);
        delete rx_buffer(Rnext)->l3_pkt_info_p;
      }
      delete rx_buffer(Rnext);
      rx_buffer(Rnext) = NULL;
      Rnext = (Rnext + 1) % seq_no_max;
    }
  }
}


} //namespace itpp

//! \endcond
