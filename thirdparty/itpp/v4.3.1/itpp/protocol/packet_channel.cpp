/*!
 * \file
 * \brief Implementation of a Packet channel class
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

#include <itpp/protocol/packet_channel.h>
#include <itpp/base/random.h>
#include <itpp/base/sort.h>
#include <itpp/base/math/min_max.h>


namespace itpp
{

Packet_Channel::Packet_Channel()
{
  parameters_ok = false;
  keep_running = false;
}

Packet_Channel::Packet_Channel(const double Pr, const Ttype Delay, const double Block_rate, const int Max_slots)
{
  set_parameters(Pr, Delay, Block_rate, Max_slots);
}


Packet_Channel::~Packet_Channel() {}

void Packet_Channel::set_parameters(const double Pr,  const Ttype Delay,  const double Block_rate, const int Max_slots)
{
  it_assert(Delay >= 0, "Packet_Channel::set_parameters(): ");
  it_assert(Pr >= 0.0 && Pr <= 1.0, "Packet_Channel::set_parameters(): ");
  it_assert(Block_rate > 0, "Packet_Channel::set_parameters(): ");
  it_assert(Max_slots >= 0, "Packet_Channel::set_parameters(): ");
  delay = Delay;
  pr = Pr;
  block_time = 1.0 / Block_rate;
  max_slots = Max_slots;
  input.forward(this, &Packet_Channel::handle_input);
  nof_inputs.forward(this, &Packet_Channel::handle_nof_inputs);
  start.forward(this, &Packet_Channel::handle_start);
  keep_running = false;
  explicit_errors = false;
  K = 0;
  k = 0;
  parameters_ok = true;
}

void Packet_Channel::handle_input(Link_Packet* M)
{
  it_assert(parameters_ok, "Packet_Channel::handle_input(): ");
  it_assert(M != NULL, "Packet_Channel::handle_input(): ");
  if (explicit_errors) {
    if (k < L) {
      lose = lost(k) == K;
      if (lose)
        k++;
    }
    K++;
  }
  else
    lose = randu() < pr;
  if (lose) {
    delete M;
  }
  else
    output(M, delay);
  lose = false;
}

void Packet_Channel::block_rate_loop()
{
  it_assert(parameters_ok, "Packet_Channel::block_rate_loop(): ");
  get_nof_inputs(NULL);
  if (keep_running)
    Event_Queue::add(new Event<Packet_Channel>(this, &Packet_Channel::block_rate_loop, block_time));
}

void Packet_Channel::handle_start(const bool run)
{
  it_assert(parameters_ok, "Packet_Channel::handle_start(): ");
  if (run && !keep_running)// Channel is in 'stop' state. Start it and keep running.
    Event_Queue::add(new Event<Packet_Channel>(this, &Packet_Channel::block_rate_loop, block_time));
  keep_running = run;
}

void Packet_Channel::handle_nof_inputs(const int Nof_ready_messages)
{
  it_assert(Nof_ready_messages >= 0, "Packet_Channel::handle_nof_inputs(): ");
  int L = 0;
  if (max_slots > 0)
    L = std::min(Nof_ready_messages, round_i(randu() * max_slots));
  else
    L = std::min(Nof_ready_messages, 1);
  if (L > 0)
    input_request(L);
}

void Packet_Channel::set_errors(const ivec &Lost)
{
  L = Lost.length();
  if (L > 0) {
    it_assert(min(Lost) >= 0, "Packet_Channel::set_errors(): ");
    lost = Lost;
    sort(lost);
    explicit_errors = true;
  }
}


// ----------------------------- Ack_Channel --------------------------------


ACK_Channel::ACK_Channel()
{
  parameters_ok = false;
}

ACK_Channel::ACK_Channel(const double Pr, const Ttype Delay)
{
  set_parameters(Pr, Delay);
}


ACK_Channel::~ACK_Channel() {}

void ACK_Channel::set_parameters(const double Pr, const Ttype Delay)
{
  it_assert(Delay >= 0, "ACK_Channel::set_parameters(): ");
  it_assert(Pr >= 0.0 && Pr <= 1.0, "ACK_Channel::set_parameters(): ");
  delay = Delay;
  pr = Pr;
  input.forward(this, &ACK_Channel::handle_input);
  explicit_errors = false;
  K = 0;
  k = 0;
  parameters_ok = true;
}

void ACK_Channel::handle_input(ACK* M)
{
  it_assert(parameters_ok, "ACK_Channel::handle_input(): ");
  it_assert(M != NULL, "ACK_Channel::handle_input(): ");
  if (explicit_errors) {
    if (k < L) {
      lose = lost(k) == K;
      if (lose)
        k++;
    }
    K++;
  }
  else
    lose = randu() < pr;
  if (lose)
    delete M;
  else
    output(M, delay);
  lose = false;
}

void ACK_Channel::set_errors(const ivec& Lost)
{
  L = Lost.length();
  if (L > 0) {
    it_assert(min(Lost) >= 0, "ACK_Channel::set_errors(): ");
    lost = Lost;
    sort(lost);
    explicit_errors = true;
  }
}

} // namespace itpp
