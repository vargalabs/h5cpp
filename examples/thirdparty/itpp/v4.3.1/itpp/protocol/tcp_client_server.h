/*!
 * \file
 * \brief Definitions of TCP Client and Server Applications
 * \author Krister Norlund
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
#ifndef TCP_CLIENT_SERVER_H
#define TCP_CLIENT_SERVER_H

#include <itpp/itexports.h>

#if (defined(_MSC_VER) && defined(ITPP_SHARED_LIB) && !defined(itpp_EXPORTS))

#ifndef ITPP_PROTOCOL_EXCLUDED
#define ITPP_PROTOCOL_EXCLUDED
#pragma message( "PROTOCOL definitions are not available for MSVC shared builds" )
#endif

#else

#include <itpp/protocol/tcp.h>


namespace itpp
{

//! \addtogroup protocol
//@{

/*! ADD DOCUMENTATION HERE

 */
class TCP_Server_Application
{
public:
  //! Default constructor
  TCP_Server_Application() {
    write.set_name("TcpServerApplicationWriteSignal");
    write.set_debug();
  }
  //! Destructor
  ~TCP_Server_Application() { }

  //! ADD DOCUMENTATION HERE
  Signal<itpp::Packet*> write;

  //! ADD DOCUMENTATION HERE
  void write_to_net(unsigned byte_size, double delta_time) {
    itpp::Packet *packet = new Packet(8*byte_size);
    write(packet, delta_time);

    std::cout << "TcpServerApplication::write_to_net,"
              << " byte_size=" << packet->bit_size() / 8
              << " ptr=" << packet
              << " time=" << Event_Queue::now() << std::endl;
  }
};

/*! ADD DOCUMENTATION HERE

 */
class TCP_Client_Application
{
public:
  //! Default constructor
  TCP_Client_Application(TCP_Sender *tcp_snd_p, TCP_Receiver *tcp_recv_p) {
    tcp_receiver_p = tcp_recv_p;
    tcp_sender_p = tcp_snd_p;
    nbr_bytes_received = 0;
    select.forward(this, &TCP_Client_Application::received_packet_indication);
    select.set_name("TcpClientApplicationSelectSlot");
    seq_num_index = 0;
  }

  //! Destructor
  ~TCP_Client_Application() { }

  //! ADD DOCUMENTATION HERE
  Slot<TCP_Client_Application, int> select;

  //! ADD DOCUMENTATION HERE
  void read_from_net(unsigned byte_size) {
    nbr_bytes_to_receive = byte_size;
    seq_num_val.set_size(10 + byte_size / 1460);
    seq_num_val.zeros();
    seq_num_time.set_size(10 + byte_size / 1460);
    seq_num_time.zeros();
    seq_num_val(0) = 0;
    seq_num_time(0) = 0;
    seq_num_index = 1;
  };

private:
  TCP_Receiver *tcp_receiver_p;
  TCP_Sender *tcp_sender_p;
  unsigned nbr_bytes_received;
  unsigned nbr_bytes_to_receive;

  vec seq_num_val;
  vec seq_num_time;
  int seq_num_index;

  void TCP_Client_Application::received_packet_indication(int label) {

    itpp::Packet &packet = tcp_receiver_p->get_user_message();
    nbr_bytes_received = nbr_bytes_received + packet.bit_size() / 8;
    delete &packet;

    if (seq_num_index >= seq_num_time.size()) {
      seq_num_time.set_size(2*seq_num_time.size(), true);
      seq_num_val.set_size(2*seq_num_val.size(), true);
    }

    seq_num_val(seq_num_index) = nbr_bytes_received;
    seq_num_time(seq_num_index) = Event_Queue::now();
    seq_num_index++;

    std::cout << "### sequence number: " << nbr_bytes_received
              << " ### time:" << Event_Queue::now() << std::endl;

    if (nbr_bytes_received >= nbr_bytes_to_receive) {
      std::cout << "###### Stop sender and receiver" << std::endl;
      tcp_receiver_p->release();
      tcp_sender_p->release();
      tcp_sender_p->save_trace("seq_num.it");
      seq_num_val.set_size(seq_num_index, true);
      seq_num_time.set_size(seq_num_index, true);
      save_to_file("seq_num.it");
    }
  }

  void TCP_Client_Application::save_to_file(string file) {

    it_file ff2(file);
    ff2 << Name("seq_num_val") << seq_num_val;
    ff2 << Name("seq_num_time") << seq_num_time;
    ff2 << Name("seq_num_index") << seq_num_index;
    ff2.flush();
    ff2.close();
  }

};

//@}

} // namespace itpp

#endif

#endif //TCP_CLIENT_SERVER_H
