/*!
 * \file
 * \brief Definition of Transport Control Protocol (TCP)
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
 *
 * Copyright (c) 2000-2004 IKR (formerly IND), University of Stuttgart
 * This file is part of the IKR (formerly IND) TCP Library.
 */

#ifndef TCP_H
#define TCP_H

#include <itpp/base/vec.h>

#if (defined(_MSC_VER) && defined(ITPP_SHARED_LIB) && !(defined(itpp_EXPORTS) || defined(itpp_debug_EXPORTS)))

#ifndef ITPP_PROTOCOL_EXCLUDED
#define ITPP_PROTOCOL_EXCLUDED
#pragma message( "PROTOCOL definitions are not available for MSVC shared builds" )
#endif

#else

#include <itpp/base/converters.h>
#include <itpp/protocol/packet.h>
#include <itpp/protocol/events.h>


namespace itpp
{

//! \addtogroup protocol
//@{

/*!
  TSequenceNumber represents a sequence number in a byte stream. As
  simulations may get quite long sequence numbers based on 32 bit integers
  may overflow. This is not a major problem as long as this overflow
  is considered when comparing sequence numbers. TSequenceNumber defines
  relational operators base on subtraction instead of using built-in
  int relational operators. This approach works well for arbitrary long
  simulations as long as the "real" sequence numbers (i.e. the numbers
  seen without overflow) compared to each other do not differ by more than
  2^31-1 (INT_MAX) which is not a serious restriction in a realistic
  TCP scenario.<p>
  @sa TTCPSegment
  @author Bodamer
*/
class Sequence_Number
{
public:
  //! Default constructor
  Sequence_Number() : seq(0) { }
  //! ADD DOCUMENTATION HERE
  Sequence_Number(const Sequence_Number &n) : seq(n.seq) { }
  //! ADD DOCUMENTATION HERE
  Sequence_Number &operator=(const Sequence_Number &n) { seq = n.seq; return *this; }
  //! ADD DOCUMENTATION HERE
  Sequence_Number &operator=(const int &rep) { seq = rep; return *this; }

  //relational operators
  //! ADD DOCUMENTATION HERE
  bool operator==(const Sequence_Number &n) const { return seq == n.seq; }
  //! ADD DOCUMENTATION HERE
  bool operator!=(const Sequence_Number &n) const { return seq != n.seq; }
  //! ADD DOCUMENTATION HERE
  bool operator>(const Sequence_Number &n) const { return (seq - n.seq) > 0; }
  //! ADD DOCUMENTATION HERE
  bool operator>=(const Sequence_Number &n) const { return (seq - n.seq) >= 0; }
  //! ADD DOCUMENTATION HERE
  bool operator<(const Sequence_Number &n) const { return (seq - n.seq) < 0; }
  //! ADD DOCUMENTATION HERE
  bool operator<=(const Sequence_Number &n) const { return (seq - n.seq) <= 0; }

  //addition and subtraction
  //! ADD DOCUMENTATION HERE
  Sequence_Number operator+(const int n) const { return Sequence_Number(seq + n); }
  //! ADD DOCUMENTATION HERE
  Sequence_Number &operator+=(const int n) { seq += n; return *this; }
  //! ADD DOCUMENTATION HERE
  Sequence_Number operator-(const int n) const { return Sequence_Number(seq - n); }
  //! ADD DOCUMENTATION HERE
  Sequence_Number &operator-=(const int n) { seq -= n; return *this; }
  //! ADD DOCUMENTATION HERE
  int operator-(const Sequence_Number &n) const { return seq - n.seq; }

  //! Access to internal representation
  int value() const { return seq; }

  //! ADD DOCUMENTATION HERE
  friend Sequence_Number operator+(const int n1, const Sequence_Number &n2) { return Sequence_Number(n1 + n2.seq); }
  //! ADD DOCUMENTATION HERE
  friend std::ostream &operator<<(std::ostream &os, const Sequence_Number &n) { os << n.seq; return os; }

protected:
  //! ADD DOCUMENTATION HERE
  Sequence_Number(int n) : seq(n) {}
  //! ADD DOCUMENTATION HERE
  int seq;
};

//! ADD DOCUMENTATION HERE
inline const Sequence_Number & min(const Sequence_Number &n1, const Sequence_Number &n2) { return (n1 < n2) ? n1 : n2; }
//! ADD DOCUMENTATION HERE
inline const Sequence_Number & max(const Sequence_Number &n1, const Sequence_Number &n2) { return (n1 > n2) ? n1 : n2; }



/*!
  TCP is a byte oriented protcol. Parts of the byte stream that is
  transmitted are called segments. They are identified by sequence numbers
  TCP_Segment contains fileds for the beginning sequence number and the
  sequence number of the first byte in the following segment (i.e. the
  sequence number of the last byte + 1). TCP_Segment provides several
  const methods to compare segments and to check whether they overlap.
  TCP_Segment is, e.g., used in TTCPPacket and TTCPReceiverBuffer.<p>
  @sa Sequence_Number
  @sa TTCPPacket
  @sa TTCPReceiverBuffer
  @author Lorang
  @author Bodamer
*/
class TCP_Segment
{
public:
  //! ADD DOCUMENTATION HERE
  TCP_Segment();
  //! ADD DOCUMENTATION HERE
  TCP_Segment(const Sequence_Number &sn_begin, const Sequence_Number &sn_end);
  //! ADD DOCUMENTATION HERE
  TCP_Segment(const TCP_Segment &segment);

  //modification
  //! ADD DOCUMENTATION HERE
  TCP_Segment &operator=(const TCP_Segment &segment);
  //! ADD DOCUMENTATION HERE
  void set_begin(const Sequence_Number &sn);
  //! ADD DOCUMENTATION HERE
  void set_end(const Sequence_Number &sn);
  //! ADD DOCUMENTATION HERE
  void combine(const TCP_Segment &segment);

  //query
  //! ADD DOCUMENTATION HERE
  bool operator==(const TCP_Segment &segment) const;
  //! ADD DOCUMENTATION HERE
  bool operator!=(const TCP_Segment &segment) const;
  //! ADD DOCUMENTATION HERE
  bool can_be_combined(const TCP_Segment &segment) const;
  //! ADD DOCUMENTATION HERE
  bool is_contained(const TCP_Segment &segment) const;
  //! ADD DOCUMENTATION HERE
  unsigned length() const;
  //! ADD DOCUMENTATION HERE
  Sequence_Number begin() const { return seq_begin; }
  //! ADD DOCUMENTATION HERE
  Sequence_Number end() const { return seq_end; }

  //! ADD DOCUMENTATION HERE
  friend std::ostream & operator<<(std::ostream &os, const TCP_Segment &segment);

protected:
  Sequence_Number  seq_begin; /**< \brief no. of first byte of segment */
  Sequence_Number  seq_end;   /**< \brief no. of last byte of segment + 1 */
};


/*!
  TCP_Packet is an IP packet with additional TCP header fields.
  Messages of this type are used for communication between TCP sender
  and receiver (data and ACKs).<p>
  The class contains a segment of type TCP_Segment with sequence numbers
  for begin and end of the segment. This is only used in data packets
  (i.e. from sender to receiver) while there are also fields containing
  sequence number for the next expected segment (used to ACK data
  packets, i.e. from receiver to sender) and for the advertised
  receiver window.<p>
  TCP_Packet has an additional field for the session id which my be compared
  with the pair (source port, destination port) in a real world
  TCP/IP packet.<p>
  Furthermore, some debug info containing state information for TCP sender
  or receiver may be attached to the message.<p>
  @sa TCP_Segment
  @sa TCP_Packet::TDebugInfo
  @author Grevent
  @author Lorang
  @author Bodamer
*/
class TCP_Packet : public itpp::Packet
{
public:
  //! ADD DOCUMENTATION HERE
  TCP_Packet();
  //! ADD DOCUMENTATION HERE
  TCP_Packet(const TCP_Packet &packet);
  //! ADD DOCUMENTATION HERE
  virtual ~TCP_Packet();
  //! ADD DOCUMENTATION HERE
  virtual TCP_Packet &clone() const;

  //TCP window mechanism
  //! ADD DOCUMENTATION HERE
  void set_segment(const TCP_Segment &seg) {  fSegment = seg; }
  //! ADD DOCUMENTATION HERE
  TCP_Segment get_segment() const { return fSegment; }
  //! ADD DOCUMENTATION HERE
  void set_wnd(unsigned val) { fWnd = val; }
  //! ADD DOCUMENTATION HERE
  unsigned get_wnd() const { return fWnd; }
  //! ADD DOCUMENTATION HERE
  void set_ACK(Sequence_Number val) { fACK = val; }
  //! ADD DOCUMENTATION HERE
  Sequence_Number get_ACK() const { return fACK; }

  //session control
  //! ADD DOCUMENTATION HERE
  void set_session_id(int val) { fSessionId = val; }
  //! ADD DOCUMENTATION HERE
  int get_session_id() const {  return fSessionId; }

  //debugging support
  //! ADD DOCUMENTATION HERE
  void set_destination_port(unsigned val) { fDestinationPort = val; }
  //! ADD DOCUMENTATION HERE
  unsigned get_destination_port() const { return fDestinationPort; }
  //! ADD DOCUMENTATION HERE
  void set_source_port(unsigned val) { fSourcePort = val; }
  //! ADD DOCUMENTATION HERE
  unsigned get_source_port() const { return fSourcePort; }
  //! ADD DOCUMENTATION HERE
  void set_info(unsigned ssThresh, unsigned recWnd, unsigned cWnd, double estRTT, Sequence_Number sndUna, Sequence_Number sndNxt, bool isRtx);
  //! ADD DOCUMENTATION HERE
  virtual void print_header(std::ostream &) const;

protected:
  //! ADD DOCUMENTATION HERE
  unsigned             fDestinationPort;
  //! ADD DOCUMENTATION HERE
  unsigned             fSourcePort;

  TCP_Segment   fSegment;     /**< \brief data segment to be transmitted */
  Sequence_Number   fACK;           /**< \brief acknowledgment (next expected sn) */
  unsigned         fWnd;           /**< \brief window size (advertised by receiver) */
  int    fSessionId;     /**< \brief session identifier */

  //Tracing

  //! ADD DOCUMENTATION HERE
  struct TDebugInfo {
    unsigned   fSSThresh; //!< ADD DOCUMENTATION HERE
    unsigned   fRecWnd; //!< ADD DOCUMENTATION HERE
    unsigned   fCWnd; //!< ADD DOCUMENTATION HERE
    double     fRTTEstimate; //!< ADD DOCUMENTATION HERE
    Sequence_Number   fSndUna; //!< ADD DOCUMENTATION HERE
    Sequence_Number   fSndNxt; //!< ADD DOCUMENTATION HERE
    bool   fRtxFlag; //!< ADD DOCUMENTATION HERE
  };

  //! ADD DOCUMENTATION HERE
  TDebugInfo  *fInfo;

  //! ADD DOCUMENTATION HERE
  friend std::ostream & operator<<(std::ostream &, TCP_Packet &);
};


/*!
  TTCPSender is an entity that models TCP flow and congestion control at
  the sender side. It is one of the key components of the TCP module.<p>
  TTCPSender communicates with its environment via three ports: <p>
  <ul>
  <li>"input": receive messages from a generator/application
  <li>"output": send TCP packets into the network
  <li>"ackinput": receive TCP ACK packets from the network originated by
  the TCP receiver
  </ul><p>
  The sender is activated when the user entity (e.g. a generator) offers
  a message, whose length indicates the amount of data to be transferred,
  or when an ACK is received from the network.
  In both cases data is only sent if there is any unsent data available,
  the sending window (determined by the congestion window and by the
  receiver advertised window) is large enough and the silly window
  syndrome avoidance algorithm is passed. Data that has been sent remains
  in the sender buffer, which is modelled in a virtual manner,
  until it is acknowledged. If data has been ACKed it is erased in the
  sender buffer and the TCP sender tries to get new data and send it.<p>
  TTCPSender can be used with different versions of congestion control:
  Tahoe, Reno, and New Reno. Moreover, various options (e.g. usage
  of Nagle/Karn/Go-Back-N algorithms) and parameters (e.g. mss, max cwnd,
  initial values, timer granularity) can be defined in the input file.<p>
  The class provides a simplified connection control via methods Setup and
  Release. Those methods are only used to reset internal state variables.
  No control messages are sent over the network, i.e. the sender assumes
  that the corresponding methods at the receiver side are called as well.<p>
  @sa TTCPSenderSet
  @sa TTCPReceiver
  @author Grevent
  @author Lorang
  @author Bodamer
*/
class TCP_Sender
{
public:
  //! ADD DOCUMENTATION HERE
  TCP_Sender(int label);

  //! ADD DOCUMENTATION HERE
  virtual ~TCP_Sender();

  //connection control
  //! ADD DOCUMENTATION HERE
  virtual void setup();
  //! ADD DOCUMENTATION HERE
  virtual void release(std::string trace_filename = "");

  //! Print support
  virtual void print_item(std::ostream &, const std::string &);

  //! ADD DOCUMENTATION HERE
  virtual void set_debug(const bool enable_debug = true);
  //! ADD DOCUMENTATION HERE
  virtual void set_debug(bool enable_debug, bool enable_signal_debug);
  //! ADD DOCUMENTATION HERE
  virtual void set_trace(const bool enable_trace = true);
  //! ADD DOCUMENTATION HERE
  virtual void save_trace(std::string filename);

  //! ADD DOCUMENTATION HERE
  Signal<itpp::Packet*> tcp_send;
  //! ADD DOCUMENTATION HERE
  Slot<TCP_Sender, itpp::Packet*> tcp_receive_ack;
  //! ADD DOCUMENTATION HERE
  Slot<TCP_Sender, itpp::Packet*> tcp_socket_write;
  //! ADD DOCUMENTATION HERE
  Slot<TCP_Sender, std::string> tcp_release;


  //Signal<itpp::Packet*> TcpSendSignal;
  //Slot<TCP_Sender, itpp::Packet*> TcpRecvAckSlot;
  //Slot<TCP_Sender, itpp::Packet*> SocketWriteSlot;
  //Slot<TCP_Sender, string> ReleaseSlot;

private:
  std::queue<itpp::Packet*> SocketWriteQueue;

  virtual void InitStatistics();       /**< \brief  reset statistic counters */
  virtual void HandleACK(TCP_Packet &);      /**< \brief  process incoming ACK */
  virtual void SendNewData(bool skipSWSA = false); /**< \brief  send new data */
  virtual void UnaRetransmit();       /**< \brief  TO or fast retransmit */
  virtual void FinishFastRecovery();       /**< \brief  actions at end of fast recovery */
  virtual void ReduceSSThresh();       /**< \brief  halving on dup ACK or TO */
  virtual void SendMsg(TCP_Packet & msg);    /**< \brief  access to network */
  virtual void HandleRtxTimeout(Ttype);       /**< \brief  what to do after Timeout */
  virtual void IdleCheck();       /**< \brief  check whether SSR after idle is done */
  virtual void HandleSWSATimeout(Ttype);      /**< \brief  handler for SWSA/Nagle timer */
  virtual unsigned GetNextSegmentSize(const Sequence_Number & begin);
  virtual unsigned SendWindow() const;   /**< \brief  min of CWnd, MaxCWnd and RecWnd */
  virtual double CalcRTOValue() const;   /**< \brief  value for rtx timer */
  virtual void SetRtxTimer();
  virtual void UpdateRTTVariables(double sampleRTT); /**< \brief  evaluate RTT measuremnt */
  virtual void TraceCWnd();
  virtual void TraceSentSeqNo(const Sequence_Number sn);
  virtual void TraceACKedSeqNo(const Sequence_Number sn);
  virtual void TraceRTTVariables(double sampleRTT);
  virtual void TraceSSThresh();
  virtual std::string GenerateFilename();

  void StopTransientPhase(); /**< \brief reset statistic counters */

  enum eTCPVersion {kTahoe, kReno, kNewReno};

  virtual void set_label(int label);

  //socket variables
  unsigned fLabel; // end point identification also used at receiver

  //message handling
  virtual void HandleUserMessageIndication(itpp::Packet *user_data);
  virtual void ReceiveMessageFromNet(itpp::Packet *msg);

  //parameters parsed from input file
  eTCPVersion  fTCPVersion;     // one of Reno, Tahoe, NewReno
  unsigned  fMSS;             // maximum segment size
  unsigned  fTCPIPHeaderLength; /**< \brief additional msg length (normally 40 bytes) */
  double       fInitialRTT;        /**< \brief RTT assumed before first measurement */
  unsigned  fInitialCWnd;       /**< \brief initial congestion window -> RFC 2581 */
  unsigned  fInitialSSThresh;   /**< \brief initial ssthresh -> RFC 2581 */
  unsigned   fMaxCWnd;           /**< \brief congestion window boundary */
  unsigned  fDupACKThreshold;   /**< \brief duplicate ACK threshold */
  double  fTimerGranularity;  /**< \brief granularity for rtx/idle timer values */
  double  fMaxRTO;     // max value of retransmission TO
  unsigned  fMaxBackoff;        /**< \brief max value of backoff value (Karn) */
  bool  fImmediateBackoffReset; /**< \brief reset backoff on first new ACK */
  bool  fKarn;              /**< \brief exclude rtx packets from RTTM */
  bool  fGoBackN;           /**< \brief use go-back-N on TO rtx */
  bool  fFlightSizeRecovery;// use flight size on fast rec. exit
  bool  fRenoConservation;  /**< \brief use cons. of packets on fast rec. (std.) */
  bool  fCarefulSSThreshReduction; /**< \brief use min(cwnd, flight size)/2 */
  bool  fIgnoreDupACKOnTORecovery; /**< \brief avoid fast rtx during TO recovery */
  bool  fCarefulMulFastRtxAvoidance; /**< \brief see RFC 2582, Section 5 */
  bool         fNagle;             /**< \brief use Nagle algorithm */
  double       fSWSATimerValue;    /**< \brief timer for silly wind. synd. avoidance */
  bool  fRestartAfterIdle;  /**< \brief perform SSR after idle period */
  bool  fTraceCWnd;     // print CWnd trace to cout
  bool  fTraceSentSeqNo;    /**< \brief print trace of sent SNs to cout */
  bool  fTraceACKedSeqNo;   /**< \brief print trace of received ACKs to cout */
  bool   fDebug;      // print additional information to cout
  bool   fTrace;      // store trace info in vectors

  //session identification
  int fSessionId;  /**< \brief is increased when Release is called */

  //TCP flow control (RFC 793)
  Sequence_Number fSndUna;       // lowest unacknowledged sn
  Sequence_Number fSndNxt;       // next byte to be sent
  Sequence_Number fSndMax;           /**< \brief highest byte that has been sent +1 */
  unsigned        fRecWnd;       // receiver advertised window
  unsigned        fMaxRecWnd;        /**< \brief maximum observed rec. window */
  Sequence_Number fUserNxt;       // next byte to be received from user

  //TCP congestion avoidance (RFC 2581, RFC 2001, RFC 2582)
  unsigned  fCWnd;       // congestion window
  unsigned  fSSThresh;      /**< \brief threshold between slow start and cong. avoid. */
  unsigned fDupACKCnt;     /**< \brief counter for duplicate ACKs */
  Sequence_Number  fRecoveryDupACK;  /**< \brief sndmax on 3rd dup ACK (see RFC 2582) */
  Sequence_Number  fRecoveryTO;  /**< \brief sndmax on TO (see RFC 2582) */

  //TCP timers
  TTimer<TCP_Sender>   fRtxTimer;     /**< \brief retransmission timer */
  Sequence_Number      fTimUna;       /**< \brief byte rtx timer is running for */
  TTimer<TCP_Sender>   fSWSATimer;    /**< \brief SWSA/Nagle timer */
  int           fBackoff;               /**< \brief backoff for rtx timer */
  bool          fPendingBackoffReset;  /**< \brief reset backoff on next new ACK */
  Ttype   fLastSendTime; /**< \brief idle detection for SSR */

  //round trip time measurement (RTTM)
  double          fSRTT;               /**< \brief smoothed mean round trip time */
  double          fRTTVar;             /**< \brief variance of RTT time */
  double          fRTTEstimate;        /**< \brief created from SRTT, RTTVar and gran. */
  Sequence_Number fRTTMByte;           /**< \brief byte for which RTTM is running */
  bool            fRTTMPending;        /**< \brief is currently RTTN running? */
  double          fRTTMStartTime;      /**< \brief beginning of a measurement */

  //statistic counters
  unsigned long       fNumberOfTimeouts;
  unsigned long       fNumberOfFastRetransmits;
  unsigned long       fNumberOfRTTMeasurements;
  unsigned long       fNumberOfReceivedACKs;
  unsigned long       fNumberOfIdleTimeouts;

  vec CWnd_val;
  vec CWnd_time;
  int CWnd_index;

  vec SSThresh_val;
  vec SSThresh_time;
  int SSThresh_index;

  ivec sent_seq_num_val;
  vec sent_seq_num_time;
  int sent_seq_num_index;

  ivec sender_recv_ack_seq_num_val;
  vec sender_recv_ack_seq_num_time;
  int sender_recv_ack_seq_num_index;

  vec RTTEstimate_val;
  vec RTTEstimate_time;
  int RTTEstimate_index;

  vec RTTsample_val;
  vec RTTsample_time;
  int RTTsample_index;
};


/*!
  TTCPReceiverBuffer is an important part of TTCPReceiver. It is
  much more complex than the buffer model at the sender side as it has to
  keep track of out of order segments. The segments received are combined if
  possible and the resulting non-contiguous segments are stored in a
  linked list. Moreover TTCPReceiverBuffer stores the smallest sequence
  number that has not (yet) been read out by the TCP receiver
  ("first byte").<p>
  Data is written to the buffer by TTCPReceiver using method Write and
  read out using method Read. Before reading data the TCP receiver has to
  check whether a data block in the stream is available by calling
  FirstBlockSize.<p>
  Furthermore, TTCPReceiver requires information from TTCPReceiverBuffer,
  e.g. about the next sequence number that is missing in the
  stream (NextExpected) or the receiver window that can be advertised
  to the sender (Window).<p>
  @sa TTCPReceiver
  @author Bodamer
  @author Kutter
*/
class TCP_Receiver_Buffer
{
public:
  //! ADD DOCUMENTATION HERE
  TCP_Receiver_Buffer();
  //! ADD DOCUMENTATION HERE
  TCP_Receiver_Buffer(const TCP_Receiver_Buffer &);
  //! ADD DOCUMENTATION HERE
  ~TCP_Receiver_Buffer();

  void reset(); /**< \brief clears internal list structure */

  void write(TCP_Segment newBlock);  /**< \brief add segment to the queue */
  void read(unsigned noOfBytes); /**< \brief read up to "noOfBytes" bytes from queue */

  unsigned first_block_size() const;        /**< \brief size of first complete block */
  Sequence_Number first_byte() const;      /**< \brief first byte stored or missing */
  Sequence_Number last_byte() const;       /**< \brief highest byte received (+1) */
  Sequence_Number next_expected() const;   /**< \brief first byte missing */
  //! ADD DOCUMENTATION HERE
  unsigned window() const;

  std::ostream &info(std::ostream &os, int detail = 0) const; /**< \brief print info */

protected:
  Sequence_Number fFirstByte;         /**< \brief first byte stored or missing */

  //! ADD DOCUMENTATION HERE
  std::list <TCP_Segment> fBufList;
};




/*!
  The TCP receiver models the receiver side of a TCP connection. It is
  connected to the network via ports "input" and "ackoutput" and to an
  entity modelling higher layers (e.g. a sink) via port "output".<p>
  Incoming TCP messages are used to update the receiver buffer and as a
  trigger to send an ACK message back to the TCP sender. If the "DelayedACK"
  option is used ACKs are only sent for every second packet (unless in the
  case of out of order packets).
  If the received packet is not out of order (i.e. the next expected
  sequence number in the byte stream has increased) new data is delivered
  to the higher layer. This delivery may be delayed in the receiver itself
  (user block processing delay) or by the subsequent entity if it blocks
  on incoming message indications.<p>
  A couple of parameters can be specified in the input file including delayed
  ACK timer value and timer granularity. The value of MSS should be equal
  to the one specified at the receiver side although it is only required
  for delayed ACK and receiver SWSA algorithms. <p>
  Like the TCP sender TTCPReceiver provides a simplified connection
  control via methods Setup and Release. Those methods are only used to
  reset internal state variables. No control messages are sent over
  the network, i.e. the receiver assumes that the corresponding methods
  at the sender side are called as well.<p>
  @sa TTCPReceiverSet
  @sa TTCPReceiverBuffer
  @sa TCP_Sender
  @author Grevent
  @author Lorang
  @author Bodamer
*/
class TCP_Receiver
{
public:

  //! ADD DOCUMENTATION HERE
  TCP_Receiver(int label);
  //! ADD DOCUMENTATION HERE
  virtual ~TCP_Receiver();

  //name  connection control
  //! ADD DOCUMENTATION HERE
  virtual void setup();
  //! ADD DOCUMENTATION HERE
  virtual void release(std::string trace_filename = "");

  //message handling
  itpp::Packet & get_user_message(); /**< \brief called by higher layer */
  bool is_user_message_available(); /**< \brief called by higher layer */

  //! ADD DOCUMENTATION HERE
  virtual void set_debug(const bool enable_debug = true);
  //! ADD DOCUMENTATION HERE
  virtual void set_debug(bool enable_debug, bool enable_signal_debug);
  //! ADD DOCUMENTATION HERE
  virtual void set_trace(const bool enable_trace = true);
  //! ADD DOCUMENTATION HERE
  virtual void save_trace(std::string filename);

  //! ADD DOCUMENTATION HERE
  Signal<itpp::Packet*> tcp_send_ack;
  //! ADD DOCUMENTATION HERE
  Slot<TCP_Receiver, itpp::Packet*> tcp_receive;
  Signal<int> tcp_new_data; /**< \brief indicate new data to higher layer */
  //! ADD DOCUMENTATION HERE
  Slot<TCP_Receiver, std::string> tcp_release;

private:
  void IndicateUserMessage(); /**< \brief indicate new data to higher layer */
  virtual void ReceiveMessageFromNet(itpp::Packet* msg); /**< \brief  receive from network */

  virtual void ReceiveDataPacket(TCP_Packet & packet); /**< \brief  receive TCP packet */
  virtual void SendACK(bool);        /**< \brief  send an ACK if necessary or enforced */
  virtual void ScheduleACKMessage(); /**< \brief  prepare ACK message for sending */
  virtual void SendACKMessage(Ttype); /**< \brief  called by ACK scheduling timer */
  virtual void DelayedACKHandler(Ttype);    /**< \brief  handler for delayed ACK timer */
  virtual void PeriodicACKHandler(Ttype);   /**< \brief  handler for periodic ACK timer */
  virtual void HandleEndOfProcessing(Ttype); /**< \brief  handler for user msg proc. */
  virtual void TraceReceivedSeqNo(const Sequence_Number &sn);
  virtual std::string GenerateFilename();

  //basic variables
  TCP_Receiver_Buffer  fReceiverBuffer;
  unsigned              fLabel;

  //parameters read by the parser
  unsigned     fTCPIPHeaderLength;  /**< \brief overhead by TCP/IP headers */
  unsigned     fMSS;            /**< \brief maximum segment size */
  unsigned  fBufferSize;       /**< \brief size of receiver buffer */
  bool         fDelayedACK;         /**< \brief if true, use delayed ACK */
  Ttype  fACKDelayTime;       /**< \brief maximum time for delaying ACKs */
  bool  fSendPeriodicACKs;   /**< \brief repeat ACKs if no data has been received */
  bool  fStrictPeriodicACKs; /**< \brief repeat ACKs independent of data received */
  Ttype  fPeriodicACKInterval;// interval after which an ACK is repeated
  Ttype  fACKSchedulingDelay; /**< \brief ACK delay due to proc. scheduling */
  bool  fACKOnBufferWrite;   /**< \brief send ACK after segment has been written */
  bool  fACKOnBufferRead;    /**< \brief send ACK when segment is read by user */
  unsigned  fMaxUserBlockSize;   /**< \brief max size of a user read block */
  unsigned  fMinUserBlockSize;   /**< \brief min size of a user read block */
  double  fUserBlockProcDelay; /**< \brief time to read data from buffer + forward */
  bool  fTrace;              /**< \brief trace end bytes of received packets */
  bool         fDebug;              /**< \brief print additional information */

  //specific TCP variables
  Sequence_Number fAdvRcvNxt;       /**< \brief advertised next expected byte (ACK) */
  unsigned        fAdvRcvWnd;       /**< \brief advertised receiver window */

  //session management
  int fSessionId;  /**< \brief increased by 1 on connection release */

  //ACK timers
  TTimer<TCP_Receiver> fDelayedACKTimer; /**< \brief timer for delayed ACK */
  TTimer<TCP_Receiver> fPeriodicACKTimer; /**< \brief timer for periodic ACKs */
  TTimer<TCP_Receiver> fACKSchedulingTimer;
  TCP_Packet *  fWaitingACKMsg;

  //user data delivery
  Packet * fUserMessage;
  TTimer<TCP_Receiver> fUserBlockProcTimer;


  //statistic counters
  ivec received_seq_num_val;
  vec received_seq_num_time;
  int received_seq_num_index;

};



// ------------------------------- Inline definitions ---------------------------------------------



inline Sequence_Number TCP_Receiver_Buffer::first_byte() const
{
  return fFirstByte;
}


inline Sequence_Number TCP_Receiver_Buffer::last_byte() const
{
  if (fBufList.empty()) {
    return fFirstByte;
  }
  else {
    return fBufList.back().end();
  }
}


inline Sequence_Number TCP_Receiver_Buffer::next_expected() const
{
  return fFirstByte + first_block_size();
}


inline void TCP_Segment::set_begin(const Sequence_Number &sn)
{
  seq_begin = sn;

  it_assert(seq_begin <= seq_end, "TCP_Segment::begin, end byte " + to_str(seq_end.value()) + " < begin byte " + to_str(seq_begin.value()));
}


inline void TCP_Segment::set_end(const Sequence_Number &sn)
{
  seq_end = sn;

  it_assert(seq_begin <= seq_end, "TCP_Segment::set_begin, end byte " + to_str(seq_end.value()) + " < begin byte " + to_str(seq_begin.value()));
}


inline bool TCP_Segment::operator==(const TCP_Segment &segment) const
{
  return (this->seq_begin == segment.seq_begin) && (this->seq_end == segment.seq_end);
}


inline bool TCP_Segment::operator!=(const TCP_Segment &segment) const
{
  return (this->seq_begin != segment.seq_begin) || (this->seq_end != segment.seq_end);
}


inline bool TCP_Segment::can_be_combined(const TCP_Segment &segment) const
{
  return (this->seq_begin <= segment.seq_end) && (segment.seq_begin <= this->seq_end);
}


inline bool TCP_Segment::is_contained(const TCP_Segment &segment) const
{
  return (segment.seq_begin <= this->seq_begin) && (this->seq_end <= segment.seq_end);
}


inline unsigned TCP_Segment::length() const
{
  return seq_end - seq_begin;
}

//@}


} // namespace itpp

#endif

#endif // #ifndef TCP_H

