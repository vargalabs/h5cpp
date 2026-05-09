/*!
 * \file
 * \brief Implementation of Transport Control Protocol (TCP)
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

#include <itpp/protocol/tcp.h>
#include <itpp/base/itfile.h>
#include <limits>
#include <cstdlib>
#include <ctime>

//! \cond

#ifdef _MSC_VER
#pragma warning(disable:4355)
#endif

namespace itpp
{

// -------------------- Default parameters ----------------------------------

// TCP sender and receiver

#define TCP_HEADERLENGTH        40

// TCP sender

#define TCP_VERSION kReno
#define TCP_SMSS 1460
#define TCP_INITIALCWNDREL 2     // related to MSS
#define TCP_INITIALSSTHRESHREL 1 // related to MaxCWnd
#define TCP_MAXCWNDREL 32        // related to MSS
#define TCP_DUPACKS 3
#define TCP_INITIALRTT 1
const double TCP_STIMERGRAN  = 0.2;
const double TCP_SWSATIMERVALUE = 0.2;
#define TCP_MAXBACKOFF 64
const double TCP_MAXRTO = std::numeric_limits<double>::max();
#define TCP_IMMEDIATEBACKOFFRESET false
#define TCP_TIMESTAMPS false
#define TCP_KARN true
#define TCP_NAGLE false
#define TCP_GOBACKN true
#define TCP_FLIGHTSIZERECOVERY false
#define TCP_RENOCONSERVATION true
#define TCP_CAREFULSSTHRESHREDUCTION true
#define TCP_IGNOREDUPACKONTORECOVERY true
#define TCP_CAREFULMULFASTRTXAVOIDANCE true
#define TCP_RESTARTAFTERIDLE true

// TCP receiver

#define TCP_RMSS 1460
const int TCP_BUFFERSIZE = std::numeric_limits<int>::max() / 4;
#define TCP_DELAYEDACK true
const double TCP_ACKDELAYTIME = 0.2;
#define TCP_SENDPERIODICACKS false
#define TCP_STRICTPERIODICACKS false
#define TCP_PERIODICACKINTERVAL 1
#define TCP_ACKSCHEDULINGDELAY 0
#define TCP_ACKBUFFERWRITE false
#define TCP_ACKBUFFERREAD true
const int TCP_MAXUSERBLOCKSIZE = std::numeric_limits<int>::max() / 4;
#define TCP_MINUSERBLOCKSIZE 1
#define TCP_USERBLOCKPROCDELAY 0

// TCP generator

#define TCPGEN_BLOCKSIZE 1460

// TCP applications

#define TCPAPP_MAXNOOFACTIVEAPPS 500
#define TCPAPP_DISTSTATARRAYSIZE 100
#define TCPAPP_DISTSTATMAXGOODPUT 1000
#define TCPAPP_DISTSTATMAXTRANSFERTIME 10000
#define TCPAPP_CONDMEANSTATARRAYSIZE 100
#define TCPAPP_CONDMEANSTATMAXREQSIZE 100000



inline int min(int opd1, int opd2)
{
  return (opd1 < opd2) ? opd1 : opd2;
}


inline int max(int opd1, int opd2)
{
  return (opd1 > opd2) ? opd1 : opd2;
}


// round is used to map a double value (e.g. RTO in TTCPSender) to the
// next higher value of a certain granularity (e.g. timer granularity).
inline double round(const double value, const double granularity)
{
  return (std::ceil(value / granularity) * granularity);
}

// -------------------- TCP_Segment ----------------------------------------

TCP_Segment::TCP_Segment() :
    seq_begin(),
    seq_end()
{
}

TCP_Segment::TCP_Segment(const Sequence_Number &sn_begin, const Sequence_Number &sn_end) :
    seq_begin(sn_begin),
    seq_end(sn_end)
{
  it_assert(seq_begin <= seq_end, "TCP_Segment::TCP_Segment, end byte " + to_str(seq_end.value()) +
            " < begin byte " + to_str(seq_begin.value()));
}


TCP_Segment::TCP_Segment(const TCP_Segment &segment) :
    seq_begin(segment.seq_begin),
    seq_end(segment.seq_end)
{
}


TCP_Segment &TCP_Segment::operator=(const TCP_Segment &segment)
{
  this->seq_begin = segment.seq_begin;
  this->seq_end   = segment.seq_end;

  return *this;
}


void TCP_Segment::combine(const TCP_Segment &segment)
{
  it_assert(can_be_combined(segment), "TCP_Segment::CombineWith, segments cannot be combined");

  seq_begin = min(seq_begin, segment.seq_begin);
  seq_end = max(seq_end, segment.seq_end);
}


std::ostream & operator<<(std::ostream &os, const TCP_Segment &segment)
{
  os << "(" << segment.seq_begin << "," << segment.seq_end << ")";
  return os;
}


// -------------------- TCP_Packet ----------------------------------------
TCP_Packet::TCP_Packet() :
    fSegment(),
    fACK(),
    fWnd(0),
    fSessionId(0),
    fInfo(0)
{
}


TCP_Packet::TCP_Packet(const TCP_Packet &packet) :
    fSegment(packet.fSegment),
    fACK(packet.fACK),
    fWnd(packet.fWnd),
    fSessionId(packet.fSessionId),
    fInfo(0)
{
  std::cout << "TCP_Packet::TCP_Packet ############" << " ";

  if (packet.fInfo != 0) {
    std::cout << "TCP_Packet::TCP_Packet rhs.fInfo ###########" << " ";
    fInfo = new TDebugInfo(*packet.fInfo);
  }
}


TCP_Packet::~TCP_Packet()
{
  delete fInfo;
}


TCP_Packet & TCP_Packet::clone() const
{
  return *new TCP_Packet(*this);
}


void TCP_Packet::set_info(unsigned ssThresh, unsigned recWnd, unsigned cWnd,
                          double estRTT, Sequence_Number sndUna,
                          Sequence_Number sndNxt, bool isRtx)
{
  if (fInfo == 0) {
    fInfo = new TDebugInfo;
  }

  fInfo->fSSThresh  = ssThresh;
  fInfo->fRecWnd       = recWnd;
  fInfo->fCWnd      = cWnd;
  fInfo->fRTTEstimate  = estRTT;
  fInfo->fSndUna    = sndUna;
  fInfo->fSndNxt    = sndNxt;
  fInfo->fRtxFlag    = isRtx;
}


void TCP_Packet::print_header(std::ostream &) const
{
  std::cout << "Hello!\n";

  std::cout << "Ses = " << get_session_id() << " ";

  std::cout << "Segment = " << get_segment() << " "
            << "ACK = "   << get_ACK()       << " "
            << "Wnd = "   << get_wnd()       << " ";

  std::cout << "DestPort = " << fDestinationPort << " "
            << "SourcePort = " << fSourcePort << " ";


  if (fInfo != 0) {
    std::cout << "SndSSThresh = " << fInfo->fSSThresh << "  ";
    std::cout << "RecWnd = "      << fInfo->fRecWnd   << "  ";
    std::cout << "SndCWnd = "     << fInfo->fCWnd     << "  ";
    std::cout << "RTTEstimate = " << fInfo->fRTTEstimate  << "  ";
    std::cout << "RtxFlag = "     << fInfo->fRtxFlag;
  }
  else
    std::cout << "fInfo = " << fInfo << " ";

  std::cout << std::endl;

}



std::ostream & operator<<(std::ostream & out, TCP_Packet & msg)
{
  msg.print_header(out);
  return out;
}


// -------------------- TCP_Sender ----------------------------------------
TCP_Sender::TCP_Sender(int label) :
    fLabel(label),
    fTCPVersion(TCP_VERSION),
    fMSS(TCP_SMSS),
    fTCPIPHeaderLength(TCP_HEADERLENGTH),
    fInitialRTT(TCP_INITIALRTT),
    fInitialCWnd(0),  // default initialization see below
    fInitialSSThresh(0),  // default initialization see below
    fMaxCWnd(0),   // default initialization see below
    fDupACKThreshold(TCP_DUPACKS),
    fTimerGranularity(TCP_STIMERGRAN),
    fMaxRTO(TCP_MAXRTO),
    fMaxBackoff(TCP_MAXBACKOFF),
    fImmediateBackoffReset(TCP_IMMEDIATEBACKOFFRESET),
    fKarn(TCP_KARN),
    fGoBackN(TCP_GOBACKN),
    fFlightSizeRecovery(TCP_FLIGHTSIZERECOVERY),
    fRenoConservation(TCP_RENOCONSERVATION),
    fCarefulSSThreshReduction(TCP_CAREFULSSTHRESHREDUCTION),
    fIgnoreDupACKOnTORecovery(TCP_IGNOREDUPACKONTORECOVERY),
    fCarefulMulFastRtxAvoidance(TCP_CAREFULMULFASTRTXAVOIDANCE),
    fNagle(TCP_NAGLE),
    fSWSATimerValue(TCP_SWSATIMERVALUE),
    fRestartAfterIdle(TCP_RESTARTAFTERIDLE),
    fDebug(false),
    fTrace(false),
    fSessionId(0),
    fRtxTimer(*this, &TCP_Sender::HandleRtxTimeout),
    fSWSATimer(*this, &TCP_Sender::HandleSWSATimeout)/*,*/
{

  // default values and parameter check for MaxCWND, InitCWND, InitSSThresh
  if (fMaxCWnd == 0) {
    fMaxCWnd = (unsigned)(TCP_MAXCWNDREL * fMSS);
  }
  else if (fMaxCWnd < fMSS) {
    //      throw (UL_CException("TCP_Sender::TCP_Sender",
    //                           "MaxCWnd must be >= MSS"));
  }

  if (fInitialCWnd == 0) {
    fInitialCWnd = (unsigned)(TCP_INITIALCWNDREL * fMSS);
  }
  else if ((fInitialCWnd < fMSS) || (fInitialCWnd > fMaxCWnd)) {
    //      throw (UL_CException("TCP_Sender::TCP_Sender",
    //                           "initial CWnd must be >= MSS and <= MaxCWnd"));
  }

  if ((fInitialSSThresh == 0) && (fMaxCWnd >= 2 * fMSS)) {
    fInitialSSThresh = (unsigned)(TCP_INITIALSSTHRESHREL * fMaxCWnd);
  }
  else if ((fInitialSSThresh < 2*fMSS) || (fInitialCWnd > fMaxCWnd)) {
    //      throw (UL_CException("TCP_Sender::TCP_Sender",
    //                           "initial CWnd must be >= 2*MSS and <= MaxCWnd"));
  }

  setup();

  InitStatistics();


  tcp_send.set_name("TCP Send");
  tcp_receive_ack.forward(this, &TCP_Sender::ReceiveMessageFromNet);
  tcp_receive_ack.set_name("TCP ACK");
  tcp_socket_write.forward(this, &TCP_Sender::HandleUserMessageIndication);
  tcp_socket_write.set_name("SocketWrite");
  tcp_release.forward(this, &TCP_Sender::release);
  tcp_release.set_name("Release");

}


TCP_Sender::~TCP_Sender()
{
}

void TCP_Sender::set_debug(const bool enable_debug)
{
  fDebug = enable_debug;
  tcp_send.set_debug(enable_debug);
}

void TCP_Sender::set_debug(bool enable_debug, bool enable_signal_debug)
{
  fDebug = enable_debug;
  tcp_send.set_debug(enable_signal_debug);
}

void TCP_Sender::set_trace(const bool enable_trace)
{
  fTrace = enable_trace;
}

void TCP_Sender::set_label(int label)
{
  fLabel = label;
}

void TCP_Sender::setup()
{
  fSndUna      = 0;
  fSndNxt      = 0;
  fSndMax      = 0;
  fMaxRecWnd   = 0;
  fRecWnd      = fMaxCWnd;
  fUserNxt     = 0;
  fCWnd      = fInitialCWnd;
  fSSThresh  = fInitialSSThresh;
  fRecoveryDupACK  = 0;
  fRecoveryTO  = 0;
  fDupACKCnt   = 0;

  // timers
  fBackoff = 1;
  fPendingBackoffReset = false;
  fLastSendTime = Event_Queue::now();

  // RTT measurement
  fTimUna  = 0;
  fSRTT        = 0;
  fRTTVar      = 0;
  fRTTEstimate = fInitialRTT;
  fRTTMPending = false;
  fRTTMByte  = 0;

  CWnd_val.set_size(1000);
  CWnd_val.zeros();
  CWnd_time.set_size(1000);
  CWnd_time.zeros();
  CWnd_val(0) = fInitialCWnd;
  CWnd_time(0) = 0;
  CWnd_index = 1;

  SSThresh_val.set_size(1000);
  SSThresh_val.zeros();
  SSThresh_time.set_size(1000);
  SSThresh_time.zeros();
  SSThresh_val(0) = fInitialSSThresh;
  SSThresh_time(0) = 0;
  SSThresh_index = 1;

  sent_seq_num_val.set_size(1000);
  sent_seq_num_val.zeros();
  sent_seq_num_time.set_size(1000);
  sent_seq_num_time.zeros();
  sent_seq_num_val(0) = 0;
  sent_seq_num_time(0) = 0;
  sent_seq_num_index = 1;

  sender_recv_ack_seq_num_val.set_size(1000);
  sender_recv_ack_seq_num_val.zeros();
  sender_recv_ack_seq_num_time.set_size(1000);
  sender_recv_ack_seq_num_time.zeros();
  sender_recv_ack_seq_num_val(0) = 0;
  sender_recv_ack_seq_num_time(0) = 0;
  sender_recv_ack_seq_num_index = 1;

  RTTEstimate_val.set_size(1000);
  RTTEstimate_val.zeros();
  RTTEstimate_time.set_size(1000);
  RTTEstimate_time.zeros();
  RTTEstimate_val(0) = fInitialRTT;
  RTTEstimate_time(0) = 0;
  RTTEstimate_index = 1;

  RTTsample_val.set_size(1000);
  RTTsample_val.zeros();
  RTTsample_time.set_size(1000);
  RTTsample_time.zeros();
  RTTsample_val(0) = 0;
  RTTsample_time(0) = 0;
  RTTsample_index = 1;

}

std::string TCP_Sender::GenerateFilename()
{
  time_t rawtime;
#ifndef _MSC_VER
  struct tm *timeinfo;
  timeinfo = localtime(&rawtime);
#else
  time(&rawtime);
  struct tm _timeinfo;
  struct tm *timeinfo = &_timeinfo;
  localtime_s(timeinfo, &rawtime);
#endif
  std::ostringstream filename_stream;
  filename_stream << "trace_tcp_sender_u" << fLabel
  << "_" << 1900 + timeinfo->tm_year
  << "_" << timeinfo->tm_mon
  << "_" << timeinfo->tm_mday
  << "__" << timeinfo->tm_hour
  << "_" << timeinfo->tm_min
  << "_" << timeinfo->tm_sec
  << "_.it";
  return filename_stream.str();
}


void TCP_Sender::release(std::string file)
{
  std::string filename;
  fSessionId++;

  fRtxTimer.Reset();
  fSWSATimer.Reset();

  if (fTrace) {
    if (file == "")
      filename = GenerateFilename();
    else
      filename = file;

    save_trace(filename);
  }
}


void TCP_Sender::InitStatistics()
{
  fNumberOfTimeouts         = 0;
  fNumberOfIdleTimeouts     = 0;
  fNumberOfFastRetransmits  = 0;
  fNumberOfRTTMeasurements  = 0;
  fNumberOfReceivedACKs     = 0;
}


void TCP_Sender::StopTransientPhase()
{
  InitStatistics();
}


void TCP_Sender::HandleUserMessageIndication(itpp::Packet *user_data_p)
{
  if (fDebug) {
    std::cout << "TCP_Sender::HandleUserMessageIndication"
              << " byte_size=" << user_data_p->bit_size() / 8
              << " ptr=" << user_data_p
              << " time=" << Event_Queue::now() << std::endl;
  }

  SocketWriteQueue.push(user_data_p);

  SendNewData();  // will call GetMessage (via GetNextSegmentSize)
  // if new data can be sent
}


void TCP_Sender::ReceiveMessageFromNet(itpp::Packet *msg)
{
  TCP_Packet & packet = (TCP_Packet &) * msg;

  if (fDebug) {
    std::cout << "TCP_Sender::ReceiveMessageFromNet"
              << " byte_size=" << msg->bit_size() / 8
              << " ptr=" << msg
              << " time=" << Event_Queue::now() << std::endl;
  }

  if ((packet.get_session_id() == fSessionId) && // ACK of current session
      (packet.get_ACK() >= fSndUna))  {    // ACK is OK
    HandleACK(packet);
  }

  delete &packet;
}


void TCP_Sender::HandleACK(TCP_Packet &msg)
{
  it_assert(msg.get_ACK() <= fSndMax, "TCP_Sender::HandleACK, received ACK > SndMax at ");

  fNumberOfReceivedACKs++;

  if (fTrace) {
    TraceACKedSeqNo(msg.get_ACK());
  }

  if (fDebug) {
    std::cout << "sender " << fLabel << ": "
              << "receive ACK: "
              << " t = " << Event_Queue::now() << ", "
              << msg << std::endl;
  }

  // update receiver advertised window size
  fRecWnd = msg.get_wnd();
  fMaxRecWnd = max(fRecWnd, fMaxRecWnd);

  if (msg.get_ACK() == fSndUna) {                  // duplicate ACK

    bool ignoreDupACK = (fSndMax == fSndUna); // no outstanding data

    if (fIgnoreDupACKOnTORecovery) {
      // don't count dupacks during TO recovery!
      if (fCarefulMulFastRtxAvoidance) {    // see RFC 2582, Section 5
        // like in Solaris
        ignoreDupACK = ignoreDupACK || (fSndUna <= fRecoveryTO);
      }
      else {
        // like in ns
        ignoreDupACK = ignoreDupACK || (fSndUna < fRecoveryTO);
      }
    }

    if (!ignoreDupACK) {
      fDupACKCnt++;   // count the number of duplicate ACKs

      if (fDupACKCnt == fDupACKThreshold) {
        // dupack threshold is reached
        fNumberOfFastRetransmits++;

        fRecoveryDupACK = fSndMax;

        ReduceSSThresh(); // halve ssthresh (in most cases)

        if ((fTCPVersion == kReno) || (fTCPVersion == kNewReno)) {
          fCWnd = fSSThresh;
        }
        else if (fTCPVersion == kTahoe) {
          fCWnd = fMSS;
        }

        if (fTCPVersion == kReno || fTCPVersion == kNewReno) {
          // conservation of packets:
          if (fRenoConservation) {
            fCWnd += fDupACKThreshold * fMSS;
          }
        }
        else if (fTCPVersion == kTahoe) {
          if (fGoBackN) {
            fSndNxt = fSndUna; // Go-Back-N (like in ns)
          }
        }

        UnaRetransmit();  // initiate retransmission
      }
      else if (fDupACKCnt > fDupACKThreshold) {
        if (fTCPVersion == kReno || fTCPVersion == kNewReno) {
          // conservation of packets
          // CWnd may exceed MaxCWnd during fast recovery,
          // however, the result of SendWindow() is always <= MaxCwnd
          if (fRenoConservation) {
            fCWnd += fMSS;
          }
        }
      }
    }
  }
  else {                                                 // new ACK
    Sequence_Number oldSndUna = fSndUna; // required for NewReno partial ACK
    fSndUna = msg.get_ACK();
    fSndNxt = max(fSndNxt, fSndUna);  // required in case of "Go-Back-N"

    // reset retransmission timer

    if ((fSndUna > fTimUna) && fRtxTimer.IsPending()) {
      // seq. no. for which rtx timer is running has been received
      fRtxTimer.Reset();
    }

    // backoff reset

    if (fImmediateBackoffReset) {
      fBackoff = 1;
    }
    else {
      if (fPendingBackoffReset) {
        fBackoff = 1;
        fPendingBackoffReset = false;
      }
      else if (fBackoff > 1) {
        // reset backoff counter only on next new ACK (this is probably
        // the way to operate intended by Karn)
        fPendingBackoffReset = true;
      }
    }

    // RTT measurement

    if ((fSndUna > fRTTMByte) && fRTTMPending) {
      UpdateRTTVariables(Event_Queue::now() - fRTTMStartTime);
      fRTTMPending = false;
    }

    // update CWnd and reset dupack counter

    if (fDupACKCnt >= fDupACKThreshold) {
      // we are in fast recovery
      if (fTCPVersion == kNewReno && fSndUna < fRecoveryDupACK) {
        // New Reno partial ACK handling
        if (fRenoConservation) {
          fCWnd = max(fMSS, fCWnd - (fSndUna - oldSndUna) + fMSS);
        }
        UnaRetransmit();  // start retransmit immediately
      }
      else {
        FinishFastRecovery();
      }
    }
    else {
      // no fast recovery
      fDupACKCnt = 0;
      if (fCWnd < fSSThresh) {
        // slow start phase
        fCWnd = min(fCWnd + fMSS, fMaxCWnd);
      }
      else {
        // congestion avoidance phase
        fCWnd += max(fMSS * fMSS / fCWnd, 1);   // RFC 2581
        fCWnd = min(fCWnd, fMaxCWnd);
      }
    }
  }  // new ACK

  SendNewData();  // try to send new data (even in the case that a retransmit
  // had to be performed)

  if (fTrace) {
    TraceCWnd();
  }
}


void TCP_Sender::SendNewData(bool skipSWSA)
{
  unsigned nextSegmentSize;

  it_assert(fSndUna <= fSndNxt, "TCP_Sender::SendNewData, SndUna > SndNxt in sender " +  to_str(fLabel) + "!");

  if (fRestartAfterIdle) {
    IdleCheck();
  }

  bool sillyWindowAvoidanceFailed = false;

  while (!sillyWindowAvoidanceFailed &&
         ((nextSegmentSize = GetNextSegmentSize(fSndNxt)) > 0)) {
    // there is new data to send and window is large enough

    // SWSA and Nagle (RFC 1122): assume PUSH to be set
    unsigned queuedUnsent = fUserNxt - fSndNxt;
    unsigned usableWindow = max(0, (fSndUna + SendWindow()) - fSndNxt);

    if (((unsigned)min(queuedUnsent, usableWindow) >= fMSS) ||
        ((!fNagle || (fSndUna == fSndNxt)) &&
         ((queuedUnsent <= usableWindow) ||  // Silly W. A.
          ((unsigned)min(queuedUnsent, usableWindow) >= fMaxRecWnd / 2)
         )
        ) ||
        skipSWSA
       ) {
      // Silly Window Syndrome Avoidance (SWSA) and Nagle passed

      TCP_Segment nextSegment(fSndNxt, fSndNxt + nextSegmentSize);
      TCP_Packet & msg = * new TCP_Packet();

      msg.set_segment(nextSegment);
      msg.set_session_id(fSessionId);
      msg.set_destination_port(fLabel); // The dest and src port are set to the same
      msg.set_source_port(fLabel);      // number for simplicity.
      msg.set_bit_size(8 * (nextSegmentSize + fTCPIPHeaderLength));

      if (fDebug) {
        std::cout << "TCP_Sender::SendNewData,"
                  << " nextSegmentSize=" << nextSegmentSize
                  << " fTCPIPHeaderLength=" << fTCPIPHeaderLength
                  << " byte_size=" << msg.bit_size() / 8
                  << " ptr=" << &msg
                  << " time=" << Event_Queue::now() << std::endl;
      }

      // no RTT measurement for retransmitted segments
      // changes on Dec. 13. 2002 (Ga, Bo, Scharf)

      if (!fRTTMPending && fSndNxt >= fSndMax) { // ##Bo##
        fRTTMStartTime = Event_Queue::now();
        fRTTMByte = nextSegment.begin();
        fRTTMPending = true;
      }

      fSndNxt += nextSegmentSize;
      fSndMax = max(fSndNxt, fSndMax);

      // reset SWSA timer if necessary
      if (skipSWSA) {
        skipSWSA = false;
      }
      else if (fSWSATimer.IsPending()) {
        fSWSATimer.Reset();
      }

      // set rtx timer if necessary
      if (!fRtxTimer.IsPending()) {
        SetRtxTimer();
      }


      if (fDebug) {
        msg.set_info(fSSThresh, fRecWnd, fCWnd, fRTTEstimate,
                     fSndUna, fSndNxt, false);
        std::cout << "sender " << fLabel
                  << ": send new data: "
                  << " t = " << Event_Queue::now() << ", "
                  << msg << std::endl;
      }

      SendMsg(msg);
    }
    else {
      sillyWindowAvoidanceFailed = true;
      // set SWSA timer
      if (!fSWSATimer.IsPending()) {
        fSWSATimer.Set(fSWSATimerValue);
      }
    }
  }

  // set timers in case that no new data could have been sent
  if (!fRtxTimer.IsPending()) {
    if (fSndMax > fSndUna) {  // there is outstanding data
      if (!fImmediateBackoffReset && fPendingBackoffReset) {
        // backoff is reset if no new data could have been sent since last
        // (successfull) retransmission; this is useful in case of
        // Reno recovery and multiple losses to avoid that in
        // the (unavoidable) series of timeouts the timer value
        // increases exponentially as this is not the intention
        // of the delayed backoff reset in Karn's algorithm
        fBackoff = 1;
        fPendingBackoffReset = false;
      }
      SetRtxTimer();
    }
  }
}


void TCP_Sender::UnaRetransmit()
{
  // resend after timeout or fast retransmit
  unsigned nextSegmentSize = GetNextSegmentSize(fSndUna);

  if (nextSegmentSize > 0) {
    TCP_Segment nextSegment(fSndUna, fSndUna + nextSegmentSize);
    TCP_Packet & msg = *new TCP_Packet();
    msg.set_segment(nextSegment);
    msg.set_session_id(fSessionId);
    msg.set_destination_port(fLabel); // The dest and src port are set to the same
    msg.set_source_port(fLabel);      // number for simplicity.
    msg.set_bit_size(8 * (nextSegmentSize + fTCPIPHeaderLength));

    fSndNxt = max(fSndNxt, fSndUna + nextSegmentSize);
    fSndMax = max(fSndNxt, fSndMax);

    // The RTT measurement is cancelled if the RTTM byte has a sequence
    // number higher or equal than the first retransmitted byte as
    // the ACK for the RTTM byte will be delayed by the rtx for at least
    // one round
    if (fKarn && (nextSegment.begin() <= fRTTMByte) && fRTTMPending) {
      fRTTMPending = false;
    }

    SetRtxTimer();

    if (fDebug) {
      msg.set_info(fSSThresh, fRecWnd, fCWnd, fRTTEstimate,
                   fSndUna, fSndNxt, true);
      std::cout << "sender " << fLabel;
      if (fDupACKCnt >= fDupACKThreshold) {
        std::cout << ": fast rtx: ";
      }
      else {
        std::cout << ": TO rtx: ";
      }
      std::cout << " t = " << Event_Queue::now() << ", "
                << msg << std::endl;
    }

    SendMsg(msg);
  }
  else {
    //      throw(UL_CException("TCP_Sender::UnaRetransmit", "no bytes to send"));
  }
}


void TCP_Sender::FinishFastRecovery()
{
  if (fTCPVersion == kTahoe) {
    fDupACKCnt = 0;
  }
  else if (fTCPVersion == kReno) {
    // Reno fast recovery
    fDupACKCnt = 0;
    if (fFlightSizeRecovery) {
      fCWnd = min(fSndMax - fSndUna + fMSS, fSSThresh);
    }
    else {
      fCWnd = fSSThresh;
    }
  }
  else if (fTCPVersion == kNewReno) {
    // New Reno fast recovery
    // "Set CWnd to ... min (ssthresh, FlightSize + MSS)
    // ... or ssthresh" (RFC 2582)
    if (fFlightSizeRecovery) {
      fCWnd = min(fSndMax - fSndUna + fMSS, fSSThresh);
    }
    else {
      fCWnd = fSSThresh;
    }
    fDupACKCnt = 0;
  }
}


void TCP_Sender::ReduceSSThresh()
{
  if (fCarefulSSThreshReduction) {
    // If Reno conservation is enabled the amount of
    // outstanding data ("flight size") might be rather large
    // and even larger than twice the old ssthresh value;
    // so this corresponds more to the ns behaviour where always cwnd is
    // taken instead of flight size.
    fSSThresh = max(2 * fMSS,
                    min(min(fCWnd, fSndMax - fSndUna), fRecWnd) / 2);
  }
  else {
    // use filght size / 2 as recommended in RFC 2581
    fSSThresh = max(2 * fMSS, min(fSndMax - fSndUna, fRecWnd) / 2);
  }

  it_assert(fSSThresh <= fMaxCWnd, "TCP_Sender::HandleACK, internal error: SndSSThresh is > MaxCWnd");

  if (fTrace) {
    TraceSSThresh();
  }
}


void TCP_Sender::SendMsg(TCP_Packet &msg)
{
  if (fTrace) {
    TraceSentSeqNo(msg.get_segment().end());
  }

  if (fRestartAfterIdle) {
    fLastSendTime = Event_Queue::now(); // needed for idle detection
  }

  tcp_send(&msg);
}


void TCP_Sender::IdleCheck()
{
  // idle detection according to Jacobson, SIGCOMM, 1988:
  // sender is currently idle and nothing has been send since RTO

  if (fSndMax == fSndUna && Event_Queue::now() - fLastSendTime > CalcRTOValue()) {
    fCWnd = fInitialCWnd; // see RFC2581

    fNumberOfIdleTimeouts++;

    if (fTrace) {
      TraceCWnd();
    }

    if (fDebug) {
      std::cout << "sender " << fLabel
                << ": idle timeout: "
                << "t = " << Event_Queue::now()
                << ", SndNxt = " << fSndNxt
                << ", SndUna = " << fSndUna
                << ", Backoff = " << fBackoff
                << std::endl;
    }
  }
}


void TCP_Sender::HandleRtxTimeout(Ttype)
{
  fNumberOfTimeouts++;

  // update backoff
  fBackoff = min(fMaxBackoff, fBackoff * 2);
  if (!fImmediateBackoffReset) {
    fPendingBackoffReset = false;
  }

  if (fDupACKCnt >= fDupACKThreshold) {
    FinishFastRecovery(); // reset dup ACK cnt and CWnd
  }
  else if (fDupACKCnt > 0) {
    fDupACKCnt = 0; // don't allow dupack action during TO recovery
  }

  // update CWnd and SSThresh
  ReduceSSThresh(); // halve ssthresh (in most cases)
  fCWnd = fMSS;               // not initial CWnd, see RFC 2581

  it_assert(fSSThresh <= fMaxCWnd, "TCP_Sender::HandleRtxTimeout, internal error: SndSSThresh is > MaxCWnd");

  fRecoveryTO = fSndMax;

  if (fGoBackN) {
    // go back N is mainly relevant in the case of multiple losses
    // which would lead to a series of timeouts without resetting sndnxt
    fSndNxt = fSndUna;
  }

  if (fDebug) {
    std::cout << "sender " << fLabel
              << ": rtx timeout: "
              << "t = " << Event_Queue::now()
              << ", SndNxt = " << fSndNxt
              << ", SndUna = " << fSndUna
              << std::endl;
  }

  if (fTrace) {
    TraceCWnd();
  }

  UnaRetransmit();    // initiate retransmission
}


void TCP_Sender::HandleSWSATimeout(Ttype)
{
  SendNewData(true);
}


unsigned TCP_Sender::GetNextSegmentSize(const Sequence_Number & begin)
{
  // try to get new user messages if available and necessary
  while ((fUserNxt < begin + fMSS) && (!SocketWriteQueue.empty())) {
    itpp::Packet *packet_p = SocketWriteQueue.front();
    SocketWriteQueue.pop();
    fUserNxt += (unsigned) packet_p->bit_size() / 8;
    delete packet_p;
  }

  Sequence_Number end = min(min(fUserNxt, begin + fMSS),
                            fSndUna + SendWindow());

  if (fDebug) {
    std::cout << "TCP_Sender::GetNextSegmentSize,"
              << " fUserNxt=" << fUserNxt
              << " begin_seq_num=" << begin
              << " fMSS=" << fMSS
              << " fSndUna=" << fSndUna
              << " SendWindow()=" << SendWindow()
              << " end_seq_num=" << end
              << " time=" << Event_Queue::now() << std::endl;
  }

  return max(0, end - begin);
}


unsigned TCP_Sender::SendWindow() const
{
  return min(fRecWnd, min(fMaxCWnd, fCWnd));
}


double TCP_Sender::CalcRTOValue() const
{
  static const double factor = 1 + 1e-8;
  // to avoid "simultaneous" TO/receive ACK events in case of const. RTT

  double rto = fBackoff * fRTTEstimate * factor;

  if (rto > fMaxRTO) {
    rto = fMaxRTO;
  }

  return rto;
}


void TCP_Sender::SetRtxTimer()
{
  double rto = CalcRTOValue();
  fRtxTimer.Set(rto);
  fTimUna = fSndUna;
  if (fDebug) {
    std::cout << "sender " << fLabel
              << ": set rtx timer: "
              << "t = " << Event_Queue::now()
              << ", RTO = " << rto
              << ", Backoff = " << fBackoff
              << ", TimUna = " << fTimUna
              << std::endl;
  }
}


void TCP_Sender::UpdateRTTVariables(double sampleRTT)
{
  if (fSRTT == 0) {
    fSRTT = sampleRTT;
    fRTTVar = sampleRTT / 2;
  }
  else {
    // see, e.g., Comer for the values used as weights
    fSRTT = 0.875 * fSRTT  + 0.125 * sampleRTT;
    fRTTVar = 0.75 * fRTTVar + 0.25 * fabs(sampleRTT - fSRTT);
  }

  fRTTEstimate = round(fSRTT + 4 * fRTTVar, fTimerGranularity);

  if (fTrace) {
    TraceRTTVariables(sampleRTT);
  }

  fNumberOfRTTMeasurements++;
}


void TCP_Sender::TraceRTTVariables(double sampleRTT)
{
  if (fDebug) {
    std::cout << "sender " << fLabel
              << ": RTT update: "
              << "t = " << Event_Queue::now()
              << ", sample = " << sampleRTT
              << ", SRTT = " << fSRTT
              << ", RTTVar = " << fRTTVar
              << ", RTTEstimate = " << fRTTEstimate
              << std::endl;
  }

  if (RTTsample_index >= RTTsample_time.size()) {
    RTTsample_time.set_size(2*RTTsample_time.size(), true);
    RTTsample_val.set_size(2*RTTsample_val.size(), true);
  }
  RTTsample_val(RTTsample_index) = sampleRTT;
  RTTsample_time(RTTsample_index) = Event_Queue::now();
  RTTsample_index++;

  if (RTTEstimate_index >= RTTEstimate_time.size()) {
    RTTEstimate_time.set_size(2*RTTEstimate_time.size(), true);
    RTTEstimate_val.set_size(2*RTTEstimate_val.size(), true);
  }
  RTTEstimate_val(RTTEstimate_index) = fRTTEstimate;
  RTTEstimate_time(RTTEstimate_index) = Event_Queue::now();
  RTTEstimate_index++;
}


void TCP_Sender::TraceCWnd()
{
  if (fDebug) {
    std::cout << "sender " << fLabel
              << " t = " << Event_Queue::now()
              << " cwnd = " << fCWnd << std::endl;
  }
  if (CWnd_index >= CWnd_time.size()) {
    CWnd_time.set_size(2*CWnd_time.size(), true);
    CWnd_val.set_size(2*CWnd_val.size(), true);
  }
  CWnd_val(CWnd_index) = fCWnd;
  CWnd_time(CWnd_index) = Event_Queue::now();
  CWnd_index++;

}

void TCP_Sender::TraceSSThresh()
{
  if (fDebug) {
    std::cout << "sender " << fLabel
              << " t = " << Event_Queue::now()
              << " cwnd = " << fSSThresh << std::endl;
  }
  if (SSThresh_index >= SSThresh_time.size()) {
    SSThresh_time.set_size(2*SSThresh_time.size(), true);
    SSThresh_val.set_size(2*SSThresh_val.size(), true);
  }
  SSThresh_val(SSThresh_index) = fSSThresh;
  SSThresh_time(SSThresh_index) = Event_Queue::now();
  SSThresh_index++;

}

void TCP_Sender::TraceSentSeqNo(const Sequence_Number sn)
{
  ////   UL_TEST_MAGIC;
  if (fDebug) {
    std::cout << "sender " << fLabel
              << " t = " << Event_Queue::now()
              << " sent = " << sn
              << std::endl;
  }
  if (sent_seq_num_index >= sent_seq_num_time.size()) {
    sent_seq_num_time.set_size(2*sent_seq_num_time.size(), true);
    sent_seq_num_val.set_size(2*sent_seq_num_val.size(), true);
  }
  sent_seq_num_val(sent_seq_num_index) = sn.value();
  sent_seq_num_time(sent_seq_num_index) = Event_Queue::now();
  sent_seq_num_index++;
}


void TCP_Sender::TraceACKedSeqNo(const Sequence_Number sn)
{
  if (fDebug) {
    std::cout << "sender " << fLabel
              << " t = " << Event_Queue::now()
              << " ACK = " << sn
              << std::endl;
  }

  if (sender_recv_ack_seq_num_index >= sender_recv_ack_seq_num_time.size()) {
    sender_recv_ack_seq_num_time.set_size(2*sender_recv_ack_seq_num_time.size(), true);
    sender_recv_ack_seq_num_val.set_size(2*sender_recv_ack_seq_num_val.size(), true);
  }
  sender_recv_ack_seq_num_val(sender_recv_ack_seq_num_index) = sn.value();
  sender_recv_ack_seq_num_time(sender_recv_ack_seq_num_index) = Event_Queue::now();
  sender_recv_ack_seq_num_index++;
}


void TCP_Sender::save_trace(std::string filename)
{

  CWnd_val.set_size(CWnd_index, true);
  CWnd_time.set_size(CWnd_index, true);

  SSThresh_val.set_size(SSThresh_index, true);
  SSThresh_time.set_size(SSThresh_index, true);

  sent_seq_num_val.set_size(sent_seq_num_index, true);
  sent_seq_num_time.set_size(sent_seq_num_index, true);

  sender_recv_ack_seq_num_val.set_size(sender_recv_ack_seq_num_index, true);
  sender_recv_ack_seq_num_time.set_size(sender_recv_ack_seq_num_index, true);

  RTTEstimate_val.set_size(RTTEstimate_index, true);
  RTTEstimate_time.set_size(RTTEstimate_index, true);

  RTTsample_val.set_size(RTTsample_index, true);
  RTTsample_time.set_size(RTTsample_index, true);

  if (fDebug) {
    std::cout << "CWnd_val" << CWnd_val << std::endl;
    std::cout << "CWnd_time" << CWnd_time << std::endl;
    std::cout << "CWnd_index" << CWnd_index << std::endl;

    std::cout << "SSThresh_val" << SSThresh_val << std::endl;
    std::cout << "SSThresh_time" << SSThresh_time << std::endl;
    std::cout << "SSThresh_index" << SSThresh_index << std::endl;

    std::cout << "sent_seq_num_val" << sent_seq_num_val << std::endl;
    std::cout << "sent_seq_num_time" << sent_seq_num_time << std::endl;
    std::cout << "sent_seq_num_index" << sent_seq_num_index << std::endl;

    std::cout << "sender_recv_ack_seq_num_val" << sender_recv_ack_seq_num_val << std::endl;
    std::cout << "sender_recv_ack_seq_num_time" << sender_recv_ack_seq_num_time << std::endl;
    std::cout << "sender_recv_ack_seq_num_index" << sender_recv_ack_seq_num_index << std::endl;

    std::cout << "RTTEstimate_val" << RTTEstimate_val << std::endl;
    std::cout << "RTTEstimate_time" << RTTEstimate_time << std::endl;
    std::cout << "RTTEstimate_index" << RTTEstimate_index << std::endl;

    std::cout << "RTTsample_val" << RTTsample_val << std::endl;
    std::cout << "RTTsample_time" << RTTsample_time << std::endl;
    std::cout << "RTTsample_index" << RTTsample_index << std::endl;

    std::cout << "TCP_Sender::saving to file: " << filename << std::endl;
  }

  it_file ff2;
  ff2.open(filename);

  ff2 << Name("CWnd_val") << CWnd_val;
  ff2 << Name("CWnd_time") << CWnd_time;
  ff2 << Name("CWnd_index") << CWnd_index;

  ff2 << Name("SSThresh_val") << SSThresh_val;
  ff2 << Name("SSThresh_time") << SSThresh_time;
  ff2 << Name("SSThresh_index") << SSThresh_index;

  ff2 << Name("sent_seq_num_val") << sent_seq_num_val;
  ff2 << Name("sent_seq_num_time") << sent_seq_num_time;
  ff2 << Name("sent_seq_num_index") << sent_seq_num_index;

  ff2 << Name("sender_recv_ack_seq_num_val") << sender_recv_ack_seq_num_val;
  ff2 << Name("sender_recv_ack_seq_num_time") << sender_recv_ack_seq_num_time;
  ff2 << Name("sender_recv_ack_seq_num_index") << sender_recv_ack_seq_num_index;

  ff2 << Name("RTTEstimate_val") << RTTEstimate_val;
  ff2 << Name("RTTEstimate_time") << RTTEstimate_time;
  ff2 << Name("RTTEstimate_index") << RTTEstimate_index;

  ff2 << Name("RTTsample_val") << RTTsample_val;
  ff2 << Name("RTTsample_time") << RTTsample_time;
  ff2 << Name("RTTsample_index") << RTTsample_index;

  ff2.flush();
  ff2.close();
}


void TCP_Sender::print_item(std::ostream &, const std::string & keyword)
{
  if (keyword == "Label") {
    std::cout << fLabel;
  }
  else if (keyword == "CWnd") {
    std::cout << fCWnd;
  }
  else if (keyword == "SSThresh") {
    std::cout << fSSThresh;
  }
  else if (keyword == "SRTT") {
    std::cout << fSRTT;
  }
  else if (keyword == "RTTvar") {
    std::cout << fRTTVar;
  }
  else if (keyword == "Backoff") {
    std::cout << fBackoff;
  }
  else if (keyword == "RTO") {
    std::cout << CalcRTOValue();
  }
  else if (keyword == "NoOfFastRets") {
    std::cout << fNumberOfFastRetransmits;
  }
  else if (keyword == "NoOfRetTOs") {
    std::cout << fNumberOfTimeouts;
  }
  else if (keyword == "NoOfIdleTOs") {
    std::cout << fNumberOfIdleTimeouts;
  }
  else if (keyword == "NoOfRTTMs") {
    std::cout << fNumberOfRTTMeasurements;
  }
  else if (keyword == "NoOfRecACKs") {
    std::cout << fNumberOfReceivedACKs;
  }
  else {
  }
}


// -------------------- TCP_Receiver_Buffer ----------------------------------------
TCP_Receiver_Buffer::TCP_Receiver_Buffer() :
    fFirstByte()
{
}


TCP_Receiver_Buffer::TCP_Receiver_Buffer(const TCP_Receiver_Buffer &  rhs) :
    fFirstByte(rhs.fFirstByte),
    fBufList(rhs.fBufList)
{
}


void TCP_Receiver_Buffer::reset()
{
  fBufList.clear();
  fFirstByte = 0;
}


TCP_Receiver_Buffer::~TCP_Receiver_Buffer()
{
}


void TCP_Receiver_Buffer::write(TCP_Segment newBlock)
{
  // error cases
  it_assert(newBlock.begin() <= newBlock.end(), "TCP_Receiver_Buffer::Write, no valid segment");

  // cut blocks beginning before fFirstByte
  if (newBlock.begin() < fFirstByte) {
    if (newBlock.end() > fFirstByte) {
      newBlock.set_begin(fFirstByte);
    }
    else {
      return; //// TODO: Is this strange?
    }
  }

  if (newBlock.length() == 0) { // empty block, nothing to do
    return;
  }

  if (fBufList.empty() || (newBlock.begin() > fBufList.back().end())) {
    // new block is behind last block in buffer
    fBufList.push_back(newBlock);
  }
  else {
    // skip list entries if beginning of newBlock > end of current one
    // (search for correct list position)
    std::list<TCP_Segment>::iterator iter;
    iter = fBufList.begin();
    while (newBlock.begin() > iter->end()) {
      iter++;
      it_assert(iter != fBufList.end(), "TCP_Receiver_Buffer::Write, internal error");
    }

    TCP_Segment & exBlock = *iter;

    if (exBlock.can_be_combined(newBlock)) {
      // overlapping or contiguous blocks -> combine
      exBlock.combine(newBlock);

      // check following blocks
      iter++;
      while ((iter != fBufList.end()) &&
             exBlock.can_be_combined(*iter)) {
        exBlock.combine(*iter);
        iter = fBufList.erase(iter);
      }
    }
    else {
      // no overlap, newBlock lies between two existing list entries
      // new list entry has to be created

      fBufList.insert(iter, newBlock);
    }
  }

  it_assert(!fBufList.empty() && fBufList.front().begin() >= fFirstByte, "TCP_Receiver_Buffer::Write, internal error");

}


// The amount of data read from the buffer is given as parameter. It has
// to be less than or equal to the size of the first block stored. This
// mean the caller of Read should first check how much data is available
// by calling FirstBlockSize.
void TCP_Receiver_Buffer::read(unsigned noOfBytes)
{
  it_assert(first_block_size() > 0, "TCP_Receiver_Buffer::Read,  No block to read");
  it_assert(noOfBytes <= first_block_size(), "TCP_Receiver_Buffer::Read, submitted block size not valid");


  if (noOfBytes < first_block_size()) {
    fBufList.front().set_begin(fBufList.front().begin() + noOfBytes);
  }
  else { // first block will be read completely
    fBufList.pop_front();
  }
  fFirstByte += noOfBytes;

  it_assert(fBufList.empty() || fBufList.front().begin() >= fFirstByte, "TCP_Receiver_Buffer::Read, internal error");
}


// FirstBlockSize returns the size of the first block stored in the
// buffer or 0 if the buffer is empty
unsigned TCP_Receiver_Buffer::first_block_size() const
{
  if (!fBufList.empty() && (fBufList.front().begin() == fFirstByte)) {
    return fBufList.front().length();
  }
  else {
    return 0;
  }
}


std::ostream & TCP_Receiver_Buffer::info(std::ostream &os, int detail) const
{
  os << "receiver buffer information" << std::endl
  << "number of blocks: " << fBufList.size() << std::endl
  << "first byte stored: " << fFirstByte << std::endl
  << "last byte stored +1: " << last_byte() << std::endl
  << "next byte expected: " << next_expected() << std::endl;

  if (detail > 0) {
    os << "segments in receiver buffer:" << std::endl;

    typedef std::list<TCP_Segment>::const_iterator LI;
    for (LI i = fBufList.begin(); i != fBufList.end(); ++i) {
      const TCP_Segment & block = *i;
      os << ". segment: " << block << std::endl;
    }

  }

  return os;
}


// -------------------- TCP_Receiver ----------------------------------------
TCP_Receiver::TCP_Receiver(int label) :
    fReceiverBuffer(),
    fLabel(label),
    fTCPIPHeaderLength(TCP_HEADERLENGTH),
    fMSS(TCP_RMSS),
    fBufferSize(TCP_BUFFERSIZE),
    fDelayedACK(TCP_DELAYEDACK),
    fACKDelayTime(TCP_ACKDELAYTIME),
    fSendPeriodicACKs(TCP_SENDPERIODICACKS),
    fStrictPeriodicACKs(TCP_STRICTPERIODICACKS),
    fPeriodicACKInterval(TCP_PERIODICACKINTERVAL),
    fACKSchedulingDelay(TCP_ACKSCHEDULINGDELAY),
    fACKOnBufferWrite(TCP_ACKBUFFERWRITE),
    fACKOnBufferRead(TCP_ACKBUFFERREAD),
    fMaxUserBlockSize(TCP_MAXUSERBLOCKSIZE),
    fMinUserBlockSize(TCP_MINUSERBLOCKSIZE),
    fUserBlockProcDelay(TCP_USERBLOCKPROCDELAY),
    fTrace(false),
    fDebug(false),
    fSessionId(0),
    fDelayedACKTimer(*this, &TCP_Receiver::DelayedACKHandler),
    fPeriodicACKTimer(*this, &TCP_Receiver::PeriodicACKHandler),
    fACKSchedulingTimer(*this, &TCP_Receiver::SendACKMessage),
    fWaitingACKMsg(0),
    fUserBlockProcTimer(*this, &TCP_Receiver::HandleEndOfProcessing)
{
  fUserMessage = NULL;


  if (!fACKOnBufferRead && !fACKOnBufferWrite) {
    //     throw(UL_CException("TCP_Receiver::TCP_Receiver",
    //                          "ACKs must be sent on buffer read or write or both"));
  }

  setup();

  tcp_receive.forward(this, &TCP_Receiver::ReceiveMessageFromNet);
  tcp_receive.set_name("TCP Receive");
  tcp_send_ack.set_name("TCP send ACK");
  tcp_new_data.set_name("TCP New Data");
  tcp_release.forward(this, &TCP_Receiver::release);
  tcp_release.set_name("TCP Release");

}


TCP_Receiver::~TCP_Receiver()
{
  delete fWaitingACKMsg;
  delete fUserMessage;
}


void TCP_Receiver::set_debug(const bool enable_debug)
{
  fDebug = enable_debug;
  tcp_send_ack.set_debug(enable_debug);
  tcp_new_data.set_debug();
}

void TCP_Receiver::set_debug(bool enable_debug, bool enable_signal_debug)
{
  fDebug = enable_debug;
  tcp_send_ack.set_debug(enable_signal_debug);
  tcp_new_data.set_debug();
}

void TCP_Receiver::set_trace(const bool enable_trace)
{
  fTrace = enable_trace;
}



void TCP_Receiver::setup()
{
  fAdvRcvWnd = 0;
  fAdvRcvNxt = 0;

  if (fSendPeriodicACKs) {
    fPeriodicACKTimer.Set(fPeriodicACKInterval);
  }

  fReceiverBuffer.reset();

  received_seq_num_val.set_size(1000);
  received_seq_num_val.zeros();
  received_seq_num_time.set_size(1000);
  received_seq_num_time.zeros();
  received_seq_num_val(0) = 0;
  received_seq_num_time(0) = 0;
  received_seq_num_index = 1;
}

std::string TCP_Receiver::GenerateFilename()
{
  time_t rawtime;
#ifndef _MSC_VER
  struct tm *timeinfo;
  timeinfo = localtime(&rawtime);
#else
  time(&rawtime);
  struct tm _timeinfo;
  struct tm *timeinfo = &_timeinfo;
  localtime_s(timeinfo, &rawtime);
#endif
  std::ostringstream filename_stream;
  filename_stream << "trace_tcp_receiver_u" << fLabel
  << "_" << 1900 + timeinfo->tm_year
  << "_" << timeinfo->tm_mon
  << "_" << timeinfo->tm_mday
  << "__" << timeinfo->tm_hour
  << "_" << timeinfo->tm_min
  << "_" << timeinfo->tm_sec
  << "_.it";
  return filename_stream.str();
}

void TCP_Receiver::release(std::string file)
{
  std::string filename;
  fSessionId++;

  if (fWaitingACKMsg != 0) {
    delete fWaitingACKMsg;
    fWaitingACKMsg = 0;
  }
  if (fUserMessage != 0) {
    delete fUserMessage;
    fUserMessage = 0;
  }

  fUserBlockProcTimer.Reset();
  fDelayedACKTimer.Reset();
  fPeriodicACKTimer.Reset();
  fACKSchedulingTimer.Reset();

  if (fTrace) {
    if (file == "")
      filename = GenerateFilename();
    else
      filename = file;

    save_trace(filename);
  }
}


void TCP_Receiver::ReceiveMessageFromNet(itpp::Packet *msg)
{
  TCP_Packet & packet = (TCP_Packet &) * msg;
  if (packet.get_destination_port() == fLabel) {
    if (packet.get_session_id() == fSessionId) {
      ReceiveDataPacket(packet);
    }
    else {
      it_warning("Received a TCP packet with wrong SessionId");
      std::cout << "TCP_Receiver::ReceiveMessageFromNet, "
                << "fLabel= " << fLabel
                << "fSessionId= " << fSessionId << std::endl;
      std::cout << "packet=" << packet
                << ", next exp. = " << fReceiverBuffer.next_expected()
                << std::endl;
      exit(0);
    }
  }
  else {
    it_warning("Received a TCP packet with label");
    exit(0);
  }
}


void TCP_Receiver::ReceiveDataPacket(TCP_Packet &msg)
{
  TCP_Segment segment = msg.get_segment();

  bool isOutOfOrder = (segment.begin() > fReceiverBuffer.next_expected()) ||
                      (segment.end() <= fReceiverBuffer.next_expected());

  if (fDebug) {
    std::cout << "TCP_Receiver::ReceiveDataPacket receiver: " << fLabel << ": "
              << "receive msg: "
              << "t = " << Event_Queue::now()
              << ", next exp. = " << fReceiverBuffer.next_expected()
              << ", " << msg << std::endl;
  }

  if (fTrace) {
    TraceReceivedSeqNo(segment.end());
  }

  it_assert(segment.end() <= fReceiverBuffer.first_byte() + fBufferSize, "TCP_Receiver::ReceiveTCPPacket, packet exceeds window at ");
  it_assert(segment.begin() < segment.end(), "TCP_Receiver::ReceiveTCPPacket, silly packet received at ");

  fReceiverBuffer.write(segment);

  if (isOutOfOrder) {
    SendACK(true);                    // create dupack conditionless
  }
  else {
    if (fACKOnBufferWrite) {
      SendACK(false);
    }
    IndicateUserMessage();
  }

  delete &msg;
}


void TCP_Receiver::IndicateUserMessage()
{
  if (fUserMessage == 0) {
    // receive a block
    unsigned noOfBytes = min(fReceiverBuffer.first_block_size(),
                             fMaxUserBlockSize);

    if (fDebug) {
      std::cout << "TCP_Receiver::IndicateUserMessage  "
                << "t = " << Event_Queue::now()
                << " noOfBytes = " << noOfBytes
                << " firstBlock = " << fReceiverBuffer.first_block_size()
                << std::endl;
    }

    if (noOfBytes >= fMinUserBlockSize) {
      fUserMessage = new Packet();
      fUserMessage->set_bit_size(8*noOfBytes);
      fUserBlockProcTimer.Set(fUserBlockProcDelay);
    }
  }
}


bool TCP_Receiver::is_user_message_available()
{
  if (fUserMessage != 0) {
    return true;
  }

  unsigned noOfBytes = min(fReceiverBuffer.first_block_size(),
                           fMaxUserBlockSize);

  if (noOfBytes >= fMinUserBlockSize) {
    fUserMessage = new Packet();
    fUserMessage->set_bit_size(8*noOfBytes);
    return true;
  }
  else {
    return false;
  }
}


itpp::Packet & TCP_Receiver::get_user_message()
{
  it_assert(fUserMessage != 0, "TCP_Receiver::GetUserMessage, no message available");
  if (fDebug) {
    std::cout << "TCP_Receiver::GetUserMessage  "
              << "receiver: " << fLabel << ": "
              << "read from buffer: "
              << "t = " << Event_Queue::now()
              << ", user msg length = " << (fUserMessage->bit_size() / 8)
              << ", first byte = " << fReceiverBuffer.first_byte()
              << ", first block size = " << fReceiverBuffer.first_block_size()
              << std::endl;
  }

  fReceiverBuffer.read(fUserMessage->bit_size() / 8);
  if (fACKOnBufferRead) {
    SendACK(false);  // send acknowledgement
  }

  itpp::Packet & msg = *fUserMessage;
  fUserMessage = 0;

  if (fReceiverBuffer.first_block_size() > 0) {
    IndicateUserMessage();
  }

  return msg;
}



void TCP_Receiver::HandleEndOfProcessing(Ttype)
{
  it_assert(fUserMessage != 0, "TCP_Receiver::HandleEndOfProcessing, no message available");


  tcp_new_data(fLabel);
}


void TCP_Receiver::DelayedACKHandler(Ttype)
{
  if (fDebug) {
    std::cout << "TCP_Receiver::DelayedACKHandler  "
              << "receiver " << fLabel
              << ": delACK TO: "
              << "t = " << Event_Queue::now() << std::endl;
  }

  SendACK(true);
}


void TCP_Receiver::PeriodicACKHandler(Ttype)
{
  if (fDebug) {
    std::cout << "TCP_Receiver::PeriodicACKHandler"
              << "receiver " << fLabel
              << ": periodicACK TO: "
              << "t = " << Event_Queue::now() << std::endl;
  }

  SendACK(true);
}


void TCP_Receiver::SendACK(bool sendConditionless)
{
  // sendConditionless is set
  // ... if packet was received out of order or
  // ... if delayed ACK timer has expired

  // Bei eingeschaltetem "delayed ACK" wird ein ACK nur
  // gesendet, wenn das Fenster um 2MSS oder 35% der
  // maximalen Fenstergroesse verschoben worden ist
  // ... oder nach delayed ACK Timeout
  // ... oder wenn es das ACK fur ein Out of Order Segment ist
  // ... oder (in der Realitat), wenn ich auch was zu senden habe.

  if (sendConditionless || !fDelayedACK ||
      (fReceiverBuffer.next_expected() - fAdvRcvNxt >= (int)(2 * fMSS)) ||
      (fReceiverBuffer.next_expected() - fAdvRcvNxt >=
       (int)(0.35 * fBufferSize))) {
    // Remark: RFC2581 recommends to acknowledge every second
    // packet conditionless (without setting this as a requirement)
    // in order to avoid excessive ack delays when the receiver MSS
    // is larger than the sender MSS. In this uni-directional
    // implementation, the receiver's MSS is not actively
    // used for sending but only for deciding when acknowledgments
    // have to be returned. Thus, the best solution to account for
    // RFC2581 is to set the receiver's MSS always equal to the
    // sender's MSS.

    // Receiver Silly Window Syndrome Avoidance:

    if (fAdvRcvNxt + fAdvRcvWnd + min(fBufferSize / 2, fMSS)
        <= fReceiverBuffer.first_byte() + fBufferSize) {
      // Die rechte Grenze des Empfangerfensters wird nur anders angezeigt
      // als beim letzten ACK, wenn sie sich seither um mindestens
      // min (BufferSize/ 2, MSS) geandert hat.
      fAdvRcvWnd = fBufferSize - fReceiverBuffer.first_block_size();
    }
    else {
      fAdvRcvWnd = fAdvRcvNxt + fAdvRcvWnd - fReceiverBuffer.next_expected();
    }

    fAdvRcvNxt = fReceiverBuffer.next_expected();

    if (fSendPeriodicACKs &&
        (!fStrictPeriodicACKs || !fPeriodicACKTimer.IsPending())) {
      fPeriodicACKTimer.Set(fPeriodicACKInterval);
    }

    if (fDelayedACK && fDelayedACKTimer.IsPending()) {
      fDelayedACKTimer.Reset();
    }

    ScheduleACKMessage();
  }
  else {
    if (!fDelayedACKTimer.IsPending()) {
      fDelayedACKTimer.Set(fACKDelayTime);
      if (fDebug) {
        std::cout << "TCP_Receiver::SendACK"
                  << "receiver " << fLabel
                  << ": set delACK timer: "
                  << "t = " << Event_Queue::now() << std::endl;
      }
    }
  }
}


void TCP_Receiver::ScheduleACKMessage()
{
  if (fWaitingACKMsg == 0) {
    fWaitingACKMsg = new TCP_Packet;
  }

  fWaitingACKMsg->set_ACK(fAdvRcvNxt);
  fWaitingACKMsg->set_wnd(fAdvRcvWnd);
  fWaitingACKMsg->set_session_id(fSessionId);
  fWaitingACKMsg->set_destination_port(fLabel);
  fWaitingACKMsg->set_source_port(fLabel);
  fWaitingACKMsg->set_bit_size(8*fTCPIPHeaderLength);

  if (fACKSchedulingDelay > 0) {
    if (!fACKSchedulingTimer.IsPending()) {
      fACKSchedulingTimer.Set(fACKSchedulingDelay);
    }
  }
  else {
    SendACKMessage(Event_Queue::now());
  }
}


void TCP_Receiver::SendACKMessage(Ttype)
{
  it_assert(fWaitingACKMsg != 0, "TCP_Receiver::SendACKMessage, no ACK message waiting");

  if (fDebug) {
    std::cout << "TCP_Receiver::SendACKMessage Ack sent"
              << "receiver " << fLabel
              << ": send ACK: "
              << "t = " << Event_Queue::now()
              << ", " << (*fWaitingACKMsg)
              << " byte_size=" << fWaitingACKMsg->bit_size() / 8
              << " ptr=" << fWaitingACKMsg << std::endl;
  }

  tcp_send_ack(fWaitingACKMsg);

  fWaitingACKMsg = 0;
}


void TCP_Receiver::TraceReceivedSeqNo(const Sequence_Number &sn)
{
  if (fDebug) {
    std::cout << "TCP_Receiver::TraceReceivedSeqNo  "
              << "receiver " << fLabel
              << " t = " << Event_Queue::now()
              << " sn = " << sn << std::endl;
  }
  if (received_seq_num_index >= received_seq_num_time.size()) {
    received_seq_num_time.set_size(2*received_seq_num_time.size(), true);
    received_seq_num_val.set_size(2*received_seq_num_val.size(), true);
  }
  received_seq_num_val(received_seq_num_index) = sn.value();
  received_seq_num_time(received_seq_num_index) = Event_Queue::now();
  received_seq_num_index++;
}


void TCP_Receiver::save_trace(std::string filename)
{

  received_seq_num_val.set_size(received_seq_num_index, true);
  received_seq_num_time.set_size(received_seq_num_index, true);

  if (fDebug) {
    std::cout << "received_seq_num_val" << received_seq_num_val << std::endl;
    std::cout << "received_seq_num_time" << received_seq_num_time << std::endl;
    std::cout << "received_seq_num_index" << received_seq_num_index << std::endl;
    std::cout << "TCP_Receiver::saving to file: " << filename << std::endl;
  }

  it_file ff2;
  ff2.open(filename);

  ff2 << Name("received_seq_num_val") << received_seq_num_val;
  ff2 << Name("received_seq_num_time") << received_seq_num_time;
  ff2 << Name("received_seq_num_index") << received_seq_num_index;

  ff2.flush();
  ff2.close();
}


} //namespace itpp

#ifdef _MSC_VER
#pragma warning(default:4355)
#endif

//! \endcond
