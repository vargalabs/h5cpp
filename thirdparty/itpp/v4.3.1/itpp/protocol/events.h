/*!
 * \file
 * \brief Definitions of an event-based simulation class
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

#ifndef EVENTS_H
#define EVENTS_H

#include <itpp/base/itassert.h>

#if (defined(_MSC_VER) && defined(ITPP_SHARED_LIB) && !(defined(itpp_EXPORTS) || defined(itpp_debug_EXPORTS)))

#ifndef ITPP_PROTOCOL_EXCLUDED
#define ITPP_PROTOCOL_EXCLUDED
#pragma message( "PROTOCOL definitions are not available for MSVC shared builds" )
#endif

#else


#include <queue>
#include <deque>

namespace itpp
{

//! \addtogroup protocol
//@{

// typedef long double Ttype; //!< 128-bit floating point time
typedef double Ttype; //!< 64-bit floating point time
// typedef long unsigned int Ttype; //!< 64-bit unsigned integer time

class Event_Queue;
class Base_Event;
class Base_Signal;

/*!
  \brief Base Event Class

  An abstract Base class of Events that can be used to derive new events. All Event classes
  need to define the exec() function which is called when the event expires. An event has an
  execution time and an id.
*/
class Base_Event
{
public:
  friend class Base_Signal;

  friend class Event_Queue;

  friend struct Compare_Base_Event_Times;

  //! Schedule an event at time \c delta_time from now
  Base_Event(const Ttype delta_time) {  // The event will occur in 'delta_time' time units from now!
    it_assert(delta_time >= 0, "Only causal simulations are possible");
    active = true;
    delta_t = delta_time;
    expire_t = 0; // Will be set correctly upon calling Event_Queue::add().
    id = global_id++;
  }

  //! Destructor
  virtual ~Base_Event() {}

  //! Cancel an event
  void cancel() { active = false; }

protected:
  //! ADD DOCUMENTATION HERE
  virtual void exec(void) = 0;
  //! ADD DOCUMENTATION HERE
  Ttype delta_t;
  //! ADD DOCUMENTATION HERE
  Ttype expire_t;
  //! ADD DOCUMENTATION HERE
  bool active;
  //! ADD DOCUMENTATION HERE
  unsigned long long int id;
  //! ADD DOCUMENTATION HERE
  static unsigned long long int global_id;
};

//! Compare to events, Returns true if expire time of event1 is larger than the expire time of event2
struct Compare_Base_Event_Times {
  //! ADD DOCUMENTATION HERE
  bool operator()(Base_Event *event1, Base_Event *event2) {
    if (event1->expire_t == event2->expire_t) // Equal expire times.
      return (event1->id > event2->id); // Base comparison on the event id.
    else
      return (event1->expire_t > event2->expire_t); // Different expire times. Regular comparison.
  }
};

/*!
  \brief Event Queue class

  A class for storing and executing events. Events can be added to the queue and when the start() is
  called all events will be executed. Observe that Events need to be created before they are added to the
  queue by calling an appropriate constructor. However, expired events are destroyed automatically (the destructor is called).

*/
class Event_Queue
{
public:
  friend class Base_Signal;

  //! Constructor
  Event_Queue() {}
  //! Destructor
  ~Event_Queue() {}

  //! Add event to Queue
  static void add(Base_Event *e);
  //! Return current time
  static Ttype now() {return t;}
  //! Start executing events
  static void start();
  //! Stop execution of events
  static void stop();
  //! Remove all events
  static void clear();
protected:
  //static void cancel_all(Base_Signal *s);
private:
  typedef std::deque<Base_Event*, std::allocator< Base_Event* > >::iterator Base_Event_Iterator;
  static void _run();
  static bool keep_running;
  static Ttype t; // Current time.
  static std::priority_queue < Base_Event*,
  std::deque<Base_Event*, std::allocator<Base_Event*> >,
  Compare_Base_Event_Times > event_queue; // Queue for the Events.
};

/*!
  \brief An Event class that executes a function when the event expires.

  Since Events are objects you need supply both a pointer to the object and the function pointer to create the Event
*/
template <class ObjectType>
class Event : public Base_Event
{
public:
  //! Construct an Event to expire delta_time from now by calling the function (*object_pointer.*object_function_pointer)()
  Event(ObjectType *object_pointer, void (ObjectType::*object_function_pointer)(), const Ttype delta_time) : Base_Event(delta_time) {
    po = object_pointer;
    pm = object_function_pointer;
  }

  //! Destructor
  virtual ~Event() {}

  //! Execute (call) the assigned function
  virtual void exec(void) {(*po.*pm)(); }

private:
  void (ObjectType::*pm)(); // Pointer to class member function to be executed on event expire.
  ObjectType *po; // Pointer to object who's member function is to be executed on event expire.
};

/*!
  \brief An Event class that executes a function with some data as input when the event expires.

  Since Events are objects you need supply both a pointer to the object and the function pointer to create the Event
*/
template <class ObjectType, class DataType> class Data_Event : public Base_Event
{
public:
  //! Construct an Event to expire delta_time from now by calling the function (*object_pointer.*object_function_pointer)(data)
  Data_Event(ObjectType *object_pointer,
             void (ObjectType::*object_function_pointer)(DataType data),
             DataType data, const Ttype delta_time) : Base_Event(delta_time) {
    po = object_pointer;
    pm = object_function_pointer;
    u = data;
  }

  //! Destructor
  virtual ~Data_Event() {}

  //! Execute (call) the assigned function with user data.
  virtual void exec(void) {
    (*po.*pm)(u);
  }

private:
  void (ObjectType::*pm)(DataType data); // Pointer to class member function to be executed on event expire.
  ObjectType* po; // Pointer to object who's member function is to be executed on event expire.
  DataType u; // User data.
};

//@}

} // namespace itpp

#endif

#endif // #ifndef EVENTS_H
