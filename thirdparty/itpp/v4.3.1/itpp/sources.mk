h_sources = $(top_srcdir)/itpp/itbase.h \
	$(top_srcdir)/itpp/itstat.h \
	$(top_srcdir)/itpp/itmex.h

if ENABLE_COMM
  h_sources += $(top_srcdir)/itpp/itcomm.h
endif
if ENABLE_FIXED
  h_sources += $(top_srcdir)/itpp/itfixed.h
endif
if ENABLE_OPTIM
  h_sources += $(top_srcdir)/itpp/itoptim.h
endif
if ENABLE_PROTOCOL
  h_sources += $(top_srcdir)/itpp/itprotocol.h
endif
if ENABLE_SIGNAL
  h_sources += $(top_srcdir)/itpp/itsignal.h
endif
if ENABLE_SRCCODE
  h_sources += $(top_srcdir)/itpp/itsrccode.h
endif
