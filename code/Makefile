CC = g++ -std=c++0x
CFLAGS = -O3 #-g -Wall
DFLAGS = -g -pg
LDFLAGS = -lm

PROG = scheduler
DEBUG_PROG = scheduler-debug
HDRS = scheduler.h   push_relabel.h   util.h
SRCS = scheduler.cpp push_relabel.cpp

## This incantation says that the object files
## have the same name as the .cpp files, but with .o
OBJS = $(SRCS:.cpp=.o)
DOBJS = $(SRCS:.cpp=.debug)

###############################################################################

## This is the first rule (the default)
default : $(PROG)
debug : $(DEBUG_PROG)
all : $(PROG) $(DEBUG_PROG)

## Build the program from the .o's
$(PROG) : $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(PROG) $(LDFLAGS)

$(DEBUG_PROG) : $(DOBJS)
	$(CC) $(DFLAGS) $(DOBJS) -o $(DEBUG_PROG) $(LDFLAGS)

## Rules for the source files
%.o : %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

%.debug : %.cpp
	$(CC) $(DFLAGS) -c -o $@ $<

## Remove all the compilation and debugging files
clean :
	rm -f core $(OBJS) $(DOBJS)

