
IRRAM ?= /usr/local

CPPFLAGS = -I$(IRRAM)/include -Iinclude
CXXFLAGS = -Wall -std=c++11 -O2
LDFLAGS  = -L$(IRRAM)/lib -Wl,-rpath,$(IRRAM)/lib
LDLIBS   = -liRRAM -lmpfr -lgmp

BIN = \
	examples/plot

LIB_OBJS = \
	src/random-core.o \
	src/random-real.o \
	src/random-matrix.o \
	src/random-complex.o \

BIN_OBJS = $(BIN:%=%.o)

OBJS = $(LIB_OBJS) $(BIN_OBJS)

LIB = librandom.a

all: $(BIN)

$(BIN): %: %.o $(LIB)

# link binaries using the C++ compiler
$(BIN): CC = $(CXX)

# make static library (archive) from $(OBJS)
$(LIB): $(patsubst %,$(LIB)(%),$(LIB_OBJS))

# all objects depend on headers (more fine-grained deps via -MMD)
$(OBJS): %.o: $(wildcard include/*.h)

clean:
	$(RM) $(BIN) $(BIN_OBJS) $(LIB) $(LIB_OBJS)

.PHONY: all clean
