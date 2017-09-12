EXECUTABLE     = foobar.x

# Compilers:
CC      = gcc
CPP	= g++
LINK    = g++

# Common flags for C and C++ compiler:
CFLAGS =-march=native -Ofast -flto -ffast-math -Wall -fno-schedule-insns -fschedule-insns2 -fstrict-aliasing -funroll-loops -fprefetch-loop-arrays
C_EXTRA =-std=gnu99

#-----------------------------
# Paths:
SOURCEDIR = src
BUILDDIR  = build

# Extra libraries and headers:
LIB	= -L$(HOME)/usr/lib  -lmpfft -lmpfr
INC	= -I$(HOME)/usr/include

# List of sources:
C_SOURCES   = $(wildcard $(SOURCEDIR)/*.c)
CPP_SOURCES = $(wildcard $(SOURCEDIR)/*.cpp)

# List of objects:
C_OBJECTS	= $(patsubst $(SOURCEDIR)/%.c,$(BUILDDIR)/%.o,$(C_SOURCES))
CPP_OBJECTS	= $(patsubst $(SOURCEDIR)/%.cpp,$(BUILDDIR)/%.o,$(CPP_SOURCES))

.PHONY: clean

all: dir $(BUILDDIR)/$(EXECUTABLE)

dir:
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/$(EXECUTABLE): $(C_OBJECTS) $(CPP_OBJECTS) 
	$(LINK) $(CFLAGS) $(LIB) $^ -o $@

$(CPP_OBJECTS): $(BUILDDIR)/%.o : $(SOURCEDIR)/%.cpp
	$(CPP) $(CFLAGS) $(INC) -c $< -o $@

$(C_OBJECTS): $(BUILDDIR)/%.o : $(SOURCEDIR)/%.c
	$(CC) $(CFLAGS) $(C_EXTRA) $(INC) -c $< -o $@

clean:
	@echo "cleaning ..."
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/$(EXECUTABLE)

#-------------------------------


