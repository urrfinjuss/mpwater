FSW     = foobar.x

#List of sources
C_SRC   = foo.c
CPP_SRC = main.cpp

# Compilers -Ofast -ffast-math -mfpmath=387
CC      = gcc
CPP	= g++
LINK    = g++
C_OPT     = -march=native -std=gnu99 -Ofast -flto -ffast-math -mfpmath=sse+387 -Wall -Wpointer-arith -Wcast-align -fno-schedule-insns -fschedule-insns2 -fstrict-aliasing -funroll-loops -fprefetch-loop-arrays

CPP_OPT     = -march=native -Ofast -flto -ffast-math -mfpmath=sse+387 -Wall -Wpointer-arith -Wcast-align -fno-schedule-insns -fschedule-insns2 -fstrict-aliasing -funroll-loops -fprefetch-loop-arrays

#-----------------------------
#generic

LIB_MPI         =
LIB_FFT         = -L$(HOME)/usr/lib  -lmpfft -lmpfr
INC_MPI         =
INC_FFT         = -I$(HOME)/usr/include
LIB_ADD         =

#-----------------------------

C_OBJ_FSW         = $(C_SRC:.c=.o) $(FC_SRC:.f=.o)
CPP_OBJ_FSW       = $(CPP_SRC:.cpp=.o) $(FCPP_SRC:.f=.o)
LIB_FSW         = $(LIB_MPI) $(LIB_FFT) $(LIB_ADD)
INC_FSW         = $(INC_MPI) $(INC_FFT)

#-----------------------------

default: foobar

.f.o:
	$(FC) $(FFLAGS) -c $<


foobar:
	$(CC) $(C_OPT) $(DEF_FSW) $(INC_FSW) -c $(C_SRC)
	$(CPP) $(CPP_OPT) $(DEF_FSW) $(INC_FSW) -c $(CPP_SRC)
	$(LINK) $(CPP_OPT) $(OBJ_FSW) $(LIB_FSW) -o $(FSW) $(C_OBJ_FSW) $(CPP_OBJ_FSW)

hostname:
	@echo $(HOSTNAME) $(INC_FFT)

clean:
	@echo "cleaning ..."
	rm -f *~ *.o

#-------------------------------


