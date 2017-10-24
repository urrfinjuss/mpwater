FSW     = fsw_dw.x

#List of sources
C_FSW   = main.c memory.c array_func.c input.c output.c mapping.c hlevel.c evolve.c pade.c

# Compilers -Ofast -ffast-math -mfpmath=387
CC      = gcc
LINK    = gcc
OPT     = -march=native -std=gnu99 -Ofast -flto -ffast-math -mfpmath=sse+387 -Wall -Wpointer-arith -Wcast-align -fno-schedule-insns -fschedule-insns2 -fstrict-aliasing -funroll-loops -fprefetch-loop-arrays
#OPT	= -march=native -O -fno-inline -Wall -std=gnu99
#OPT	= -std=gnu99 -mfpmath=387 -Wall

#-----------------------------
#generic

LIB_MPI         =
LIB_FFT         = -L$(HOME)/usr/lib  -lfftw3l_threads -lfftw3l -lm -lpthread 
INC_MPI         =
INC_FFT         = -I$(HOME)/usr/include
LIB_ADD         =

#-----------------------------

OBJ_FSW         = $(C_FSW:.c=.o) $(F_FSW:.f=.o)
LIB_FSW         = $(LIB_MPI) $(LIB_FFT) $(LIB_ADD)
INC_FSW         = $(INC_MPI) $(INC_FFT)

#-----------------------------

default: fsw

.f.o:
	$(FC) $(FFLAGS) -c $<


fsw:
	$(CC) $(OPT) $(DEF_FSW) $(INC_FSW) -c $(C_FSW)
	$(LINK) $(OPT) $(OBJ_FSW) $(LIB_FSW) -o $(FSW)
	cp fsw_dw.x ./debug/fsw_dw.x

hostname:
	@echo $(HOSTNAME) $(INC_FFT)

clean:
	@echo "cleaning ..."
	rm -f *~ *.o
	rm -f ./debug/fsw_dw.x

#-------------------------------


