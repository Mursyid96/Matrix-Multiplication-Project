CC 	= gcc
CFLAGS	= -Wall -O3 -fopenmp
LDFLAGS	= -lm

CC2 = mpicc
CFLAGS2 = -g -Wall

all: TaskB-1 TaskB-2

TaskB-1:  TaskB-1.c
	$(CC2) $(CFLAGS) -o $@ $? $(LDFLAGS)

TaskB-2:  TaskB-2.c
	$(CC2) $(CFLAGS) -o $@ $? $(LDFLAGS)

.SUFFIXES:	.o .c

%.o : %.c
	$(CC) $(CFLAGS) -c $<
%.s : %.c
	$(CC) $(CFLAGS) -S $<

clean:
	/bin/rm -f *.o *~ TaskB-1 TaskB-2

