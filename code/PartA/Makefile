CC 	= gcc
CFLAGS	= -Wall -O3 -fopenmp
CFLAGS2 = -Wall -O3
LDFLAGS	= -lm
LDFLAGS2 = -lm -lpthread

CC2 = mpicc
CFLAGS2 = -g -Wall

all: Task1 Task2 Task3 Task4A Task4B

Task1:	Task1.c
	$(CC) $(CFLAGS) -o $@ $? $(LDFLAGS)

Task2:  Task2.c
	$(CC2) $(CFLAGS2) -o $@ $? $(LDFLAGS)

Task3:  Task3.c
	$(CC2) $(CFLAGS2) -o $@ $? $(LDFLAGS)

Task4A:	Task4A.c
	$(CC) $(CFLAGS2) -o $@ $? $(LDFLAGS2)

Task4B:	Task4B.c
	$(CC) $(CFLAGS2) -o $@ $? $(LDFLAGS2)

.SUFFIXES:	.o .c

%.o : %.c
	$(CC) $(CFLAGS) -c $<
%.s : %.c
	$(CC) $(CFLAGS) -S $<

clean:
	/bin/rm -f *.o *~ Task1 Task2 Task3 Task4A

