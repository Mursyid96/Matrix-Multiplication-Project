CC 	= gcc
CFLAGS	= -Wall -O3 -fopenmp
CFLAGS2 = -Wall -O3
LDFLAGS	= -lm
LDFLAGS2 = -lm -lpthread

CC2 = mpicc
CFLAGS2 = -g -Wall

all: Task1-1 Task1-2 Task2-1 Task2-2 Task3-1 Task3-2 Task4A-1 Task4A-2 Task4B-1 Task4B-2

Task1-1:	Task1-1.c
	$(CC) $(CFLAGS) -o $@ $? $(LDFLAGS)

Task1-2:	Task1-2.c
	$(CC) $(CFLAGS) -o $@ $? $(LDFLAGS)

Task2-1:  Task2-1.c
	$(CC2) $(CFLAGS2) -o $@ $? $(LDFLAGS)

Task2-2:  Task2-2.c
	$(CC2) $(CFLAGS2) -o $@ $? $(LDFLAGS)

Task3-1:  Task3-1.c
	$(CC2) $(CFLAGS2) -o $@ $? $(LDFLAGS)

Task3-2:  Task3-2.c
	$(CC2) $(CFLAGS2) -o $@ $? $(LDFLAGS)

Task4A-1:	Task4A-1.c
	$(CC) $(CFLAGS2) -o $@ $? $(LDFLAGS2)

Task4A-2:	Task4A-2.c
	$(CC) $(CFLAGS2) -o $@ $? $(LDFLAGS2)

Task4B-1:	Task4B-1.c
	$(CC) $(CFLAGS2) -o $@ $? $(LDFLAGS2)

Task4B-2:	Task4B-2.c
	$(CC) $(CFLAGS2) -o $@ $? $(LDFLAGS2)

.SUFFIXES:	.o .c

%.o : %.c
	$(CC) $(CFLAGS) -c $<
%.s : %.c
	$(CC) $(CFLAGS) -S $<

clean:
	/bin/rm -f *.o *~ Task1-1 Task1-2 Task2-1 Task2-2 Task3-1 Task3-2 Task4A-1 Task4A-2 Task4B-1 Task4B-2

