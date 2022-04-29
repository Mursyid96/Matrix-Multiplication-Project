CC 	= gcc
CFLAGS	= -Wall -O3 -fopenmp
LDFLAGS	= -lm

all: Task1

Task1:	Task1.c
	$(CC) $(CFLAGS) -o $@ $? $(LDFLAGS)

.SUFFIXES:	.o .c

%.o : %.c
	$(CC) $(CFLAGS) -c $<
%.s : %.c
	$(CC) $(CFLAGS) -S $<

clean:
	/bin/rm -f *.o *~ Task1

