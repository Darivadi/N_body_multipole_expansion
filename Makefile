CC = gcc
CFLAGS = -O3 -I/home/$(USER)/libs/include/
CFLAGSDEBUG = -g -Wall -I/home/$(USER)/libs/include/
LFLAGS = -L/home/$(USER)/libs/lib
PROGRAM = main_N_body_multipolar

$(PROGRAM):
	echo "Compiling"
	$(CC) -c $@.c $(CFLAGS)
	$(CC) $@.o -lm $(LFLAGS) -o $@.x

debug:
	$(CC) $(PROGRAM).c $(CFLAGSDEBUG)
	$(CC) $(PROGRAM).o -lm $(LFLAGS) -o $(PROGRAM).x

clean:
	rm -rf *.out
	rm -rf *-
	rm -rf *.out
	rm -rf *#
	rm -rf *.o
	rm -rf *.a
	rm -rf *.so
	rm *.x
