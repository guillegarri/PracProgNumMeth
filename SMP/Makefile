CFLAGS += -Wall -O -std=gnu1x `gsl-config --cflags` -fopenmp
LDLIBS += `gsl-config --libs` -lgomp

.PHONEY: all
all: comparison.txt error.svg results.out

comparison.txt: main mainnoSMP Makefile
	echo "Both parallel and not parallel take longer than my normal not parallel, since they need to use same RNG, which must be thread-safe.\n" > $@
	echo "SMP:" >> $@
	time --append --output $@ ./main 1>/dev/null 2>/dev/null
	echo "\nnot SMP:" >> $@
	time --append --output $@ ./mainnoSMP 1>/dev/null 2>/dev/null

error.svg: plot.gpi error.out
	gnuplot $<

results.out error.out: main
	./$< 1> results.out 2> error.out

main: main.o montecarlo.o

mainnoSMP: mainnoSMP.o montecarlonoSMP.o

clean:
	find . -type f -executable -delete
	rm -f *.o
	rm -f *.out
