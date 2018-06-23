CFLAGS += -Wall -O -std=gnu1x `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

.PHONY: all
all: test1.out test2.out times.svg true.svg type2.svg

true.svg: true.gpi true2.out
	gnuplot $<

type2.svg: int2.gpi int2.out
	gnuplot $<

true2.out int2.out: test2test
	./$<

test2.out: test2
	./$< > $@ 100

test1.out: test1
	./$< > $@ 100

times.svg: plot.gpi times1.out times2.out
	gnuplot $<

times1.out: test1 Makefile
	rm -f times1.out
	for n in `seq 1000000 1000000 10000000`; do time --format "$$n %U" --append --output $@ ./test1 $$n 1>/dev/null 2>/dev/null && echo $$n; done

times2.out: test2 Makefile
	rm -f times2.out
	for n in `seq 1000000 1000000 10000000`; do time --format "$$n %U" --append --output $@ ./test2 $$n 1>/dev/null 2>/dev/null && echo $$n; done

test1: test1.o bilinear.o

test2: test2.o bilin_type2.o

test2test: test2test.o bilin_type2.o



clean:
	find . -type f -executable -delete
	rm -f *.o
	rm -f *.out
	rm -f *.svg
