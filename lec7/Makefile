CFLAGS = -Wall `gsl-config --cflags`
LDLIBS = -lm `gsl-config --libs`
MADEFILES = airydata.dat airycalc airycalc.o linsolv linsolv.o linres.txt

.PHONY: all clean

all: airyplot.svg linres.txt

airyplot.svg: airyplot.gpi airydata.dat
	gnuplot airyplot.gpi

airydata.dat: airycalc
	./$< > $@

airycalc: airycalc.o

airycalc.o: airycalc.c

linres.txt: linsolv
	./$< > $@

linsolv: linsolv.o

linsolv.o: linsolv.c

clean:
	$(RM) $(MADEFILES)
