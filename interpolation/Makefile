CFLAGS = -Wall -std=gnu1x
LDLIBS = -lm

linintplot.svg: linint.gpi ldat.txt
	gnuplot $<

linplot.svg: linplot.gpi ldat.txt
	gnuplot $<

ldat.txt: lintest
	./$< > ldat.txt

lintest: lintest.o search.o lspline.o

clean:
	find . -type f -executable -delete
	rm *.o
	rm *.txt
	rm *.svg
