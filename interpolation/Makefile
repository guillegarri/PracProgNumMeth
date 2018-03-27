CFLAGS = -Wall -std=gnu1x
LDLIBS = -lm
allplots = linintplot.svg linplot.svg qintplot.svg qderivplot.svg qplot.svg ycompplot.svg intcompplot.svg dercompplot.svg

.PHONY: plots
plots: $(allplots)

qintplot.svg: qintplot.gpi qdat.txt
	gnuplot $<

qderivplot.svg: qderivplot.gpi qdat.txt
		gnuplot $<

qplot.svg: qplot.gpi qdat.txt
	gnuplot $<

linintplot.svg: linint.gpi ldat.txt
	gnuplot $<

linplot.svg: linplot.gpi ldat.txt
	gnuplot $<

ycompplot.svg: ycompplot.gpi alltestdat.txt
	gnuplot $<

intcompplot.svg: intcompplot.gpi alltestdat.txt
	gnuplot $<

dercompplot.svg: dercompplot.gpi alltestdat.txt
	gnuplot $<

ldat.txt: lintest
	./$< > $@

qdat.txt: qtest
	./$< > $@

cdat.txt: ctest
	./$< > $@

qtest.txt: qtest2
	./$< > $@

alltestdat.txt: alltest
	./$< > $@

alltest: cspline.o qspline.o lspline.o search.o

qtest2: qtest2.o qspline.o

ctest: ctest.o cspline.o

qtest: qtest.o qspline.o

lintest: lintest.o search.o lspline.o

clean:
	find . -type f -executable -delete
	rm *.o
	rm *.txt
	rm *.svg
