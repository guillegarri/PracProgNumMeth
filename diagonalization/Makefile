CFLAGS += -Wall -O -std=gnu1x `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

.PHONY: all
all: times.txt sweepsout.txt eigenbyeigenout.txt times.svg

times.txt: main Makefile
	rm -f times.txt
	for n in `seq 200 2 300`; do time --format "$$n %U" --append --output $@ ./main $$n 1>/dev/null; done

eigenbyeigenout.txt: eigenbyeigen
	./$< > $@

sweepsout.txt: main
	./$< > $@

eigenbyeigen: eigenbyeigen.o jacobi.o

main: main.o jacobi.o

times.svg: timeplot.gpi times.txt
	gnuplot $<

clean:
	find . -type f -executable -delete
	rm -f *.o
	rm -f *.txt
	rm -f *.svg
