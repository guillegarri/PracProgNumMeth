CFLAGS += -Wall -O -std=gnu1x `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

out.txt: main
	./$< > $@

main: main.o gramschmidtfuncs.o

clean:
	find . -type f -executable -delete
	rm -f *.o
