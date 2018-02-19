CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
MADEFILES = main-stdio main-stdio.o main-cdmline main-cmdline.o plot.svg
TEXTS = out1.txt out.io.txt out.cmd.txt input.data.txt

.PHONEY: all clean progs
progs : main-cmdline main-stdio

all              : out.io.txt out.cmd.txt plot.svg         ;
	cat out.io.txt
	cat out.cmd.txt

out.io.txt       : main-stdio input.data.txt  ;
	cat input.data.txt | ./main-stdio > out.io.txt
out1.txt         : main-stdio              ; ./$< > $@
main-stdio       : main-stdio.o
main-stdio.o     : main-stdio.c         # built-in compiling rule is used here

out.cmd.txt         : main-cmdline input.data.txt   ;
	./main-cmdline `cat input.data.txt` > out.cmd.txt
main-cmdline		 : main-cmdline.o
main-cmdline.o   : main-cmdline.c

input.data.txt   : ;
	echo `seq 0 0.2 5` > input.data.txt

plot.svg         : plot.gpi out.io.txt out.cmd.txt;
	gnuplot plot.gpi

clean : ; $(RM) $(MADEFILES) $(TEXTS)
