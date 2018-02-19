CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
.PHONEY: all clean
all              : out.txt           ; cat $<
out.txt          : main              ; ./$< > $@
main             : main.o nvector.o  # built-in linking rule is used here
main.o nvector.o : nvector.h         # built-in compiling rule is used here
clean            :                   ; $(RM) main main.o nvector.o out.txt
