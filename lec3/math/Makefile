CFLAGS = -std=gnu99 -Wall -O

.PHONEY:default
default: out.txt
	cat out.txt

out.txt: math
	./math > out.txt

math: math.c
	gcc -o math math.c -lm

clean:
	$(RM) out.txt math
