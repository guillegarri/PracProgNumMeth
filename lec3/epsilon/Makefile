

.PHONEY:default
default: out.txt
	cat out.txt

.PHONEY:fun
fun: equaltester digitwriter


madefiles = out.txt minint maxint maxint.txt minint.txt epsilon epsilon.txt sums sums.txt equaltests equaltests.txt equaltester digitwriter
texts = maxint.txt minint.txt epsilon.txt sums.txt equaltests.txt

out.txt: $(texts)
	cat maxint.txt > out.txt
	cat minint.txt >> out.txt
	cat epsilon.txt >> out.txt
	cat sums.txt >> out.txt
	cat equaltests.txt >> out.txt

maxint.txt: maxint
	./maxint > maxint.txt

minint.txt: minint
	./minint > minint.txt

epsilon.txt: epsilon
	./epsilon > epsilon.txt

sums.txt: sums
	./sums > sums.txt

equaltests.txt: equaltests
	./equaltests > equaltests.txt

maxint: maxint.c
	gcc maxint.c -o maxint

minint: minint.c
	gcc minint.c -o minint

epsilon: epsilon.c
	gcc epsilon.c -o epsilon

sums: sums.c
	gcc sums.c -o sums

equaltests: equaltests.c equal.c
	gcc equal.c equaltests.c -o equaltests

clean:
	$(RM) $(madefiles)

#Just for fun
equaltester: equaltester.c equal.c
	gcc equal.c equaltester.c -o equaltester

digitwriter: digitwriter.c
	gcc digitwriter.c -o digitwriter
