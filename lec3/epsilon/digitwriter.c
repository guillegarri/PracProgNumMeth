#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"

void name_digit(int i){
  switch (i) {
    case 1:
    printf("one");
    break;
    case 2:
    printf("two");
    break;
    case 3:
    printf("three");
    break;
    case 4:
    printf("four");
    break;
    case 5:
    printf("five");
    break;
    case 6:
    printf("six");
    break;
    case 7:
    printf("seven");
    break;
    case 8:
    printf("eight");
    break;
    case 9:
    printf("nine");
    break;
    case 0:
    printf("zero");
    break;
  }

}

int main(int argc, char *argv[]) {
  if (argc ==2 && atoi(argv[1])<10 && atoi(argv[1])>= 0) {
    printf("You have provided ");
    name_digit(atoi(argv[1]));
    printf(".\n");
  } else {
    printf("You must provide exactly one digit.\n");
  }

  return 0;
}
