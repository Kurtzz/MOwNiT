main:
	gcc -Wall -g main.c `pkg-config --cflags --libs gsl` -o main

mainO:
	gcc -O2 -Wall -g main.c `pkg-config --cflags --libs gsl` -o main

clean:
	rm -f main

testt:
	./main 100
	./main 150
	./main 200
	./main 250
	./main 300
	./main 350
	./main 400
	./main 450
	./main 500
	./main 550
	./main 600
	./main 650
	./main 700
