gsl: clean
	gcc interpol_gsl.c `pkg-config --cflags --libs gsl` -o gsl
	gcc interpol_times.c `pkg-config --cflags --libs gsl` -o times

clean:
	rm -f gsl times
