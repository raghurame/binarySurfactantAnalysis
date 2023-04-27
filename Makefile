compile:
	gcc -o multistageOrdering multistageOrdering.c -lm -fopenmp
run:
	./multistageOrdering surfactantVisual.lammpstrj latest.data
