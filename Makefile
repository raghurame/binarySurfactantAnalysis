compile:
	gcc -o multistageOrdering multistageOrdering.c -lm -fopenmp -Wall
	gcc -o assignMolID assignMolID.c -lm -Wall
run:
	./multistageOrdering surfactantVisual.lammpstrj latest.data
