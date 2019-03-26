# FFCFLAGS = -O -r quadrature -fsplit -l dolfin

FFCFLAGS = -l dolfin

all:
ifeq (, $(shell which ffc))
	$(error "No ffc in $(PATH), consider installing it.")
else
	ffc $(FFCFLAGS) *.ufl 
endif
