CFLAGS  = `pkg-config --cflags dolfin`
LIBS    = `pkg-config --libs dolfin`
CXX     = `pkg-config --variable=compiler dolfin` -g

DEST    = demo
OBJECTS = main.o ufc/Bloch_Torrey3D.o ufc/Bloch_Torrey_NoTime3D.o ufc/Comp_Sig3D.o ufc/InitialCondition3D.o 

all: $(DEST)

install:

clean:
	-rm -rf results *.o core *.core $(OBJECTS) $(DEST)

$(DEST): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(CFLAGS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -c $< -o $@
