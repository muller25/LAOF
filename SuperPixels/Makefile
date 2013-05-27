INCDIR = -I./vlib/include -I./
DBG    = -g
OPT    = -O3
CPP    = g++
CFLAGS = $(DBG) $(OPT) $(INCDIR) -fPIC
OPENCVFLAGS = `pkg-config --cflags opencv`
LINK   = -lm
OPENCVLIBS = `pkg-config --libs opencv`
OBJS	= graph.o maxflow.o SuperPixels.o
TARGET  = superpixels

all: $(TARGET)
	echo "Done"

superpixels: $(OBJS) test_superpixels.o
	$(CPP) -o superpixels $^ $(LINK) $(OPENCVLIBS)

graph.o: graph.cpp graph.h
	$(CPP) $(CFLAGS) -c graph.cpp -o graph.o

maxflow.o: maxflow.cpp 
	$(CPP) $(CFLAGS) -c maxflow.cpp -o maxflow.o

SuperPixels.o: SuperPixels.cpp SuperPixels.h graph.h energy.h
	$(CPP) $(CFLAGS) $(OPENCVFLAGS) -c $< -o $@

test_superpixels.o: test_SuperPixels.cpp
	$(CPP) $(CFLAGS) $(OPENCVFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) *.o cscope.out
