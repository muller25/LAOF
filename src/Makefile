CXX				= g++
CXXFLAGS  		= -O3 -Wall 
INCDIR			= `pkg-config --cflags opencv`

LINK			= g++
LFLAGS			= 
LIBDIR 			= .

DEL				= rm -rf

GC				= ../gc
OPENCVLIBS		= `pkg-config --libs opencv`
OBJS			= main.o OpticalFlow.o GaussianPyramid.o Flow2Color.o \
MotionLayers.o LAOF.o
TARGET			= run

all: $(TARGET)
	@echo "DONE"

run: $(OBJS) $(GC)/GCoptimization.o $(GC)/LinkedBlockList.o
	$(LINK) $(LFAGS) -o run $^ $(OPENCVLIBS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) -c main.cpp -o main.o

OpticalFlow.o: OpticalFlow.cpp OpticalFlow.h Maths.h Image.h GaussianPyramid.h ImageProcess.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

Flow2Color.o: Flow2Color.cpp Flow2Color.h Image.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

MotionLayers.o: MotionLayers.cpp MotionLayers.h ML.h Image.h Maths.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(GC) -c $< -o $@

LAOF.o: LAOF.cpp LAOF.h Image.h OpticalFlow.h MotionLayers.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) $(OBJS) $(TARGET) cscope.out *.yml
