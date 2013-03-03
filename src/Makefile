CXX				= g++
CXXFLAGS  		= -O3 -Wall
INCDIR			= .

LINK			= g++
LFLAGS			= 
LIBDIR 			= .

DEL				= rm -rf

OPENCVLIBS		= `pkg-config --libs opencv`
OBJS			= OpticalFlow.o GaussianPyramid.o Flow2Color.o main.o
TARGET			= run

all: $(TARGET)
	@echo "DONE"

run: $(OBJS)
	$(LINK) $(LFAGS) -o run $(OBJS) $(OPENCVLIBS)

main.o: main.cpp
	$(CXX) $(CFLAGS) -c main.cpp -o main.o

OpticalFlow.o: OpticalFlow.cpp OpticalFlow.h Maths.h Image.h GaussianPyramid.h ImageProcess.h
	$(CXX) $(CFLAGS) -c $< -o $@

Flow2Color.o: Flow2Color.cpp Flow2Color.h Image.h
	$(CXX) $(CFLAGS) -c $< -o $@

%.o: %.cpp %.h
	$(CXX) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) $(OBJS) $(TARGET) cscope.out
