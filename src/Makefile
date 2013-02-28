CXX				= g++
CXXFLAGS  		= -O3 -Wall
INCDIR			= .

LINK			= g++
LFLAGS			= 
LIBDIR 			= .

DEL				= rm -rf

OPENCVLIBS		= `pkg-config --libs opencv`
OBJS			= Maths.o OpticalFlow.o main.o Flow2Color.o GaussianPyramid.o
TARGET			= run

all: $(TARGET)
	@echo "DONE"

run: $(OBJS)
	$(LINK) $(LFAGS) -o run $(OBJS) $(OPENCVLIBS)

main.o: main.cpp
	$(CXX) $(CFLAGS) -c main.cpp -o main.o

%.o: %.cpp %.h
	$(CXX) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) $(OBJS) $(TARGET)
