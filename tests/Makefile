CXX				= g++
CXXFLAGS  		= -O3 -Wall
INCDIR			= ../include

LINK			= g++
LFLAGS			= 
LIBDIR 			= ../libs
OPENCVLIBS		= `pkg-config --libs opencv`

DEL				= rm -rf

LAOF			= ../src
OBJS			= test_Maths.o test_OpticalFlow.o
TARGET 			= run_all_tests flow2color pyramid

all: $(TARGET)
	@echo "DONE"

run_all_tests: $(OBJS) $(LAOF)/OpticalFlow.o $(LAOF)/Maths.o $(LAOF)/GaussianPyramid.o
	$(LINK) $(LFLAGS) -o $@ $^ -L$(LIBDIR) -lgtest -lgtest_main -lpthread $(OPENCVLIBS)

flow2color: test_Flow2Color.o $(LAOF)/Maths.o $(LAOF)/Flow2Color.o
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

pyramid: test_GaussianPyramid.o $(LAOF)/GaussianPyramid.o
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -I$(LAOF) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) *.o $(TARGET)
