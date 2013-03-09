CXX				= g++
CXXFLAGS  		= -O3 -Wall
INCDIR			= -I../include `pkg-config --cflags opencv`

LINK			= g++
LFLAGS			= 
LIBDIR 			= ../libs
OPENCVLIBS		= `pkg-config --libs opencv`

DEL				= rm -rf

LAOF			= ../src
OBJS			= test_Maths.o test_OpticalFlow.o
TARGET 			= flow2color pyramid imop kmeans #run_all_tests

all: $(TARGET)
	@echo "DONE"

run_all_tests: $(OBJS) $(LAOF)/OpticalFlow.o $(LAOF)/GaussianPyramid.o
	$(LINK) $(LFLAGS) -o $@ $^ -L$(LIBDIR) -lgtest -lgtest_main -lpthread $(OPENCVLIBS)

flow2color: test_Flow2Color.o $(LAOF)/Flow2Color.o
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

pyramid: test_GaussianPyramid.o $(LAOF)/GaussianPyramid.o
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

test_ImageOP.o: test_ImageOP.cpp $(LAOF)/ImageProcess.h $(LAOF)/ImageIO.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(LAOF) -c $< -o $@

imop: test_ImageOP.o 
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

test_kmeans.o: test_kmeans.cpp $(LAOF)/ML.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(LAOF) -c $< -o $@

kmeans: test_kmeans.o 
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(LAOF) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) *.o $(TARGET) *.jpg *.yml
