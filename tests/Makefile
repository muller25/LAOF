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
TARGET 			= run_all_tests flow2color pyramid imop

all: pyramid imop flow2color
	@echo "DONE"

run_all_tests: $(OBJS) $(LAOF)/OpticalFlow.o $(LAOF)/Maths.o $(LAOF)/GaussianPyramid.o
	$(LINK) $(LFLAGS) -o $@ $^ -L$(LIBDIR) -lgtest -lgtest_main -lpthread $(OPENCVLIBS)

flow2color: test_Flow2Color.o $(LAOF)/Maths.h $(LAOF)/Flow2Color.o
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS) $(INCDIR)

pyramid: test_GaussianPyramid.o $(LAOF)/GaussianPyramid.o
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

imop: test_ImageOP.o $(LAOF)/ImageProcess.h $(LAOF)/ImageIO.h $(LAOF)/Image.h
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS) $(INCDIR)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(LAOF) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) *.o $(TARGET) *.jpg *.yml
