CXX				= g++
CXXFLAGS  		= -O3 -Wall
INCDIR			= ../include

LINK			= g++
LFLAGS			= 
LIBDIR 			= ../libs
OPENCVLIBS		= `pkg-config --libs opencv`

DEL				= rm -rf

LAOF			= ../src
OBJS			= test_Maths.o test_Flow2Color.o
TARGET 			= run_all_tests flow2color

all: $(TARGET)
	@echo "DONE"

run_all_tests: test_Maths.o $(LAOF)/Maths.o
	$(LINK) $(LFLAGS) -o $@ $^ -L$(LIBDIR) -lgtest -lgtest_main -lpthread $(OPENCVLIBS)

flow2color: test_Flow2Color.o $(LAOF)/Maths.o $(LAOF)/Flow2Color.o
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -I$(LAOF) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) $(OBJS) $(TARGET)
