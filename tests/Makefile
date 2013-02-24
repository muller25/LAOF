CXX				= g++
CXXFLAGS  		= -O3 -Wall
INCDIR			= ../include

LINK			= g++
LFLAGS			= 
LIBDIR 			= ../libs
OPENCVLIBS		= `pkg-config --libs opencv`

DEL				= rm -rf

LAOF			= ../src
OBJS			= test_Maths.o
TARGET 			= run_all_tests

all: $(TARGET)
	@echo "DONE"

run_all_tests: $(OBJS)
	$(LINK) $(LFLAGS) -o $@ $(OBJS) -L$(LIBDIR) -lgtest -lgtest_main -lpthread -L$(LAOF) $(OPENCVLIBS)

test_Maths.o: test_Maths.cpp $(LAOF)/Maths.o
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -I$(LAOF) -c $< -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) $(OBJS) $(TARGET)
