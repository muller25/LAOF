CXX				= g++
CXXFLAGS  		= -O3 -Wall -fpermissive

LINK			= g++
LFLAGS			= 
LIBDIR 			= .

DEL				= rm -rf

OBJS			= example.o GCoptimization.o graph.o LinkedBlockList.o maxflow.o

TARGET			= gc

all: $(TARGET)
	@echo "DONE"

gc: $(OBJS)
	$(LINK) $(LFAGS) -o gc $(OBJS) 

example.o: example.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) $(OBJS) $(TARGET) cscope.out
