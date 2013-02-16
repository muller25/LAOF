CXX				= g++
CFLAGS  		= -O3
INCPATH			= -I$(LIBDIR)
LINK			= g++
LFLAGS			= 
OPENCVLIBS		= `pkg-config --libs opencv`
LIBS			= $(OPENCVLIBS) -L$(LIBDIR)
DEL				= rm -rf

LIBDIR 			= .
OBJS			= Maths.o toy.o OpticalFlow.o
TARGET 			= toy

all: $(TARGET)
	@echo "Done"

$(TARGET): $(OBJS)
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

toy.o: toy.cpp
	$(CXX) $(CFLAGS) -c toy.cpp -o toy.o

%.o: %.cpp %.h
	$(CXX) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) $(OBJS) $(TARGET)
