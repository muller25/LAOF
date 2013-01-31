CXX				= g++
CFLAGS  		= -O3
INCPATH			= -I$(LIBDIR)
LINK			= g++
LFLAGS			= 
OPENCVLIBS		= `pkg-config --libs opencv`
LIBS			= $(OPENCVLIBS) -L$(LIBDIR)
DEL				= rm -rf

LIBDIR 			= .
OBJS			= Maths.o toy.o
TARGET 			= toy

all: $(TARGET)
	@echo "Done"

$(TARGET): Maths.o toy.o
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

clean:
	$(DEL) $(OBJS) $(TARGET)

%.o: %.cpp
	$(CXX) -c $(CFLAGS) $(INCPATH) -o $@ $^