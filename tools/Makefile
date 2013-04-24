################ Compiler, tools and options

CXX			= g++
CXXFLAGS		= -g
INCPATH			= $(OPENCVINC) 
OPENCVINC		= `pkg-config --cflags opencv`
LINK			= g++
OPENCVLIBS		= `pkg-config --libs opencv`
LIBS			= $(OPENCVLIBS)
DEL_FILE		= rm -f

################ Output directory


############### Files

SOURCES			= $(wildcard *.cpp)
OBJECTS			= $(patsubst %.cpp,%.o,$(SOURCES))
EXECS			= $(patsubst %.cpp,%,$(SOURCES))

############### Build rules

all: $(EXECS) $(OBJECTS)
	@echo "Done!"

$(EXECS):

clean:
	-$(DEL_FILE) *.o
	-$(DEL_FILE) $(EXECS)

################ Compile

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $^

%: %.o
	$(LINK) -o $@ $^ $(LIBS)
