CXX				= g++
CXXFLAGS  		= -O3 -Wall 
INCDIR			= `pkg-config --cflags opencv`

LINK			= g++
LFLAGS			= 
LIBDIR 			= .

DEL				= rm -rf

OPENCVLIBS		= `pkg-config --libs opencv`
OBJS			= NormalMap.o Lighting.o MakeBrushMap.o LightingEffect.o \
EdgeOptimization.o ImportantMap.o PainterlyService.o PainterlyStyle.o \
ColorModels.o SplineStroke.o
TARGET			= pr vr

all: $(TARGET)
	@echo "DONE"

vr: $(OBJS) vr.o
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

pr: $(OBJS) pr.o
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) *.o $(TARGET) *.jpg *.out
