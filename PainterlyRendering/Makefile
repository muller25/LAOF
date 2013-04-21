CXX				= g++
CXXFLAGS  		= -O3 -Wall 
INCDIR			= `pkg-config --cflags opencv`

LINK			= g++
LFLAGS			= 
LIBDIR 			= .

DEL				= rm -rf

OPENCVLIBS		= `pkg-config --libs opencv`
OBJS			= RenderingImage.o GaussianBlur.o NormalMap.o Lighting.o \
MakeBrushMap.o LightingEffect.o EdgeOptimization.o ImportantMap.o \
PainterlyService.o PainterlyStyle.o ColorModels.o SplineStroke.o CubicBSpline.o \
Queue.o main.o
TARGET			= pr

all: $(TARGET)
	@echo "DONE"

pr: $(OBJS)
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) *.o $(TARGET)
