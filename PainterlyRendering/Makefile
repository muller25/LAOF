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
ColorModels.o SplineStroke.o RBF.o StupidPaint.o
TARGET			= pr vr rbf stupid

all: $(TARGET)
	@echo "DONE"

rbf: RBF.o test_rbf.o
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

stupid: StupidPaint.o Brush.h test_StupidPaint.o RBF.o
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

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
