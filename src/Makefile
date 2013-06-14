CXX				= g++
CXXFLAGS  		= -O3 -Wall 
INCDIR			= `pkg-config --cflags opencv`

LINK			= g++
LFLAGS			= 
LIBDIR 			= .

DEL				= rm -rf

GC				= ../gc
OPENCVLIBS		= `pkg-config --libs opencv`
OBJS			= OpticalFlow.o GaussianPyramid.o Flow2Color.o MotionLayers.o LAOF.o OFPara.o
TARGET			= biof adof stof superpixels layers #em

all: $(TARGET)
	@echo "DONE"

OpticalFlow.o: OpticalFlow.cpp OpticalFlow.h Maths.h Image.h GaussianPyramid.h ImageProcess.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

Flow2Color.o: Flow2Color.cpp Flow2Color.h Image.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

MotionLayers.o: MotionLayers.cpp MotionLayers.h ML.h types.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(GC) -c $< -o $@

#LAOF.o: LAOF.cpp LAOF.h Image.h OpticalFlow.h MotionLayers.h
#	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

OFPara.o: OFPara.cpp OFPara.h Maths.h ImageProcess.h Image.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

graph.o: graph.cpp graph.h
	$(CXX) $(CFLAGS) -c graph.cpp -o graph.o

maxflow.o: maxflow.cpp 
	$(CXX) $(CFLAGS) -c maxflow.cpp -o maxflow.o

SuperPixels.o: SuperPixels.cpp SuperPixels.h
	$(CXX) $(CFLAGS) $(OPENCVFLAGS) -I$(GC) -c $< -o $@

biof: main_biof.o OpticalFlow.o GaussianPyramid.o Flow2Color.o
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

adof: main_adof.o OpticalFlow.o GaussianPyramid.o Flow2Color.o
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

stof: main_stof.o OpticalFlow.o GaussianPyramid.o Flow2Color.o
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

#em: main_em.o $(OBJS) $(GC)/GCoptimization.o $(GC)/LinkedBlockList.o
#	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

layers: main_motionlayers.o MotionLayers.o SuperPixels.o graph.o maxflow.o $(GC)/GCoptimization.o $(GC)/LinkedBlockList.o
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

superpixels: graph.o maxflow.o SuperPixels.o main_superpixels.o
	$(CXX) -o superpixels $^ $(OPENCVLIBS)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) *.o $(TARGET) cscope.out *.yml *.jpg
