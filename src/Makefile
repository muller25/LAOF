CXX				= g++
CXXFLAGS  		= -O3 -Wall 
INCDIR			= `pkg-config --cflags opencv`

LINK			= g++
LFLAGS			= 
LIBDIR 			= .

DEL				= rm -rf

GC				= ../gc
SP				= ../SuperPixels
OPENCVLIBS		= `pkg-config --libs opencv`
OBJS			= OpticalFlow.o GaussianPyramid.o Flow2Color.o MotionLayers.o LAOF.o OFPara.o
TARGET			= biof adof stof #layers em

all: $(TARGET) MotionLayers.o
	@echo "DONE"

biof: main_biof.o OpticalFlow.o GaussianPyramid.o Flow2Color.o
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

adof: main_adof.o OpticalFlow.o GaussianPyramid.o Flow2Color.o
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

stof: main_stof.o OpticalFlow.o GaussianPyramid.o Flow2Color.o
	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

#em: main_em.o $(OBJS) $(GC)/GCoptimization.o $(GC)/LinkedBlockList.o
#	$(LINK) $(LFAGS) -o $@ $^ $(OPENCVLIBS)

#layers: main_motionlayers.o MotionLayers.o Flow2Color.o $(GC)/GCoptimization.o $(GC)/LinkedBlockList.o
#	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

#main_motionlayers.o: main_motionlayers.cpp ImageProcess.h Utils.h
#	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

OpticalFlow.o: OpticalFlow.cpp OpticalFlow.h Maths.h Image.h GaussianPyramid.h ImageProcess.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

Flow2Color.o: Flow2Color.cpp Flow2Color.h Image.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

MotionLayers.o: MotionLayers.cpp MotionLayers.h ML.h types.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(GC) -I$(SP) -c $< -o $@

#LAOF.o: LAOF.cpp LAOF.h Image.h OpticalFlow.h MotionLayers.h
#	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

OFPara.o: OFPara.cpp OFPara.h Maths.h ImageProcess.h Image.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) *.o $(TARGET) cscope.out *.yml *.jpg
