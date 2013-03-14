CXX				= g++
CXXFLAGS  		= -O3 -Wall
INCDIR			= -I../include `pkg-config --cflags opencv`

LINK			= g++
LFLAGS			= 
LIBDIR 			= ../libs
OPENCVLIBS		= `pkg-config --libs opencv`

DEL				= rm -rf

GC				= ../gc
LAOF			= ../src
TARGET 			= flow2color pyramid imop kmeans layers tests gc seg

all: $(TARGET)
	@echo "DONE"

test_Image.o: test_Image.cpp $(LAOF)/Image.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(LAOF) -c $< -o $@

test_Maths.o: test_Maths.cpp $(LAOF)/Maths.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(LAOF) -c $< -o $@

test_hist.o: test_hist.cpp $(LAOF)/Histogram.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(LAOF) -c $< -o $@

tests: test_Maths.o test_Image.o test_hist.o
	$(LINK) $(LFLAGS) -o $@ $^ -L$(LIBDIR) -lgtest -lgtest_main -lpthread

flow2color: test_Flow2Color.o $(LAOF)/Flow2Color.o
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

pyramid: test_GaussianPyramid.o $(LAOF)/GaussianPyramid.o
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

test_ImageOP.o: test_ImageOP.cpp $(LAOF)/ImageProcess.h $(LAOF)/ImageIO.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(LAOF) -c $< -o $@

imop: test_ImageOP.o 
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

test_kmeans.o: test_kmeans.cpp $(LAOF)/ML.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(LAOF) -c $< -o $@

kmeans: test_kmeans.o 
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

test_MotionLayers.o: test_MotionLayers.cpp $(LAOF)/MotionLayers.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(LAOF) -c $< -o $@

layers: test_MotionLayers.o $(LAOF)/MotionLayers.o $(LAOF)/Flow2Color.o
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

test_gc.o: test_gc.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(LAOF) -I$(GC) -c $< -o $@

gc: test_gc.o $(GC)/GCoptimization.o $(GC)/LinkedBlockList.o
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

test_segment.o: test_segment.cpp $(LAOF)/SegmentImage.h $(LAOF)/SegmentGraph.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(LAOF) -c $< -o $@

seg: test_segment.o $(LAOF)/Flow2Color.o
	$(LINK) $(LFLAGS) -o $@ $^ $(OPENCVLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) -I$(LAOF) -c $< -o $@

.PHONY: clean
clean:
	-$(DEL) *.o $(TARGET) *.jpg *.yml cscope.out
