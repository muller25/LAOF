#include "SuperPixels.h"

#include <vector>
#include <cstdio>
using namespace std;

#include <cv.h>
#include <highgui.h>
using namespace cv;

void check_error(int boolE, const char *error_message)
{ 
    if (boolE) 
    {
        printf("%s\n", error_message);
        exit(1);
    }
}

void check_input_arguments(int argc)
{
	check_error((argc < 4), "superpixels input.pgm output.ppm patch_size [NUM_ITERATIONS] [TYPE] [lambda] [h_edges.pgm v_edges.pgm d_b_edges.pgm d_f_edges.pgm]");
	check_error((argc > 7 && argc < 11), "If providing edge images, provide all 4 of them");
}

int main(int argc, char **argv)
{
	check_input_arguments(argc);

    Mat im = imread(argv[1]);
    check_error(im.data == NULL, "Cannot load input image");

	int patch_size = atoi(argv[3]);
	
	int numIter = 2;
	if (argc >= 5) numIter = (Value) atoi(argv[4]);

	int TYPE = 0;
	if (argc >= 6) TYPE = atoi(argv[5]);
	
	Value lambda = 10;
	if (argc >= 7) lambda = (Value) atoi(argv[6]);
	
    SuperPixels sp(patch_size, numIter, TYPE, lambda);
    sp.setSourceImage(im);
	if (argc > 7) {
		sp.loadHEdge(argv[7]);
		sp.loadVEdge(argv[8]);
		sp.loadD1Edge(argv[9]);
		sp.loadD2Edge(argv[10]);
	}

    clock_t start = clock();
    sp.generate();
//	int numSegm = saveSegmentationColor(argv[2]);
    clock_t end = clock();
    double duration = (double)(end - start) / CLOCKS_PER_SEC;

	int numSegm = sp.saveSegmentationEdges(argv[2]);
	printf("total used: %f sec\nNumber of superpixels is %d\n", duration, numSegm);
    return 0;
}
