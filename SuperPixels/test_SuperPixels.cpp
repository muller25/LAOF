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
	
	int numIter = 2;
	if (argc >= 5) numIter = (Value) atoi(argv[4]);
	
	Value lambda = 10;
	if (argc >= 7) lambda = (Value) atoi(argv[6]);
	
	int TYPE = 0;
	if (argc >= 6) TYPE = atoi(argv[5]);
		
	int PATCH_SIZE = atoi(argv[3]);
    Mat I = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
    check_error(I.data == NULL, "Cannot load input image");

	srand(clock());
    int width = I.cols, height = I.rows;
    int num_pixels = width*height; 
	float variance = computeImageVariance(I);
	
	// Initialize and place seeds
	vector<int> Seeds(num_pixels);
	int numSeeds = 0;
	PlaceSeeds(Seeds, numSeeds, I, PATCH_SIZE);

	vector<int> horizWeights(num_pixels,lambda);
	vector<int> vertWeights(num_pixels,lambda);
	vector<int> diag1Weights(num_pixels,lambda);
	vector<int> diag2Weights(num_pixels,lambda);
	if (argc <= 7)
	{
		computeWeights(horizWeights,I,lambda,variance,-1,0,TYPE);
		computeWeights(vertWeights,I,lambda,variance,0,-1,TYPE);
		computeWeights(diag1Weights,I,lambda,variance,-1,-1,TYPE);
		computeWeights(diag2Weights,I,lambda,variance,1,-1,TYPE);
	}
	else{
		loadEdges(horizWeights,I,lambda,argv[7]);
		loadEdges(vertWeights,I,lambda,argv[8]);
		loadEdges(diag1Weights,I,lambda,argv[9]);
		loadEdges(diag2Weights,I,lambda,argv[10]);
	}

	vector<int> labeling(num_pixels);
	initializeLabeling(labeling, Seeds, I, numSeeds, PATCH_SIZE);
		
	Value oldEnergy,newEnergy;
	vector<int> changeMask(num_pixels, 1), changeMaskNew(num_pixels, 0);
    vector<int> order(numSeeds);
	for (int i = 0; i < numSeeds; i++)
		order[i] = i;

	int j = 0;
	//purturbSeeds(order,numSeeds);

    while (1) {
		newEnergy = computeEnergy(I, labeling, horizWeights,vertWeights,
			                      diag1Weights,diag2Weights,Seeds, TYPE);

		if (j == 0){
			oldEnergy = newEnergy+1;
			printf("Initial Energy: %d,  %f sec", newEnergy, ((float)clock())/CLOCKS_PER_SEC);
		}
		//else {
		//	printf("\nAfter iteration %d: ", j-1);
		//	printf("Energy: %d,  %f sec", newEnergy, ((float)clock())/CLOCKS_PER_SEC );
		//}

		if (newEnergy == oldEnergy || j >= numIter) break;

		oldEnergy = newEnergy;
		
		for (int i = 0; i < numSeeds; i++){
			expandOnLabel(order[i], labeling, changeMask, changeMaskNew, I, Seeds,numSeeds,horizWeights,
                          vertWeights, diag1Weights, diag2Weights, lambda,PATCH_SIZE,TYPE,variance);
		}
		for (int i = 0; i < num_pixels; i++){
			changeMask[i] = changeMaskNew[i];
			changeMaskNew[i] = 0;
		}

		//purturbSeeds(order,numSeeds);
		j++;
	}
	
	printf("\nFinal Energy:   %d,",newEnergy);

//	int numSegm = saveSegmentationColor(I,labeling,numSeeds,argv[2]);

    Mat im = imread(argv[1]);
	int numSegm = saveSegmentationEdges(im,labeling,numSeeds,argv[2]);
	printf("%f sec\nNumber of superpixels is %d\n",((float)clock())/CLOCKS_PER_SEC,numSegm);
    return 0;
}
