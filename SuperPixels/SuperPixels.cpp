////////////////////////////////////////////////////////
////////// First version of faster iterations///////////
////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <image.h>
#include <misc.h>
#include <matrix.h>
#include <misc.h>
#include <draw.h>
#include <errno.h>

#include "SuperPixels.h"
#include "graph.h"
#include "energy.h"

#include <iostream>
#include <fstream>
using namespace std;
using namespace vlib;

#include <cv.h>
#include <highgui.h>
using namespace cv;

#define NUM_COLORS 255
#define MULTIPLIER_VAR 1.5
#define sq(x) ((x)*(x))

void PlaceSeeds(vector<int> &Seeds, int &numSeeds, const Mat &I, int PATCH_SIZE)
{
    assert(I.type() == CV_8U && numSeeds >= 0);

    int width = I.cols, height = I.rows;
	int bSize = PATCH_SIZE / 2, delta = PATCH_SIZE / 4 - 1;
    
    for (int h = bSize; h < height; h += bSize)
        for (int w = bSize; w < width; w += bSize){
            
            int bestX = w, bestY = h;
            int startX = w - delta, endX = w + delta;
            int startY = h - delta, endY = h + delta;

            if (startX < 1) startX = 1;
            if (endX > width-2) endX = width-2;
            if (startY < 1) startY = 1;
            if (endY > height-2) endY = height-2;

            int bestScore = 255 * sq(2*delta+1);
            //int bestScore = 0;

            for (int y = startY; y <= endY; y++)
                for (int x = startX; x <= endX; x++){
                    int currScore = 0;

                    for (int yy = -1; yy <= 1; ++yy)
                        for (int xx = -1; xx <= 1; ++xx)
                        {
                            int c1 = I.at<uchar>(y, x);
                            int c2 = I.at<uchar>(y+yy, x+xx);
                            currScore += abs(c1 - c2);
                        }

                    if (currScore < bestScore){
                        bestScore = currScore;
                        bestX = x;
                        bestY = y;
                    }
                }
            Seeds[numSeeds++] = bestX + bestY * width;
        }
}

Value computeEnergy(const Mat &I, const vector<int> &labeling,
					const vector<Value> &horizWeights, const vector<Value> &vertWeights,
					const vector<Value> &diag1Weights, const vector<Value> &diag2Weights,
					const vector<int> &Seeds, int TYPE)
{
    assert(I.type() == CV_8U);

    int height = I.rows, width = I.cols;
	TotalValue engSmooth = 0, engData = 0;
	float SQRT_2= 1/sqrt(2.0);

	if (TYPE == 1)
	{
		for (int y = 0; y < height; y++)
			for (int x = 0; x < width; x++){
				int label = labeling[x+y*width];
				int seedY = Seeds[label]/width;
				int seedX = Seeds[label] - seedY*width;
				int scolor = I.at<uchar>(seedY, seedX);
                int color = I.at<uchar>(y, x);
				int diff = abs(color - scolor);
				int maxD = 15;
				if (diff > maxD) diff = maxD;
				engData = engData + diff;
			}
	}

	for (int y = 1; y < height; y++){
		for (int x = 0; x < width; x++)
			if (labeling[x+y*width] != labeling[x+(y-1)*width])
				engSmooth = engSmooth + vertWeights[x+(y-1)*width];
    }

	for (int y = 0; y < height; y++){
		for (int x = 1; x < width; x++)
			if (labeling[x+y*width] != labeling[(x-1)+y*width]){
				engSmooth = engSmooth + horizWeights[(x-1)+y*width];
			}
    }

	for (int y = 1; y < height; y++){
		for (int x = 1; x < width; x++)
			if (labeling[x+y*width] != labeling[x-1+(y-1)*width])
				engSmooth = engSmooth + SQRT_2*diag1Weights[x-1+(y-1)*width];
    }

	for (int y = 1; y < height; y++){
		for (int x = 0; x < width-1; x++)
			if (labeling[x+y*width] != labeling[x+1+(y-1)*width])
				engSmooth = engSmooth + SQRT_2*diag2Weights[x+1+(y-1)*width];
	}

	//printf("\nDeng %d ",engSmooth);
	return (engSmooth+engData);
}

void getBounds(int &seedX, int &seedY, int &startX, int &startY, int &endX, int &endY,
               const Mat &I, const vector<int> &Seeds, int label, int PATCH_SIZE)
{
    int width = I.cols, height = I.rows;
    
	seedY = Seeds[label]/width;
	seedX = Seeds[label]-seedY*width;

	startX = seedX-PATCH_SIZE;
	endX   = seedX+PATCH_SIZE;

	startY = seedY-PATCH_SIZE;
	endY   = seedY+PATCH_SIZE;

	if (startX < 0) startX = 0;
	if (startY < 0) startY = 0;
	if (endX >= width) endX = width-1;
	if (endY >= height) endY = height-1;
}

// variable 0 corresponds to the old label, variable 1 to the new label, which
// is the first input parameter into the procedure
// does alpha-expansion in a block of size BLOCK_SIZE by BLOCK_SIZE
// the border stays fixed to the old label values
void expandOnLabel(int label, vector<int> &labeling,
                   vector<int> &changeMask, vector<int> &changeMaskNew,
                   const Mat &I, const vector<int> &Seeds, int numSeeds, 
				   const vector<Value> &horizWeights, const vector<Value> &vertWeights, 
				   const vector<Value> &diag1Weights, const vector<Value> &diag2Weights,
                   Value lambda, int PATCH_SIZE, int TYPE, float variance)
{
	int seedX,seedY,startX,startY,endX,endY,numVars,blockWidth;
	int width = I.cols, height = I.rows, somethingChanged = 0;

	getBounds(seedX,seedY,startX,startY,endX,endY,I,Seeds,label,PATCH_SIZE);

	for (int y = startY; y <= endY; y++)
		for (int x = startX; x <= endX; x++)
			if (changeMask[x+y*width] == 1)
			{
				somethingChanged = 1;
				break;
			}
	
	if (somethingChanged == 0) return;
		
	blockWidth = endX-startX+1;
	numVars = (endY-startY+1)*blockWidth;
	
	vector<Var> variables(numVars);
	Energy<int,int,int> *e = new Energy<int,int,int>(numVars,numVars*3);

	for (int i = 0; i < numVars; i++)
        variables[i] = e->add_variable();

	Value LARGE_WEIGHT = lambda*NUM_COLORS*8;

	// First fix the border to old labels, except the edges of the image
	for (int y = startY; y <= endY; y++){
		if (startX != 0)
			e->add_term1(variables[(y-startY)*blockWidth],0,LARGE_WEIGHT);
		else if (y == startY || y == endY)
			e->add_term1(variables[(y-startY)*blockWidth],0,LARGE_WEIGHT);

		if(endX != width - 1)
			e->add_term1(variables[(endX-startX)+(y-startY)*blockWidth],0,LARGE_WEIGHT);
		else if (y == startY || y == endY)
			e->add_term1(variables[(endX-startX)+(y-startY)*blockWidth],0,LARGE_WEIGHT);
	}

	for (int x = startX+1; x < endX; x++){
		if (startY != 0)
			e->add_term1(variables[(x-startX)],0,LARGE_WEIGHT);
		if (endY != height - 1)
			e->add_term1(variables[(x-startX)+(endY-startY)*blockWidth],0,LARGE_WEIGHT);
	}

	// add links to center of the patch for color constant superpixels
	if (TYPE == 1)
	{
		int scolor = I.at<uchar>(seedY, seedX);
		for (int y = startY+1; y < endY; y++)
			for (int x = startX+1; x < endX; x++){
				Value E00=0,E01=0,E10=LARGE_WEIGHT,E11=0;

				if (seedX != x && seedY != y)
					e->add_term2(variables[(x-startX)+(y-startY)*blockWidth],
                                 variables[(seedX-startX)+(seedY-startY)*blockWidth],E00,E01,E10,E11);

                int color = I.at<uchar>(y, x);
                int diff = abs(color - scolor);
				int maxD = (int) variance * MULTIPLIER_VAR;
				if (diff > maxD) diff = maxD;

				int oldLabel = labeling[x+y*width];
				int oldY = Seeds[oldLabel]/width;
				int oldX = Seeds[oldLabel]-oldY*width;
				int oldColor = I.at<uchar>(oldY, oldX);
				int oldDiff = abs(color - oldColor);
				if (oldDiff > maxD) oldDiff = maxD;

				if (oldDiff > diff)
                    e->add_term1(variables[(x-startX)+(y-startY)*blockWidth],oldDiff-diff,0);
				else
                    e->add_term1(variables[(x-startX)+(y-startY)*blockWidth],0,diff-oldDiff);
			}
	}

	// First set up horizontal links 
	for (int y = startY; y <= endY; y++)
		for (int x = startX+1; x <= endX; x++){
			int oldLabelPix       = labeling[x+y*width];
			int oldLabelNeighbPix = labeling[x-1+y*width];
			Value E00,E01,E10,E11 = 0;

			if (oldLabelPix != oldLabelNeighbPix)
                E00 = horizWeights[x-1+y*width];
			else
                E00 = 0;
			if (oldLabelNeighbPix != label)
                E01 = horizWeights[x-1+y*width];
			else
                E01 = 0;
			if (label != oldLabelPix)
                E10 = horizWeights[x-1+y*width];
			else
                E10 = 0;
			
			e->add_term2(variables[(x-startX)-1+(y-startY)*blockWidth],variables[(x-startX)+(y-startY)*blockWidth],E00,E01,E10,E11);
		}

	// Next set up vertical links
	for (int y = startY+1; y <= endY; y++)
		for (int x = startX; x <=endX; x++){
			int oldLabelPix       = labeling[x+y*width];
			int oldLabelNeighbPix = labeling[x+(y-1)*width];
			Value E00,E01,E10,E11=0;

			if (oldLabelPix != oldLabelNeighbPix)
                E00 = vertWeights[x+(y-1)*width];
			else
                E00 = 0;
			if (oldLabelNeighbPix != label)
                E01 = vertWeights[x+(y-1)*width];
			else
                E01 = 0;
			if (label != oldLabelPix)
                E10 = vertWeights[x+(y-1)*width];
			else
                E10 = 0;
			
			e->add_term2(variables[(x-startX)+(y-startY-1)*blockWidth],variables[(x-startX)+(y-startY)*blockWidth],E00,E01,E10,E11);
		}

	// Next set up diagonal links 
	float SQRT_2= 1/sqrt(2.0);
	for (int y = startY+1; y <= endY; y++)
		for (int x = startX+1; x <= endX; x++){
			int oldLabelPix       = labeling[x+y*width];
			int oldLabelNeighbPix = labeling[x-1+(y-1)*width];
			Value E00,E01,E10,E11=0;
			
			if (oldLabelPix != oldLabelNeighbPix) 
				E00 = SQRT_2*diag1Weights[x-1+(y-1)*width];
			else E00 = 0;
			if (oldLabelNeighbPix != label) 
				E01 = SQRT_2*diag1Weights[x-1+(y-1)*width];
			else E01 = 0;
			if (label != oldLabelPix)
				E10 = SQRT_2*diag1Weights[x-1+(y-1)*width];
			else E10 = 0;
			
			e->add_term2(variables[(x-startX)-1+(y-startY-1)*blockWidth],variables[(x-startX)+(y-startY)*blockWidth],E00,E01,E10,E11);
		}
	
	// More diagonal links
	for (int y = startY+1; y <= endY; y++)
		for (int x = startX; x <=endX-1; x++){
			int oldLabelPix       = labeling[x+y*width];
			int oldLabelNeighbPix = labeling[(x+1)+(y-1)*width];
			Value E00,E01,E10,E11=0;

			if (oldLabelPix != oldLabelNeighbPix) 
				E00 = SQRT_2*diag2Weights[(x+1)+(y-1)*width];
			else
                E00 = 0;
			if (oldLabelNeighbPix != label) 
				E01 = SQRT_2*diag2Weights[(x+1)+(y-1)*width];
			else
                E01 = 0;
			if (label != oldLabelPix)
				E10 = SQRT_2*diag2Weights[(x+1)+(y-1)*width];
			else
                E10 = 0;

			e->add_term2(variables[(x-startX+1)+(y-startY-1)*blockWidth],variables[(x-startX)+(y-startY)*blockWidth],E00,E01,E10,E11);
		}

	e->minimize();

	for (int y = startY; y <= endY; y++)
		for (int x = startX; x <= endX; x++){
			if (e->get_var(variables[(x-startX)+(y-startY)*blockWidth]) != 0)
			{
				if (labeling[x+y*width] != label){
					labeling[x+y*width] = label; 
					changeMaskNew[x+y*width] = 1;
					changeMask[x+y*width] = 1;
				}
			}
		}

	delete e;
}

void initializeLabeling(vector<int> &labeling, const vector<int> &Seeds, const Mat &I,
                        int numSeeds, int PATCH_SIZE)
{
    int width = I.cols, height = I.rows;
    int seedX, seedY, startX, startY, endX, endY;
	for (int i = 0; i < numSeeds; i++)
    {
		seedY = Seeds[i]/width;
		seedX = Seeds[i]-seedY*width;

		startX = seedX-PATCH_SIZE/2-1;
		endX   = seedX+PATCH_SIZE/2+1;

		startY = seedY-PATCH_SIZE/2-1;
		endY   = seedY+PATCH_SIZE/2+1;

		if (startX < 0)     startX  = 0;
		if (startY < 0)     startY  = 0;
		if (endX >= width)  endX    = width-1;
		if (endY >= height) endY    = height-1;

		for (int y = startY; y <= endY; y++)
			for (int x = startX; x <= endX; x++)
				labeling[x+y*width] = i;
	}
}

float computeImageVariance(const Mat &I)
{
    assert(I.type() == CV_8U);
    
    int width = I.cols, height = I.rows, total = 0;
	float v = (float) 0.0;

	for (int y = 1; y < height; y++)
		for (int x = 1; x < width; x++){
            int off = I.at<uchar>(y, x);
            int offx_1 = I.at<uchar>(y, x-1);
            int offy_1 = I.at<uchar>(y-1, x);

            v += abs(off - offx_1) + abs(off - offy_1);
			total += 2;
        }
    
	return (v/total);
}

void purturbSeeds(vector<int> &order,int numSeeds)
{
	for (int i = 0; i < 3*numSeeds; i++ )
	{
		int first  = (rand()*rand())%numSeeds;
		int second = (rand()*rand())%numSeeds;
		int temp = order[first];
		order[first] = order[second];
		order[second] = temp;
	}
}

void computeWeights(vector<Value> &weights, const Mat &I, Value lambda, float variance,
					int incrX, int incrY, int TYPE)
{
    assert(I.type() == CV_8U);

    int height = I.rows, width = I.cols, startX = 0, startY = 0;
	float sigma = 2.0f;

	if (incrX != 0) startX  = abs(incrX);
	if (incrY != 0) startY = abs(incrY);

	Value smallPenalty;
	if (TYPE == 1) smallPenalty = (MULTIPLIER_VAR*variance)/8+1;
	else smallPenalty = 1;

	for (int y = startY; y < height; y++)
		for (int x = startX; x < width; x++){
            int c1 = I.at<uchar>(y, x), c2 = I.at<uchar>(y+incrY, x+incrX);
			int difference = sq((c1 - c2));
            weights[(x+incrX)+(y+incrY)*width] = (Value) (lambda*exp((-difference/(sigma*sq(variance))))+smallPenalty);
		}
	
	//image<uchar> *e = new image<uchar>(width,height);
		
	//for ( int y = 0; y < height; y++ )
	//	for (int  x = 0; x < width; x++ ){
	//		imRef(e,x,y) = weights[x+y*width];
	//	}

	//savePGM(e,name);  
}

void loadEdges(vector<Value> &weights, const Mat &I, Value lambda, const char *name)
{
    int height = I.rows, width = I.cols;
	Mat edges = imread(name, CV_LOAD_IMAGE_GRAYSCALE);
    assert(height == edges.rows && width == edges.cols);
    
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++){
            int c = edges.at<uchar>(y, x);
			weights[x+y*width] = (Value) lambda * c;
		}
}

int saveSegmentationColor(const Mat &I, const vector<int> &labeling, int numSeeds, const char *name)
{
    int width = I.cols, height = I.rows;
    Mat out(height, width, CV_8UC3);

	vector<Vec3b> colorLookup(numSeeds);
	vector<int> counts(numSeeds,0);

    // Vec3b: BGR
	for (int i = 0; i < numSeeds; i++){
		colorLookup[i][2] = rand()%NUM_COLORS;
		colorLookup[i][1] = rand()%NUM_COLORS;
		colorLookup[i][0] = rand()%NUM_COLORS;
	}

	for (int y = 0; y < height; y++)
		for (int  x = 0; x < width; x++ ){
			out.at<Vec3b>(y, x) = colorLookup[labeling[x+y*width]];
			counts[labeling[x+y*width]]++;
		}

    imwrite(name, out);

/*    
      image<uchar> *e = new image<uchar>(width,height);
      e->init(0);
      for ( int y = 1; y < height; y++ )
      for (int  x = 1; x < width; x++ ){
      if ( labeling[x+y*width] != labeling[x-1+y*width])
      imRef(e,x-1,y) = 255;
      if ( labeling[x+y*width] != labeling[x+(y-1)*width])
      imRef(e,x,y-1) = 255;
      }

      savePGM(e,"edges.pgm");  
*/
	int num_superpixels = 0;
	for (int i = 0; i < numSeeds; i++)
		if (counts[i] > 0 ) num_superpixels++;

	return (num_superpixels);
}

int saveSegmentationEdges(const Mat &I, const vector<int> &labeling, int numSeeds, const char *name)
{
    assert(I.type() == CV_8UC3);

 	vector<int> counts(numSeeds,0);
    int height = I.rows, width = I.cols;
    Mat dst;
    I.copyTo(dst);

    Vec3b green(0, 255, 0);
    for (int h = 1; h < height; ++h)
        for (int w = 1; w < width; ++w)
        {
            int off = w + h * width;
            int offx_1 = (w-1) + h * width;
            int offy_1 = w + (h-1) * width;
            if (labeling[off] != labeling[offx_1]) {
                dst.at<Vec3b>(h, w) = green;
                dst.at<Vec3b>(h, w-1) = green;
            }
            if (labeling[off] != labeling[offy_1]) {
                dst.at<Vec3b>(h, w) = green;
                dst.at<Vec3b>(h-1, w) = green;
            }
 			counts[labeling[off]]++;
        }

    imwrite(name, dst);

	int num_superpixels = 0;
    for (int i = 0; i < numSeeds; ++i)
        if (counts[i] > 0) num_superpixels++;

    return num_superpixels;
}
