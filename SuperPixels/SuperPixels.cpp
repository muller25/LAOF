////////////////////////////////////////////////////////
////////// First version of faster iterations///////////
////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "SuperPixels.h"
#include "graph.h"
#include "energy.h"

#include <iostream>
using namespace std;

#include <cv.h>
#include <highgui.h>
using namespace cv;

#define NUM_COLORS 255
#define MULTIPLIER_VAR 1.5
#define sq(x) ((x)*(x))
#define DEBUG

#include "types.h"

void SuperPixels::setSourceImage(Mat &im)
{
    im.copyTo(m_im);
    m_width = m_im.cols, m_height = m_im.rows;
    if (m_im.channels() == 1) m_im.copyTo(m_gray);
    else cvtColor(m_im, m_gray, CV_BGR2GRAY);

    m_label.create(m_height, m_width, CV_LABEL);
    m_seeds.clear();
    m_init = true;
 	m_variance = computeImageVariance();
    computeWeights();
}

void SuperPixels::PlaceSeeds()
{
    assert(m_init);
    m_seeds.clear();
    
	int bSize = m_patch_size / 2, delta = m_patch_size / 4 - 1;
    for (int h = bSize; h < m_height; h += bSize)
        for (int w = bSize; w < m_width; w += bSize)
        {
            int bestX = w, bestY = h;
            int startX = w - delta, endX = w + delta;
            int startY = h - delta, endY = h + delta;
            int bestScore = 255 * sq(2*delta+1);
            
            if (startX < 1) startX = 1;
            if (endX > m_width-2) endX = m_width-2;
            if (startY < 1) startY = 1;
            if (endY > m_height-2) endY = m_height-2;

            for (int y = startY; y <= endY; y++)
                for (int x = startX; x <= endX; x++){
                    int currScore = 0;

                    for (int yy = -1; yy <= 1; ++yy)
                        for (int xx = -1; xx <= 1; ++xx)
                        {
                            int c1 = m_gray.at<PERCISION>(y, x);
                            int c2 = m_gray.at<PERCISION>(y+yy, x+xx);
                            currScore += abs(c1 - c2);
                        }

                    if (currScore < bestScore){
                        bestScore = currScore;
                        bestX = x, bestY = y;
                    }
                }
            m_seeds.push_back(Point(bestX, bestY));
        }
}

Value SuperPixels::computeEnergy()
{
    assert(m_init);

	TotalValue engSmooth = 0, engData = 0;
	float SQRT_2 = 1./sqrt(2.0);

	if (m_TYPE == 1)
	{
		for (int y = 0; y < m_height; y++)
			for (int x = 0; x < m_width; x++){
				int label = m_label.at<LABEL>(y, x);
				int seedY = m_seeds[label].y;
				int seedX = m_seeds[label].x;
				int scolor = m_gray.at<PERCISION>(seedY, seedX);
                int color = m_gray.at<PERCISION>(y, x);
				int diff = abs(color - scolor);
				int maxD = 15;
				if (diff > maxD) diff = maxD;
				engData += diff;
			}
	}

	for (int y = 1; y < m_height; y++)
		for (int x = 0; x < m_width; x++)
			if (m_label.at<LABEL>(y, x) != m_label.at<LABEL>((y-1), x))
				engSmooth += m_vertWeights.at<Value>(y-1, x);

	for (int y = 0; y < m_height; y++)
		for (int x = 1; x < m_width; x++)
            if (m_label.at<LABEL>(y, x) != m_label.at<LABEL>(y, (x-1)))
				engSmooth += m_horizWeights.at<Value>(y, x-1);

	for (int y = 1; y < m_height; y++)
		for (int x = 1; x < m_width; x++)
            if (m_label.at<LABEL>(y, x) != m_label.at<LABEL>(y-1, x-1))
				engSmooth += SQRT_2 * m_diag1Weights.at<Value>(y-1, x-1);

	for (int y = 1; y < m_height; y++)
		for (int x = 0; x < m_width-1; x++)
            if (m_label.at<LABEL>(y, x) != m_label.at<LABEL>(y-1, x+1))
				engSmooth += SQRT_2 * m_diag2Weights.at<Value>(y-1, x+1);

	return (engSmooth+engData);
}

void SuperPixels::getBounds(int &seedX, int &seedY, int &startX, int &startY,
                            int &endX, int &endY, int label)
{
    assert(m_init);
    
	seedY = m_seeds[label].y;
	seedX = m_seeds[label].x;

	startX = seedX - m_patch_size;
	endX   = seedX + m_patch_size;

	startY = seedY - m_patch_size;
	endY   = seedY + m_patch_size;

	if (startX < 0) startX = 0;
	if (startY < 0) startY = 0;
	if (endX >= m_width) endX = m_width-1;
	if (endY >= m_height) endY = m_height-1;
}

// variable 0 corresponds to the old label, variable 1 to the new label, which
// is the first input parameter into the procedure
// does alpha-expansion in a block of size BLOCK_SIZE by BLOCK_SIZE
// the border stays fixed to the old label values
void SuperPixels::expandOnLabel(int label,
                                vector<int> &changeMask, vector<int> &changeMaskNew)
{
    assert(m_init);
    
	int seedX,seedY,startX,startY,endX,endY,numVars,blockWidth;
	int somethingChanged = 0;

	getBounds(seedX,seedY,startX,startY,endX,endY,label);

	for (int y = startY; y <= endY; y++)
		for (int x = startX; x <= endX; x++)
			if (changeMask[x+y*m_width] == 1)
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

	Value LARGE_WEIGHT = m_lambda*NUM_COLORS*8;

	// First fix the border to old labels, except the edges of the image
	for (int y = startY; y <= endY; y++){
		if (startX != 0)
			e->add_term1(variables[(y-startY)*blockWidth],0,LARGE_WEIGHT);
		else if (y == startY || y == endY)
			e->add_term1(variables[(y-startY)*blockWidth],0,LARGE_WEIGHT);

		if(endX != m_width - 1)
			e->add_term1(variables[(endX-startX)+(y-startY)*blockWidth],0,LARGE_WEIGHT);
		else if (y == startY || y == endY)
			e->add_term1(variables[(endX-startX)+(y-startY)*blockWidth],0,LARGE_WEIGHT);
	}

	for (int x = startX+1; x < endX; x++){
		if (startY != 0)
			e->add_term1(variables[(x-startX)],0,LARGE_WEIGHT);
		if (endY != m_height - 1)
			e->add_term1(variables[(x-startX)+(endY-startY)*blockWidth],0,LARGE_WEIGHT);
	}

	// add links to center of the patch for color constant superpixels
	if (m_TYPE == 1)
	{
		int scolor = m_gray.at<uchar>(seedY, seedX);
		for (int y = startY+1; y < endY; y++)
			for (int x = startX+1; x < endX; x++){
				Value E00=0,E01=0,E10=LARGE_WEIGHT,E11=0;

				if (seedX != x && seedY != y)
					e->add_term2(variables[(x-startX)+(y-startY)*blockWidth],
                                 variables[(seedX-startX)+(seedY-startY)*blockWidth],E00,E01,E10,E11);

                int color = m_gray.at<uchar>(y, x);
                int diff = abs(color - scolor);
				int maxD = (int) m_variance * MULTIPLIER_VAR;
				if (diff > maxD) diff = maxD;

				int oldLabel = m_label.at<int>(y, x);
                int oldY = m_seeds[oldLabel].y;
				int oldX = m_seeds[oldLabel].x;
				int oldColor = m_gray.at<uchar>(oldY, oldX);
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
			int oldLabelPix       = m_label.at<int>(y, x);
			int oldLabelNeighbPix = m_label.at<int>(y, x-1);
			Value E00,E01,E10,E11 = 0;

			if (oldLabelPix != oldLabelNeighbPix)
                E00 = m_horizWeights.at<Value>(y, x-1);
			else
                E00 = 0;
			if (oldLabelNeighbPix != label)
                E01 = m_horizWeights.at<Value>(y, x-1);
			else
                E01 = 0;
			if (label != oldLabelPix)
                E10 = m_horizWeights.at<Value>(y, x-1);
			else
                E10 = 0;
			
			e->add_term2(variables[(x-startX)-1+(y-startY)*blockWidth],variables[(x-startX)+(y-startY)*blockWidth],E00,E01,E10,E11);
		}

	// Next set up vertical links
	for (int y = startY+1; y <= endY; y++)
		for (int x = startX; x <= endX; x++){
			int oldLabelPix       = m_label.at<int>(y, x);
			int oldLabelNeighbPix = m_label.at<int>(y-1, x);
			Value E00,E01,E10,E11 = 0;

			if (oldLabelPix != oldLabelNeighbPix)
                E00 = m_vertWeights.at<Value>(y-1, x);
			else
                E00 = 0;
			if (oldLabelNeighbPix != label)
                E01 = m_vertWeights.at<Value>(y-1, x);
			else
                E01 = 0;
			if (label != oldLabelPix)
                E10 = m_vertWeights.at<Value>(y-1, x);
			else
                E10 = 0;
			
			e->add_term2(variables[(x-startX)+(y-startY-1)*blockWidth],variables[(x-startX)+(y-startY)*blockWidth],E00,E01,E10,E11);
		}

	// Next set up diagonal links 
	float SQRT_2= 1/sqrt(2.0);
	for (int y = startY+1; y <= endY; y++)
		for (int x = startX+1; x <= endX; x++){
			int oldLabelPix       = m_label.at<int>(y, x);
			int oldLabelNeighbPix = m_label.at<int>(y-1, x-1);
			Value E00,E01,E10,E11=0;
			
			if (oldLabelPix != oldLabelNeighbPix) 
				E00 = SQRT_2*m_diag1Weights.at<Value>(y-1, x-1);
			else E00 = 0;
			if (oldLabelNeighbPix != label) 
				E01 = SQRT_2*m_diag1Weights.at<Value>(y-1, x-1);
			else E01 = 0;
			if (label != oldLabelPix)
				E10 = SQRT_2*m_diag1Weights.at<Value>(y-1, x-1);
            else E10 = 0;
			
			e->add_term2(variables[(x-startX)-1+(y-startY-1)*blockWidth],variables[(x-startX)+(y-startY)*blockWidth],E00,E01,E10,E11);
		}
	
	// More diagonal links
	for (int y = startY+1; y <= endY; y++)
		for (int x = startX; x <=endX-1; x++){
			int oldLabelPix       = m_label.at<int>(y, x);
			int oldLabelNeighbPix = m_label.at<int>(y-1, x+1);
			Value E00,E01,E10,E11=0;

			if (oldLabelPix != oldLabelNeighbPix) 
				E00 = SQRT_2*m_diag2Weights.at<Value>(y-1, x+1);
			else
                E00 = 0;
			if (oldLabelNeighbPix != label) 
				E01 = SQRT_2*m_diag2Weights.at<Value>(y-1, x+1);
			else
                E01 = 0;
			if (label != oldLabelPix)
				E10 = SQRT_2*m_diag2Weights.at<Value>(y-1, x+1);
			else
                E10 = 0;

			e->add_term2(variables[(x-startX+1)+(y-startY-1)*blockWidth],variables[(x-startX)+(y-startY)*blockWidth],E00,E01,E10,E11);
		}

	e->minimize();

	for (int y = startY; y <= endY; y++)
		for (int x = startX; x <= endX; x++){
			if (e->get_var(variables[(x-startX)+(y-startY)*blockWidth]) != 0)
			{
				if (m_label.at<int>(y, x) != label){
                    m_label.at<int>(y, x) = label; 
					changeMaskNew[x+y*m_width] = 1;
					changeMask[x+y*m_width] = 1;
				}
			}
		}

	delete e;
}

void SuperPixels::initializeLabeling()
{
    assert(m_init && m_label.type() == CV_32S);
    
    int seedX, seedY, startX, startY, endX, endY;
	for (size_t i = 0; i < m_seeds.size(); i++)
    {
		seedY = m_seeds[i].y;
		seedX = m_seeds[i].x;

		startX = seedX-m_patch_size/2-1;
		endX   = seedX+m_patch_size/2+1;

		startY = seedY-m_patch_size/2-1;
		endY   = seedY+m_patch_size/2+1;

		if (startX < 0)       startX  = 0;
		if (startY < 0)       startY  = 0;
		if (endX >= m_width)  endX    = m_width-1;
		if (endY >= m_height) endY    = m_height-1;

		for (int y = startY; y <= endY; y++)
			for (int x = startX; x <= endX; x++)
				m_label.at<int>(y, x) = i;
	}
}

float SuperPixels::computeImageVariance()
{
    assert(m_init);
    
    int total = 0;
	float v = (float) 0.0;

	for (int y = 1; y < m_height; y++)
		for (int x = 1; x < m_width; x++){
            int off = m_gray.at<uchar>(y, x);
            int offx_1 = m_gray.at<uchar>(y, x-1);
            int offy_1 = m_gray.at<uchar>(y-1, x);

            v += abs(off - offx_1) + abs(off - offy_1);
			total += 2;
        }
    
	return (v/total);
}

void SuperPixels::generate()
{
    assert(m_init);
    
	// Initialize and place seeds
	PlaceSeeds();
	initializeLabeling();
		
	Value oldEnergy, newEnergy;
    int num_pixels = m_width * m_height;
    int numSeeds = m_seeds.size();
	vector<int> changeMask(num_pixels, 1), changeMaskNew(num_pixels, 0);
    vector<int> order(numSeeds);
	for (int i = 0; i < numSeeds; i++)
		order[i] = i;

	//purturbSeeds(order,numSeeds);

    for (int j = 0; j < m_numIter; ++j){
		newEnergy = computeEnergy();

		if (j == 0){
			oldEnergy = newEnergy+1;
#ifdef DEBUG
			printf("Initial Energy: %d\n", newEnergy);
#endif
		}
#ifdef DEBUG
		else {
			printf("After iteration %d: ", j-1);
			printf("Energy: %d,  %f sec\n", newEnergy, ((float)clock())/CLOCKS_PER_SEC);
		}
#endif
		if (newEnergy == oldEnergy) break;

		oldEnergy = newEnergy;
		
		for (int i = 0; i < numSeeds; i++){
			expandOnLabel(order[i], changeMask, changeMaskNew);
		}
		for (int i = 0; i < num_pixels; i++){
			changeMask[i] = changeMaskNew[i];
			changeMaskNew[i] = 0;
		}

		//purturbSeeds(order,numSeeds);
	}

#ifdef DEBUG	
	printf("Final Energy:   %d\n",newEnergy);
#endif
}

void SuperPixels::purturbSeeds()
{
    assert(m_init);
    
    int numSeeds = m_seeds.size();
    srand(time(NULL));
	for (int i = 0; i < 3*numSeeds; i++)
	{
		int first  = (rand()*rand())%numSeeds;
		int second = (rand()*rand())%numSeeds;
		Point temp = m_seeds[first];
		m_seeds[first] = m_seeds[second];
		m_seeds[second] = temp;
	}
}

void SuperPixels::computeWeights()
{
    computeWeights(m_horizWeights, -1, 0);
    computeWeights(m_vertWeights, 0, -1);
    computeWeights(m_diag1Weights, -1, -1);
    computeWeights(m_diag2Weights, 1, -1);
}

void SuperPixels::computeWeights(Mat &weights, int incrX, int incrY)
{
    assert(m_init);

    weights.create(m_height, m_width, CV_VALUE);

    int startX = 0, startY = 0;
	float sigma = 2.0f;
    float var2 = sq(m_variance);
    
	if (incrX != 0) startX  = abs(incrX);
	if (incrY != 0) startY = abs(incrY);

	Value smallPenalty;
	if (m_TYPE == 1) smallPenalty = (MULTIPLIER_VAR*m_variance)/8+1;
	else smallPenalty = 1;

	for (int y = startY; y < m_height; y++)
		for (int x = startX; x < m_width; x++){
            int c1 = m_gray.at<uchar>(y, x);
            int c2 = m_gray.at<uchar>(y+incrY, x+incrX);
			int difference = sq((c1 - c2));
            weights.at<Value>(y+incrY, x+incrX) = (Value) (m_lambda*exp(-difference/(sigma*var2))+smallPenalty);
		}
}

void SuperPixels::loadHEdge(const char *name)
{
    loadEdge(m_horizWeights, name);
}

void SuperPixels::loadVEdge(const char *name)
{
    loadEdge(m_vertWeights, name);
}

void SuperPixels::loadD1Edge(const char *name)
{
    loadEdge(m_diag1Weights, name);
}

void SuperPixels::loadD2Edge(const char *name)
{
    loadEdge(m_diag2Weights, name);
}

void SuperPixels::loadEdge(Mat &weights, const char *name)
{
    assert(m_init && weights.type() == CV_VALUE);
    
	Mat edges = imread(name, CV_LOAD_IMAGE_GRAYSCALE);
    assert(m_height == edges.rows && m_width == edges.cols);
    
	for (int y = 0; y < m_height; y++)
		for (int x = 0; x < m_width; x++){
            int c = edges.at<uchar>(y, x);
			weights.at<Value>(y, x) = (Value) m_lambda * c;
		}
}

int SuperPixels::countLabels()
{
    assert(m_init);

    int numSeeds = m_seeds.size();
    vector<int> counts(numSeeds, 0);
    for (int h = 0; h < m_height; ++h)
        for (int w = 0; w < m_width; ++w)
        {
            int label = m_label.at<int>(h, w);
 			counts[label]++;
        }

    // re-arrange labels, to make label id continous
    for (size_t i = 0; i < counts.size(); ++i)
    {
        if (counts[i] > 0) continue;

        for (int size = counts.size() - 1; size > 0; --size)
        {
            if (counts[size] == 0) counts.pop_back();
            else
            {
                for (int h = 0; h < m_height; ++h)
                    for (int w = 0; w < m_width; ++w)
                    {
                        if (size == m_label.at<int>(h, w))
                            m_label.at<int>(h, w) = i;
                    }
                counts.pop_back();
                break;
            }
        }
    }

    return counts.size();
}

int SuperPixels::saveSegmentationColor(const char *name)
{
    assert(m_init);
    
    Mat out(m_height, m_width, CV_8UC3);

    int numSeeds = m_seeds.size();
	vector<Vec3b> colorLookup(numSeeds);
	vector<int> counts(numSeeds,0);

    // Vec3b: BGR
    srand(time(NULL));
	for (int i = 0; i < numSeeds; i++){
		colorLookup[i][2] = rand()%NUM_COLORS;
		colorLookup[i][1] = rand()%NUM_COLORS;
		colorLookup[i][0] = rand()%NUM_COLORS;
	}

	for (int y = 0; y < m_height; y++)
		for (int  x = 0; x < m_width; x++ ){
            int label = m_label.at<int>(y, x);
			out.at<Vec3b>(y, x) = colorLookup[label];
			counts[label]++;
		}

    imwrite(name, out);

	int num_superpixels = 0;
	for (int i = 0; i < numSeeds; i++)
		if (counts[i] > 0 ) num_superpixels++;

	return (num_superpixels);
}

int SuperPixels::saveSegmentationEdges(const char *name)
{
    assert(m_init);

    Mat dst;
    if (m_im.channels() == 3)
        m_im.copyTo(dst);
    else{
        Mat arr[3] = {m_im, m_im, m_im};
        merge(arr, 3, dst);
    }

    int num_superpixels = countLabels();

    Vec3b green(0, 255, 0);
    for (int h = 1; h < m_height; ++h)
        for (int w = 1; w < m_width; ++w)
        {
            int label = m_label.at<int>(h, w);
            if (label != m_label.at<int>(h, w-1)) {
                dst.at<Vec3b>(h, w) = green;
                dst.at<Vec3b>(h, w-1) = green;
            }
            if (label != m_label.at<int>(h-1, w)) {
                dst.at<Vec3b>(h, w) = green;
                dst.at<Vec3b>(h-1, w) = green;
            }
        }

    imwrite(name, dst);

    return num_superpixels;
}
