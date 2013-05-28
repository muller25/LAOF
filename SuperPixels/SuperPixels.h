#ifndef _SuperPixels_H
#define _SuperPixels_H

#include <vector>
using std::vector;

#include <cv.h>
using cv::Mat;

typedef int Value;
typedef int TotalValue;
typedef int Var;

void PlaceSeeds(vector<int> &Seeds, int &numSeeds, const Mat &I, int PATCH_SIZE);

Value computeEnergy(const Mat &I, const vector<int> &labeling,
					const vector<Value> &horizWeights, const vector<Value> &vertWeights,
					const vector<Value> &diag1Weights, const vector<Value> &diag2Weights,
					const vector<int> &Seeds, int TYPE);

void getBounds(int &seedX, int &seedY, int &startX, int &startY, int &endX, int &endY,
               const Mat &I, const vector<int> &Seeds, int label, int PATCH_SIZE);

void expandOnLabel(int label, vector<int> &labeling,
                   vector<int> &changeMask, vector<int> &changeMaskNew,
                   const Mat &I, const vector<int> &Seeds, int numSeeds, 
				   const vector<Value> &horizWeights, const vector<Value> &vertWeights, 
				   const vector<Value> &diag1Weights, const vector<Value> &diag2Weights,
                   Value lambda, int PATCH_SIZE, int TYPE, float variance);

void initializeLabeling(vector<int> &labeling, const vector<int> &Seeds, const Mat &I,
                        int numSeeds, int PATCH_SIZE);

float computeImageVariance(const Mat &I);

void purturbSeeds(vector<int> &order,int numSeeds);

void computeWeights(vector<Value> &weights, const Mat &I, Value lambda, float variance,
					int incrX, int incrY, int TYPE);

void loadEdges(vector<Value> &weights, const Mat &I, Value lambda, const char *name);

int saveSegmentationColor(const Mat &I, const vector<int> &labeling, int numSeeds,
                          const char *name);

int saveSegmentationEdges(const Mat &I, const vector<int> &labeling, int numSeeds, const char *name);

#endif
