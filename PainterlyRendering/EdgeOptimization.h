#ifndef _Edge_Optimization_H
#define _Edge_Optimization_H

#include<cv.h>

class EdgeOptimization{

public:
	EdgeOptimization(void);
	~EdgeOptimization(void);
	static IplImage * OptimizationEdge(IplImage * src,IplImage * ref,int para);

};

#endif
