#ifndef _Edge_Optimization_H
#define _Edge_Optimization_H

#include<cv.h>

class Edge_Optimization{

public:
	Edge_Optimization(void);
	~Edge_Optimization(void);
	static IplImage * Optimization_edge(IplImage * src,IplImage * ref,int para);

};

#endif
