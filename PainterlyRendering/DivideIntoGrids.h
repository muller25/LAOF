#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include <cmath>
#include<queue>

#define THRESHOLD 500 
#define R_MAX 16
#define R_MIN 4
#define δ 0.01

struct node{
	CvPoint top;//左上角
	CvPoint bottom;//右下角
	int areas;
	friend bool operator < (node  a, node b){
		return a.areas < b.areas; //结构体中，x小的优先级高
	}
};

class DivideIntoGrids{
public:
	DivideIntoGrids(void);
	~DivideIntoGrids(void);
	static double * energy;
	static int w;
	static int h;
	static void setData(IplImage* src);
	static double calculate_sum_energy(CvPoint start,CvPoint end);
	static std::priority_queue<node> make_grids(CvPoint start,CvPoint end);
};

