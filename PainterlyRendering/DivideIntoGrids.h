#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include <cmath>
#include<queue>

#define THRESHOLD 500 
#define R_MAX 16
#define R_MIN 4
#define �� 0.01

struct node{
	CvPoint top;//���Ͻ�
	CvPoint bottom;//���½�
	int areas;
	friend bool operator < (node  a, node b){
		return a.areas < b.areas; //�ṹ���У�xС�����ȼ���
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

