#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include <cmath>  
class bilateral_filter{
public:
	bilateral_filter(){}
	~bilateral_filter(){}
	static IplImage* bilateral_Filters(char *filename);
};