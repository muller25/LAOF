#include<cv.h>
#include<highgui.h>

class ColorAdjust{

public:
	ColorAdjust();
	~ColorAdjust();
	static IplImage* hsv_ad(IplImage *srcimg, int h, int s,int v,float nPercent);
};