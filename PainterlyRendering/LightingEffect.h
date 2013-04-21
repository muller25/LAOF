#ifndef _Light_Effect_H
#define _Light_Effect_H

#include<cv.h>

class LightingEffect{

public:
	LightingEffect(void);
	~LightingEffect(void);
	static IplImage * AddLight(IplImage * src,IplImage* height_map,IplImage * dst,float para);

};

#endif
