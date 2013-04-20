#ifndef _Light_Effect_H
#define _Light_Effect_H

#include<cv.h>

class Lighting_Effect{

public:
	Lighting_Effect(void);
	~Lighting_Effect(void);
	static IplImage * Add_Light(IplImage * src,IplImage* height_map,IplImage * dst,float para);

};

#endif
