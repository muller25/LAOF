#ifndef _NormalMap_H
#define _NormalMap_H

#include<cv.h>

class NormalMap
{
public:
	NormalMap(void);
	~NormalMap(void);
	static IplImage * computeNormalMap(IplImage * hMap);
};

#endif
