#ifndef _MakeBrushMap_H
#define _MakeBrushMap_H

#include<cv.h>

class MakeBrushMap{
public:
    MakeBrushMap();
    ~MakeBrushMap();
    static IplImage * NewBrushMap(int brushSize);


    
};

#endif
