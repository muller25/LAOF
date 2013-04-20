#ifndef RENDERING_H
#define RENDERING_H

#include <iostream>
#include <fstream>
#include  <io.h>  
#include <direct.h>  
#include <string>
#include  "time.h"  
#include <vector>  
#include <iomanip>  
#include <ctime> 
#include <stdio.h>
#include <stdlib.h>
#include <pshpack2.h>
#include "GaussianBlur.h"
#include "PainterlyService.h"
#include <queue>
#include <windows.h>
#include "RenderService.h"
#include "Lighting_Effect.h"
#include "Edge_Optimization.h"
#include<QImage>
#include"ImageOperate.h"
#include"MainWindow.h"
using namespace std;


class Rendering{

public:
	Rendering();
	~Rendering();
     static IplImage* Processing(const QImage * src,QImage* res,IplImage* edgeImage);
	 static IplImage* operateLight( QImage* res );
	 static IplImage* opereateEdge( QImage* res );
public:

	//static int max_len;
	static int max_stroke_length;
	static int min_stroke_length;
	static int threshold;
	static float grid_size;
	static int R[6];

};
#endif