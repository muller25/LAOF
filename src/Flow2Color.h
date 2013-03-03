#ifndef _Flow2Color_H
#define _Flow2Color_H

#include <cstdio>

#include "Image.h"
#include "Maths.h"

void flow2color(UCImage &flowImg, UCImage &idxImg, const DImage &u, const DImage &v);
void computeColor(UCImage &im, const DImage &u, const DImage &v);
DImage makeColorWheel();

#endif
