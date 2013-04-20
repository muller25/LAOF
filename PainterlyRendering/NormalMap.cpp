#include "NormalMap.h"
#include <highgui.h>

NormalMap::NormalMap(void){}
NormalMap::~NormalMap(void){}

IplImage * NormalMap::computeNormalMap(IplImage * hmap)//ehgraymap
{
//    int q = hmap->nChannels;
    CvSize csize;
    csize.width = hmap->width;
    csize.height = hmap->height;
 
    IplImage *nP = cvCreateImage(csize, 8, 3);  //创建一个带有RBG三个通道的深度为8的法线图
 
    CvScalar s;
    CvScalar r;
    int Hg,Ha,Hr,Hb,Hl;
    for(int i = 0;i<csize.height;i++)
    {
        for(int j = 0 ;j<csize.width;j++)
        {
            s=cvGet2D(hmap,i,j);
            Hg = s.val[0];
 
            if(i-1<0)   Ha = 0;
            else
            {
                s = cvGet2D(hmap,i-1,j);
                Ha = s.val[0];
            }
            if(i+1>=csize.height)   Hb = 0;
            else
            {
                s = cvGet2D(hmap,i+1,j);
                Hb = s.val[0];
            }
            if(j+1>=csize.width)   Hr = 0;
            else
            {
                s = cvGet2D(hmap,i,j+1);
                Hr = s.val[0];
            }
            if(j-1<0)   Hl = 0;
            else
            {
                s = cvGet2D(hmap,i,j-1);
                Hl = s.val[0];
            }
 
            r = cvGet2D(nP, i, j);
            r.val[2] = (Hb-Ha)*0.1+128;     // 为了保存法线图为位图，这里将法线转换为[0,255]之间的整数
            r.val[1] = (Hl-Hr)*0.1+128;//(Hb-Ha)*0.2中系数是控制厚度
            r.val[0] = 255;
            cvSet2D(nP,i,j,r);
        }
    }
    return nP;
}
