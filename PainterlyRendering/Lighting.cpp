#include "Lighting.h"
#include<highgui.h>

Lighting::Lighting(void){}
Lighting::~Lighting(void){}

IplImage* Lighting::phong(IplImage* srcImage, IplImage* nmap)
{
	int w=srcImage->width;
	int h=srcImage->height;
	IplImage* resImage=cvCreateImage(cvSize(w,h),8,3);
	float ka,kd,ks;//系数
	ka=0.5;
	ks=0.7;
	kd=0.2;
	float cosi,cosj;//cosi反射角。cosj镜面反射角。
	vector3i sun(0,0,100);//光源位置
	vector3i point;//当前点
	vector3i n,l,r;//法线，入射，反射
	vector3i v(0,0,1);//视点
	CvScalar src,nor,res;//原图，法线，结果
	for(int i=0;i<h;++i)
	{
		for(int j=0;j<w;++j)
		{
			nor=cvGet2D(nmap,i,j);
			src=cvGet2D(srcImage,i,j);

			point=vector3i(i,j,0);
			//l=sun-point;
			l=vector3i(100,100,1);
			n=vector3i(nor.val[0],nor.val[1],nor.val[2]);
			cosi=abs(l*n);
			r=n*(2*cosi)-l;
			cosj=abs(r*v);

			res.val[0]=ka*src.val[0]+kd*cosi*src.val[0]+ks*cosj*src.val[0];
			res.val[1]=ka*src.val[1]+kd*cosi*src.val[1]+ks*cosj*src.val[1];
			res.val[2]=ka*src.val[2]+kd*cosi*src.val[2]+ks*cosj*src.val[2];
			cvSet2D(resImage,i,j,res);
		}
	}
	return resImage;
}
