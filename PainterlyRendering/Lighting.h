#ifndef _Lighting_H
#define _Lighting_H

#include<cv.h>
#include<math.h>

class Lighting
{
public:
	Lighting(void);
	~Lighting(void);
	static IplImage* phong(IplImage* src, IplImage* nmap);
};

class vector3i
{
public:
	int x,y,z;
	float operator *(vector3i& a)
	{
		return (x*a.x+y*a.y+z+a.z)/((sqrt((float)x*x+y*y+z*z))*(sqrt((float)a.x*a.x+a.y*a.y+a.z*a.z)));

	}
	vector3i operator -(vector3i& a)
	{
		vector3i res(x-a.x,y-a.y,z-a.z);
		return res;
	}
	vector3i& operator *(float a)
	{
		x*=a;
		y*=a;
		z*=a;
		return *this;
	}
	vector3i(int a,int b,int c):x(a),y(b),z(c)
	{
	}
	vector3i():x(0),y(0),z(0)
	{
	}
};

#endif
