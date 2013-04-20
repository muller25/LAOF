#include "GaussianBlur.h"
/******************************************************************************/
/*** C Headers                                                              ***/
/******************************************************************************/
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include "highgui.h"
#include <omp.h>
using namespace std;

void GaussianBlur::gaussian_blur(IplImage * dst_image, const IplImage * src_image, float sigma)
{
    clock_t   start,   finish;   
    double     duration;
    start=clock();
	int n = (int) (sigma + 0.5) * 2 + 1;    // width and height of the gaussian kernel
    //int n=5;
	int n_h = n / 2;                        // n divided by 2 for optimizations
	float * matrix = (float *) malloc(sizeof(float) * n * n);   // gaussian kernel
	float sum = 0;
	int i, j, x, y;

	//#pragma omp parallel for private(y) reduction(+:sum)
	for (x = -n_h; x <= n_h; x++)
	{
		for (y = -n_h; y <= n_h; y++)
		{
			matrix[(x + n_h) * n + (y + n_h)] = (float) (/*1.0 / (2 * 3.141592 * sigma * sigma)**/  exp(-(x * x + y * y) / (2 * sigma * sigma)));
			sum += matrix[(x + n_h) * n + (y + n_h)];
		}
	}
	//#pragma omp parallel for private(y)
	for (x = -n_h; x <= n_h; x++)
	{
		for (y = -n_h; y <= n_h; y++)
		{
			matrix[(x + n_h) * n + (y + n_h)] /= sum;
		}
	}
	
// #pragma omp parallel for private(j)//这里改并行
	for (i = 0; i < src_image->width; i++)
	{
		for (j = 0; j < src_image->height; j++)
		{
			int ii, jj;
			float r = 0.0, g = 0.0, b = 0.0;
			int base_index;
			//#pragma omp parallel for private(jj,base_index) lastprivate(r,g,b)
			for (ii = 0; ii < n; ii++)
			{
				for (jj = 0; jj < n; jj++)
				{
					int iii = i - n_h + ii, jjj = j - n_h + jj;
					int kernel_index = ii * n + jj;

					if (iii < 0)
						iii = i;
					else if (iii >= src_image->width)
						iii = i;

					if (jjj < 0)
						jjj = j;
					else if (jjj >= src_image->height)
						jjj = j;

					base_index = (jjj) * src_image->widthStep + (iii)*3;//一定用widthStep；
					r += (uchar)src_image->imageData[base_index + 0] * matrix[kernel_index];//char类型一定要转为uchar
					g += (uchar)src_image->imageData[base_index + 1] * matrix[kernel_index];
					b += (uchar)src_image->imageData[base_index + 2] * matrix[kernel_index];
				}
			}
			
			base_index = j * dst_image->widthStep + i*3;
			dst_image->imageData[base_index] = r;
			dst_image->imageData[base_index + 1] = g;
			dst_image->imageData[base_index + 2] = b;
		}
#ifdef SHOWMODE
        cvShowImage("show",dst_image);
        cvWaitKey(1);
#endif

	}

	free(matrix);

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout<<"time of running gaussian blur(): "<< duration <<endl;
}
