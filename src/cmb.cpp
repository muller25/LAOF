#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <stdlib.h>
#include <stdio.h>

using namespace cv;

/// 全局变量
Mat src, dst;
int top, bottom, left, right;
int borderType;
Scalar value;
RNG rng(12345);

/** @函数 main  */
int main( int argc, char** argv )
{

  int c;

  /// 装载图像
  src = imread("car1.jpg");

  if( !src.data )
  { return -1;
    printf(" No data entered, please enter the path to an image file \n");
  }

  /// 初始化输入参数
  top = (int) (0.05*src.rows); bottom = (int) (0.05*src.rows);
  left = (int) (0.05*src.cols); right = (int) (0.05*src.cols);
  dst = src;

  namedWindow( "src", CV_WINDOW_AUTOSIZE );
  imshow( "src", dst );

  value = Scalar( rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255) );
  copyMakeBorder( src, dst, top, bottom, left, right, BORDER_CONSTANT, value );
  namedWindow( "constant", CV_WINDOW_AUTOSIZE );
  imshow( "constant", dst );

  copyMakeBorder( src, dst, top, bottom, left, right, BORDER_REPLICATE, value );
  namedWindow( "replicate", CV_WINDOW_AUTOSIZE );
  imshow( "replicate", dst );

  copyMakeBorder( src, dst, top, bottom, left, right, BORDER_REFLECT, value );
  namedWindow( "reflect", CV_WINDOW_AUTOSIZE );
  imshow( "reflect", dst );

  copyMakeBorder( src, dst, top, bottom, left, right, BORDER_REFLECT_101, value );
  namedWindow( "reflect 101", CV_WINDOW_AUTOSIZE );
  imshow( "reflect 101", dst );

  copyMakeBorder( src, dst, top, bottom, left, right, BORDER_WRAP, value );
  namedWindow( "wrap", CV_WINDOW_AUTOSIZE );
  imshow( "wrap", dst );

  copyMakeBorder( src, dst, top, bottom, left, right, BORDER_DEFAULT, value );
  namedWindow( "default", CV_WINDOW_AUTOSIZE );
  imshow( "default", dst );

  waitKey(0);

  return 0;
}
