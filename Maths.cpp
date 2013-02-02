#include "Maths.h"

/*
  description:
  use central difference to approximate
  first order of derivative x of matrix src
  params:
  src - input matrix (in)
  return:
  matrix of dx
*/
Mat Maths::dx(const Mat &src)
{
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    double kernel_arr[] = {1./12, -8./12, 0, 8./12, -1./12};
    Mat kernel(1, 5, CV_64FC1, kernel_arr), dst;

//    cout << kernel << endl;
    filter2D(src, dst, CV_64FC(channels), kernel);
    return dst;
}

/*
  description:
  use central difference to approximate
  first order of derivative y of matrix src
  params:
  src - input matrix (in)
  return:
  matrix of dy
*/
Mat Maths::dy(const Mat &src)
{
    int r, c, k, offset;
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    double kernel_arr[5][1] = {{1./12},{-8./12},{0},{8./12},{-1./12}};
    Mat dst, kernel(5, 1, CV_64FC1, kernel_arr);
    
//    cout << kernel << endl;
    filter2D(src, dst, CV_64FC(channels), kernel);

    return dst;
}

/*
  description:
  use central difference to approximate
  second order of derivative x of matrix src
  params:
  src - input matrix (in)
  return:
  matrix of dxx
*/
Mat Maths::dxx(const Mat &src)
{
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    double kernel_arr[] = {-1./12,16./12,-30./12,16./12,-1./12};
    Mat dst, kernel(1, 5, CV_64FC1, kernel_arr);

    cout << kernel << endl;
    filter2D(src, dst, CV_64FC(channels), kernel);

    return dst;
}

/*
  description:
  use central difference to approximate
  second order of derivative y of matrix src
  params:
  src - input matrix (in)
  return:
  matrix of dyy
*/
Mat Maths::dyy(const Mat &src)
{
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    double kernel_arr[5][1] = {{-1./12},{16./12},{-30./12},{16./12},{-1./12}};
    Mat dst, kernel(5, 1, CV_64FC1, kernel_arr);

    cout << kernel << endl;
    filter2D(src, dst, CV_64FC(channels), kernel);
    
    return dst;
}

/*
  description:
  use taylor series to approximate
  derivative x-y of matrix src
  params:
  src - input matrix (in)
  return:
  matrix of dxy
*/
Mat Maths::dxy(const Mat &src)
{
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    double kernel_arr[3][3] = {
        {.25, 0, -.25},
        {0,   0,    0},
        {-.25,0,  .25}
    };
    Mat dst, kernel(3, 3, CV_64FC1, kernel_arr);

    cout << kernel << endl;
    filter2D(src, dst, CV_64FC(channels), kernel);
    
    return dst;
}

Mat Maths::laplace3D(const Mat &im0, const Mat &im1, const Mat &im2)
{
    int rows = im1.rows;
    int cols = im1.cols;
    int channels = im1.channels();

    assert(im0.rows == rows && rows == im2.rows &&
           im0.cols == cols && cols == im2.cols &&
           im0.channels() == channels && channels == im2.channels());

    double kernel_arr[3][3] = {
        {0,  1.,  0},
        {1.,-6., 1.},
        {0,  1.,  0}
    };
    Mat ctmp, tmp, dst, kernel(3, 3, CV_64FC1, kernel_arr);
    
    filter2D(im1, dst, CV_64FC(channels), kernel);
    cout << "***** dst *****" << endl << dst << endl;
    
    tmp = im0 + im2;
    tmp.convertTo(ctmp, CV_64FC(channels));

    cout << "***** ctmp *****" << endl << ctmp << endl;
    dst += ctmp;
    
    return dst;
}

/*
  desciption:
  Successive Over-Relaxation solver for linear equation system Ax = b
  params:
  A - cofficient matrix (in)
  b - column vector (in)
  x - variables column vector, it may contain init values (in, out)
  nIters - number of iterations (in)
  w - relation rate, should be in (0, 2) for convergence (in)
  return:
  0 for success, others for failure
*/
int Maths::sor_solver(const Mat &A, const Mat &b, Mat &x,
                      int nIters, double w)
{
    // sanity check
    if (nIters <= 0 || w >= 2 || w <= 0) return PARAMS_CHECK_FAIL;

    // match dimension
    if (A.rows != b.rows || A.cols != x.rows || x.cols > 1 || b.cols > 1)
        return PARAMS_CHECK_FAIL;
    
    int iter, r, c;
    double tmp;
    const double *a_ptr;
    for (iter = 0; iter < nIters; iter++)
    {
        for (r = 0; r < x.rows; r++)
        {
            a_ptr = A.ptr<double>(r);
            for (tmp = b.at<double>(r, 0), c = 0; c < A.cols; c++)
                tmp -= a_ptr[c] * x.at<double>(c, 0);
            
            x.at<double>(r, 0) = x.at<double>(r, 0) + w * tmp / A.at<double>(r, r);
        }
    }

    return SUCCESS;
}
