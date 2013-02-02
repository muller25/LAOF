#include "Maths.h"

/*
  description:
  use forward difference to approximate
  first order of derivative x of matrix src
  params:
  src - input matrix (in)
  return:
  matrix of dx
*/
Mat Maths::dx(const Mat &src)
{
    int r, c, k, offset;
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    Mat tmp, dst(rows, cols, CV_64FC(channels));
    int step = dst.step / sizeof(double);
    double *ptr, *dptr;

    src.convertTo(tmp, CV_64FC(channels));
    ptr = (double *)tmp.data;
    dptr = (double *)dst.data;
    for (r = 0; r < rows; r++)
        for (c = 0; c < cols; c++)
            for (k = 0; k < channels; k++)
            {
                offset = r * step + c * channels + k;
                // use reflect to extrapolate pixels out of boundary
                if (c == cols - 1)
                    dptr[offset] = ptr[offset-channels] - ptr[offset];
                else
                    dptr[offset] = ptr[offset+channels] - ptr[offset];
            }

    return dst;
}

/*
  description:
  use forward difference to approximate
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
    Mat tmp, dst(rows, cols, CV_64FC(channels));
    int step = dst.step / sizeof(double);
    double *ptr, *dptr;
    
    src.convertTo(tmp, CV_64FC(channels));
    ptr = (double *)tmp.data;
    dptr = (double *)dst.data;
    for (r = 0; r < rows; r++)
        for (c = 0; c < cols; c++)
            for (k = 0; k < channels; k++)
            {
                offset = r * step + c * channels + k;
                // use reflect to extrapolate pixels out of boundary
                if (r == rows - 1)
                    dptr[offset] = ptr[offset-step] - ptr[offset];
                else
                    dptr[offset] = ptr[offset+step] - ptr[offset];
            }

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
    int r, c, k, offset;
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    Mat tmp, dst(rows, cols, CV_64FC(channels));
    int step = dst.step / sizeof(double);
    double *ptr, *dptr;

    src.convertTo(tmp, CV_64FC(channels));
    ptr = (double *)tmp.data;
    dptr = (double *)dst.data;
    for (r = 0; r < rows; r++)
        for (c = 0; c < cols; c++)
            for (k = 0; k < channels; k++)
            {
                offset = r * step + c * channels + k;
                
                // use reflect to extrapolate pixels out of boundary
                if (c == 0)
                    dptr[offset] = 2 * ptr[offset+channels] - 2 * ptr[offset];
                else if (c == cols - 1)
                    dptr[offset] = 2 * ptr[offset-channels] - 2 * ptr[offset];
                else
                    dptr[offset] = ptr[offset+channels] + ptr[offset-channels] - 2 * ptr[offset];
            }

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
    int r, c, k, offset;
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    Mat tmp, dst(rows, cols, CV_64FC(channels));
    int step = dst.step / sizeof(double);
    double *ptr, *dptr;

    src.convertTo(tmp, CV_64FC(channels));
    ptr = (double *)tmp.data;
    dptr = (double *)dst.data;
    for (r = 0; r < rows; r++)
        for (c = 0; c < cols; c++)
            for (k = 0; k < channels; k++)
            {
                offset = r * step + c * channels + k;
                
                // use reflect to extrapolate pixels out of boundary
                if (r == 0)
                    dptr[offset] = 2 * ptr[offset+step] - 2 * ptr[offset];
                else if (r == rows - 1)
                    dptr[offset] = 2 * ptr[offset-step] - 2 * ptr[offset];
                else
                    dptr[offset] = ptr[offset+step] + ptr[offset-step] - 2 * ptr[offset];
            }

    return dst;
}

/*
  description:
  use central difference to approximate
  derivative x-y of matrix src
  params:
  src - input matrix (in)
  return:
  matrix of dxy
*/
Mat Maths::dxy(const Mat &src)
{
    int r, c, k, offset;
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    Mat tmp, dst(rows, cols, CV_64FC(channels));
    int step = dst.step / sizeof(double);
    double *ptr, *dptr;

    src.convertTo(tmp, CV_64FC(channels));
    ptr = (double *)tmp.data;
    dptr = (double *)dst.data;
    for (r = 0; r < rows; r++)
        for (c = 0; c < cols; c++)
            for (k = 0; k < channels; k++)
            {
                offset = r * step + c * channels + k;
                
                // use reflect to extrapolate pixels out of boundary
                if (r == 0)
                    dptr[offset] = 0;
                else if (r == rows - 1)
                    dptr[offset] = 0;
                else
                {
                    if (c == 0 || c == cols-1)
                        dptr[offset] = 0;
                    else
                    {
                        dptr[offset] = ptr[offset+step+channels] + ptr[offset-step-channels]
                            - ptr[offset+step-channels] - ptr[offset-step+channels];
                        dptr[offset] /= 4;
                    }
                }
            }

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
