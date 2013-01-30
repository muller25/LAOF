#include "OpticalFlow.h"

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
int OpticalFlow::SOR_solver(const Mat &A, const Mat &b, Mat &x,
                            int nIters, double w)
{
    // sanity check
    if (nIters <= 0 || w >= 2 || w <= 0) return -1;

    // match dimension
    if (A.rows != b.rows || A.cols != x.rows || x.cols > 1 || b.cols > 1)
        return -1;
    
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

    return 0;
}
