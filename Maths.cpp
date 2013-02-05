#include "Maths.h"

/*
  [description]
  use central difference to approximate
  first order of derivative x of matrix src
  [params]
  src - input matrix (in)
  [return]
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
  [description]
  use central difference to approximate
  first order of derivative y of matrix src
  [params]
  src - input matrix (in)
  [return]
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
  [description]
  use central difference to approximate
  second order of derivative x of matrix src
  [params]
  src - input matrix (in)
  [return]
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
  [description]
  use central difference to approximate
  second order of derivative y of matrix src
  [params]
  src - input matrix (in)
  [return]
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
  [description]
  use taylor series to approximate
  derivative x-y of matrix src
  [params]
  src - input matrix (in)
  [return]
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

/*
  [description]
  first apply backward difference, then apply forward difference
  to approximate weighted laplacian
  [params]
  flow - flow matrix, it may contain multi-channel (in)
  weight - weight matrix, it only contains one-channel (in)
  [return]
  weighted laplacian matrix
 */
Mat Maths::weighted_laplacian(const Mat &flow, const Mat &weight)
{
    int rows = flow.rows;
    int cols = flow.cols;
    int channels = flow.channels();
    int step = flow.step / get_step(flow.depth());
   
//    assert(flow.size() == weight.size() && flow.depth() == CV_64F &&
//           weight.channels() == 1 && weight.depth() == CV_64F);
    
    Mat lap = Mat::zeros(rows, cols, CV_32SC(channels));
    
    int r, c, k, offset;
    int woffset, wstep = weight.step / get_step(weight.depth());
    int *fptr, *wptr, *lptr;

    fptr = (int *)flow.data;
    wptr = (int *)weight.data;
    lptr = (int *)lap.data;
    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            woffset = r * wstep + c;
            for (k = 0; k < channels; k++)
            {
                offset = r * step + c * channels + k;
                if (c < cols-1)
                    lptr[offset] += wptr[woffset+1] * (fptr[offset+channels] - fptr[offset]);
                if (c > 0)
                    lptr[offset] -= wptr[woffset] * (fptr[offset] - fptr[offset-channels]);
                if (r < rows-1)
                    lptr[offset] += wptr[woffset+wstep] * (fptr[offset+step] - fptr[offset]);
                if (r > 0)
                    lptr[offset] -= wptr[woffset] * (fptr[offset] - fptr[offset-step]);
            }
        }
    }

    return lap;
}

/*
 */
Mat Maths::weighted_laplacian3D(const Mat &pflow, const Mat &flow, const Mat &nflow,
                                const Mat &weight, const Mat &nweight)
{
    // assert(match(pflow, flow, true) && match(flow, nflow, true) && flow.channels() == 2 &&
    //        match(weight, nweight, true) && weight.channels() == 1 && 
    //        flow.size() == weight.size() &&
    //        weight.depth() == CV_64F && flow.depth() == CV_64F);
    
    int rows = flow.rows;
    int cols = flow.cols;
    int channels = flow.channels();
    int step = flow.step / get_step(flow.depth());
    int wstep = weight.step / get_step(weight.depth());
    int r, c, k, offset, nr, nc, noffset, woffset;
    double *pfptr, *fptr, *nfptr, *wptr, *nwptr, *l3ptr;

    Mat lap3d = weighted_laplacian(flow, weight);
    l3ptr = (double *)lap3d.data;
    pfptr = (double *)pflow.data;
    fptr = (double *)flow.data;
    nfptr = (double *)nflow.data;
    wptr = (double *)weight.data;
    nwptr = (double *)nweight.data;

    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            offset = r * step + c * channels;
            nr = r + pfptr[offset];
            nc = c + pfptr[offset+1];
            if (nr < 0 || nr >= rows || nc >= cols || nc < 0) continue;

            woffset = r * wstep + c;
            for (k = 0; k < channels; k++)
            {
                noffset = nr * step + nc * channels + k;
                l3ptr[noffset] -= wptr[woffset] * (fptr[noffset] - pfptr[offset+k]);
            }
        }
    }

    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            offset = r * step + c * channels;
            nr = fptr[offset] + r;
            nc = fptr[offset+1] + c;
            if (nr < 0 || nr >= rows || nc >= cols || nc < 0) continue;

            noffset = nr * step + nc * channels;
            woffset = nr * wstep + nc;
            for (k = 0; k < channels; k++)
                l3ptr[offset+k] += nwptr[woffset] * (nfptr[noffset+k] - fptr[offset+k]);
        }
    }
    
    return lap3d;
}

/*
  [desciption]
  Successive Over-Relaxation solver for linear equation system Ax = b
  [params]
  A - cofficient matrix (in)
  b - column vector (in)
  x - variables column vector, it may contain init values (in, out)
  nIters - number of iterations (in)
  w - relation rate, should be in (0, 2) for convergence (in)
  [return]
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

/*
  [description]
  match two matrix size, channel and depth
  [params]
  m1, m2 - matrix to match (in)
  depth - whether match depth or not (in)
  [return]
  true if two matrix match in all aspect, otherwise returns false
 */
bool match(const Mat &m1, const Mat &m2, bool depth)
{
    bool size_match = (m1.rows == m2.rows && m1.cols == m2.cols &&
                       m1.channels() == m2.channels());

    if (depth)
        return (size_match && m1.depth() == m2.depth());

    return size_match;
}

int get_step(const int depth)
{
    int type = -1;
    switch (depth)
    {
    case CV_64F:
        type = sizeof(double);
        break;
    case CV_32F:
        type = sizeof(float);
        break;
    case CV_32S:
        type = sizeof(int);
        break;
    case CV_16S:
    case CV_16U:
        type = sizeof(short);
        break;
    case CV_8S:
    case CV_8U:
        type = sizeof(char);
        break;
    }

    assert(type >= 0);
    return type;
}
