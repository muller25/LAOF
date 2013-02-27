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
Mat dx(const Mat &src)
{
    int channels = src.channels();
    double kernel_arr[] = {1./12, -8./12, 0, 8./12, -1./12};
    Mat kernel(1, 5, CV_64F, kernel_arr), dst;

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
Mat dy(const Mat &src)
{
    int r, c, k, offset;
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    double kernel_arr[5][1] = {{1./12},{-8./12},{0},{8./12},{-1./12}};
    Mat kernel(5, 1, CV_64F, kernel_arr), dst;

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
Mat dxx(const Mat &src)
{
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    double kernel_arr[] = {-1./12,16./12,-30./12,16./12,-1./12};
    Mat kernel(1, 5, CV_64F, kernel_arr), dst;

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
Mat dyy(const Mat &src)
{
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    double kernel_arr[5][1] = {{-1./12},{16./12},{-30./12},{16./12},{-1./12}};
    Mat kernel(5, 1, CV_64F, kernel_arr), dst;

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
Mat dxy(const Mat &src)
{
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    double kernel_arr[3][3] = {
        {.25, 0, -.25},
        {0,   0,    0},
        {-.25,0,  .25}
    };
    Mat kernel(3, 3, CV_64F, kernel_arr), dst;

    filter2D(src, dst, CV_64FC(channels), kernel);
    return dst;
}

void collapse(Mat &src, Mat &dst)
{
    int rows = src.rows;
    int cols = src.cols;
    int channels = src.channels();
    double *ps = (double *)src.data, *pd = (double *)dst.data;
    int sstep = src.step / sizeof(double), dstep = dst.step / sizeof(double);
    int r, c, k, soffset, doffset;
    double tmp;
    
    if (channels == 1)
    {
        src.copyTo(dst);
        return;
    }
        
    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            soffset = r * sstep + c * channels;
            doffset = r * dstep + c;
            for (tmp = 0, k = 0; k < channels; k++)
                tmp += ps[soffset];
            
            pd[doffset] = tmp / channels;
        }
    }
}

/*
  [description]
  match two matrix in 2 dimension: row and col
  [params]
  m1, m2 - matrix to match (in)
  [return]
  true if two matrix match in size(row, col), otherwise returns false
 */
bool match2D(const Mat &m1, const Mat &m2)
{
    return (m1.rows == m2.rows && m1.cols == m2.cols);
}

/*
  [description]
  match two matrix in 3 dimension: row ,col and channel
  [params]
  m1, m2 - matrix to match (in)
  [return]
  true if two matrix match in row, col and channel, otherwise returns false
 */
bool match3D(const Mat &m1, const Mat &m2)
{
    return (match2D(m1, m2) && m1.channels() == m2.channels());
}

/*
  [description]
  match two matrix in all aspects: row, col, channel and depth
  [params]
  m1, m2 - matrix to match (in)
  depth - whether match depth or not (in)
  [return]
  true if two matrix match in all aspect, otherwise returns false
 */
bool matchAll(const Mat &m1, const Mat &m2)
{
    return (match3D(m1, m2) && m1.depth() == m2.depth());
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
void weighted_lap(const Mat &flow, const Mat &weight, Mat &lap)
{
    int r, c, k, offset, woffset, loffset;
    int rows = flow.rows;
    int cols = flow.cols;
    int channels = flow.channels();
    int step = flow.step / sizeof(double);
    double *fptr = (double *)flow.data;

    int wstep = weight.step / sizeof(double);
    double *wptr = (double *)weight.data;
    
    int lstep = lap.step / sizeof(double);
    double *lptr = (double *)lap.data;

    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            woffset = r * wstep + c;
            for (k = 0; k < channels; k++)
            {
                offset = r * step + c * channels + k;
                loffset = r * lstep + c * channels + k;
                if (c < cols-1)
                    lptr[loffset] += wptr[woffset+1] * (fptr[offset+channels] - fptr[offset]);
                if (c > 0)
                    lptr[loffset] -= wptr[woffset] * (fptr[offset] - fptr[offset-channels]);
                if (r < rows-1)
                    lptr[loffset] += wptr[woffset+wstep] * (fptr[offset+step] - fptr[offset]);
                if (r > 0)
                    lptr[loffset] -= wptr[woffset] * (fptr[offset] - fptr[offset-step]);
            }
        }
    }
}

/*
 */
/*void weighted_lap3D(const Mat &pflow, const Mat &flow, const Mat &nflow,
                    const Mat &weight, const Mat &nweight)
{
    // assert(match(pflow, flow, true) && match(flow, nflow, true) && flow.channels() == 2 &&
    //        match(weight, nweight, true) && weight.channels() == 1 && 
    //        flow.size() == weight.size() &&
    //        weight.depth() == CV_64F && flow.depth() == CV_64F);
    
    int rows = flow.rows;
    int cols = flow.cols;
    int channels = flow.channels();
    int step = 0;//flow.step / get_step(flow.depth());
    int wstep = 0;//weight.step / get_step(weight.depth());
    int r, c, k, offset, nr, nc, noffset, woffset;
    double *pfptr, *fptr, *nfptr, *wptr, *nwptr, *l3ptr;

    Mat lap3d = weighted_lap<double, double>(flow, weight);
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
*/
