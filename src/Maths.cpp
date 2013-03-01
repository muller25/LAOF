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
    double kernel_arr[] = {1./12, -8./12, 0, 8./12, -1./12};
    Mat kernel(1, 5, CV_64F, kernel_arr), dst;

    filter2D(src, dst, src.type(), kernel, Point(-1, -1), 0, BORDER_REPLICATE);
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
    double kernel_arr[5][1] = {{1./12},{-8./12},{0},{8./12},{-1./12}};
    Mat kernel(5, 1, CV_64F, kernel_arr), dst;

    filter2D(src, dst, src.type(), kernel, Point(-1, -1), 0, BORDER_REPLICATE);
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
    double kernel_arr[] = {-1./12, 16./12, -30./12, 16./12, -1./12};
    Mat kernel(1, 5, CV_64F, kernel_arr), dst;

    filter2D(src, dst, src.type(), kernel, Point(-1, -1), 0, BORDER_REPLICATE);
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
    double kernel_arr[5][1] = {{-1./12},{16./12},{-30./12},{16./12},{-1./12}};
    Mat kernel(5, 1, CV_64F, kernel_arr), dst;

    filter2D(src, dst, src.type(), kernel, Point(-1, -1), 0, BORDER_REPLICATE);
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
    double kernel_arr[3][3] = {
        {.25, 0, -.25},
        {0,   0,    0},
        {-.25,0,  .25}
    };
    Mat kernel(3, 3, CV_64F, kernel_arr), dst;

    filter2D(src, dst, src.type(), kernel, Point(-1, -1), 0, BORDER_REPLICATE);
    return dst;
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
    assert(matchAll(flow, weight) && flow.type() == CV_64F);
    
    int rows = flow.rows, cols = flow.cols, offset;
    int step = flow.step / sizeof(double);
    double *pf = (double *)flow.data, *pw = (double *)weight.data;

    lap.create(rows, cols, flow.type());
    lap.setTo(0);    
    double *pl = (double *)lap.data;

    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            offset = r * step + c;
            if (c < cols-1)
                pl[offset] -= pw[offset] * (pf[offset+1]-pf[offset]);
            if (c > 0)
                pl[offset] += pw[offset-1] * (pf[offset]-pf[offset-1]);
            if (r < rows-1)
                pl[offset] -= pw[offset] * (pf[offset+step]-pf[offset]);
            if (r > 0)
                pl[offset] += pw[offset-step] * (pf[offset]-pf[offset-step]);

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
