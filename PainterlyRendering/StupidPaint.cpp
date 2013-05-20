#include "StupidPaint.h"
#include "RBF.h"

#include <highgui.h>
using namespace cv;

Mat StupidPaint::m_height_map;
Mat StupidPaint::m_alpha_map;

void rotate(Mat &dst, const Mat &src, double radian)
{
    int height = src.rows, width = src.cols;
    Point2f center((width + 1.) / 2., (height + 1.) / 2.);
    double degree = radian * 180 / CV_PI;
    Mat rotate = getRotationMatrix2D(center, degree, 1);
    double sina = fabs(sin(radian)), cosa = fabs(cos(radian));
    int nwidth = height * sina + width * cosa, nheight = width * sina + height * cosa;
    Size dsize(nwidth, nheight);
    rotate.at<double>(0, 2) += (nwidth - width) / 2.;
    rotate.at<double>(1, 2) += (nheight - height) / 2.;
    warpAffine(src, dst, rotate, dsize, INTER_LINEAR, BORDER_CONSTANT, Scalar::all(255));
}

void StupidPaint::loadTexture()
{
    m_height_map = imread("height.bmp", CV_LOAD_IMAGE_GRAYSCALE);
    m_alpha_map = imread("alpha.bmp", CV_LOAD_IMAGE_GRAYSCALE);
}

void StupidPaint::generate_strokes(vector<Brush> &strokes, const Mat &src, int radius)
{
    assert(src.type() == CV_8UC3);
    strokes.clear();

    // resize brush texture
    int bwidth = radius * 4 + 1, bheight = radius * 2 + 1;
    Size bsize(bwidth, bheight);
    Mat alpha_map;
    resize(m_alpha_map, alpha_map, bsize);
    
    // calc global orientation
    Mat smooth, im, orient;
    int ksize = 2 * radius + 1;
    double factor = 0.75 / radius;
    GaussianBlur(src, smooth, Size(ksize, ksize), radius, radius);
    resize(smooth, im, Size(0, 0), factor, factor);
    stroke_orient(orient, im);

    // extract strokes
    Mat mask;
    double scale = 1. / factor, len = scale / 2.;
    int width = im.cols, height = im.rows;
    for (int h = 0; h < height; ++h)
        for (int w = 0; w < width; ++w)
        {
            PERCISION theta = orient.at<PERCISION>(h, w);
//            rotate(mask, alpha_map, theta);

            int nw = w * scale + len;
            int nh = h * scale + len;
//            Vec3b color = extract_color(src, mask, nw, nh);
            Vec3b color = extract_color(src, nw, nh, radius);
            Brush brush(Point(nw, nh), color, bsize, theta);
            strokes.push_back(brush);
        }
}

/*
  mask: 0 for no taken, otherwise for taken
  x, y: center point
*/
Vec3b StupidPaint::extract_color(const Mat &src, const Mat &mask, int x, int y)
{
    assert(src.type() == CV_8UC3 && mask.type() == CV_8U);
    
    Vec3f sum(0, 0, 0);
    Vec3b color;
    int width = src.cols, height = src.rows, count = 0;
    int wsize = (width + 1) / 2., hsize = (height + 1) / 2.;
    for (int h = -hsize; h <= hsize; ++h)
        for (int w = -wsize; w <= wsize; ++w)
        {
            int nw = std::min(std::max(x + w, 0), width-1);
            int nh = std::min(std::max(y + h, 0), height-1);
            if (mask.at<uchar>(nh, nw) == 255) continue; // not used

            sum += src.at<Vec3b>(nh, nw);;
            ++count;
        }

    if (count <= 0) count = 1;
    color[0] = sum[0] / count, color[1] = sum[1] / count, color[2] = sum[2] / count;

    return color;
}

Vec3b StupidPaint::extract_color(const Mat &src, int x, int y, int radius)
{
    int width = src.cols, height = src.rows;
    int r2 = radius * radius, count = 0;
    Vec3f sum(0, 0, 0);
    Vec3b color;
    for (int h = -radius; h <= radius; ++h)
    {
        for (int w = -radius; w <= radius; ++w)
        {
            if (h * h + w * w > r2) continue;
            int hh = std::min(std::max(y + h, 0), height-1);
            int ww = std::min(std::max(x + w, 0), width-1);
            sum += src.at<Vec3b>(hh, ww);
            count++;
        }
    }

    if (count <= 0) count = 1;
    color[0] = sum[0] / count, color[1] = sum[1] / count, color[2] = sum[2] / count;

    return color;
}

void StupidPaint::stroke_orient(Mat &orient, const Mat &src)
{
    Mat gx, gy;
    vector<Point> centers;
    RBF::rbf_interpolate(gx, gy, centers, src);

    // plot
    Mat rbfres;
    double factor = 0.75 / 8;
    RBF::plot(rbfres, centers, gx, gy, src, 1 / factor);
    imshow("rbf", rbfres);

//    RBF::rbf_interpolate(gx, gy, src);

    int width = src.cols, height = src.rows;
    orient.create(height, width, CV_PERCISION);
    for (int h = 0; h < height; ++h)
        for (int w = 0; w < width; ++w)
        {
            PERCISION dx = gx.at<PERCISION>(h, w);
            PERCISION dy = gy.at<PERCISION>(h, w);
            orient.at<PERCISION>(h, w) = atan2(dy, dx) + CV_PI / 2.;
        }
}

void StupidPaint::strokes_placement(Mat &dst, const Mat &src, const vector<Brush> &strokes)
{
    assert(src.type() == CV_8UC3 && !strokes.empty());

    printf("strokes: %d\n", strokes.size());
    
    int width = src.cols, height = src.rows;    

    // resize brush texture
    Size bsize = strokes[0].getSize();
    Mat alpha_map;
    resize(m_alpha_map, alpha_map, bsize);
    
    // paint
    Mat amask;
    float halpha = 0.9;
    dst = Mat::zeros(height, width, src.type());
    src.copyTo(dst);
    for (size_t i = 0; i < strokes.size(); ++i)
    {
        double angle = strokes[i].getAngle();
        Point center = strokes[i].getCenter();
        Vec3b bcolor = strokes[i].getColor();

        rotate(amask, alpha_map, angle);
        int mwidth = amask.cols, mheight = amask.rows;
        int wsize = mwidth / 2., hsize = mheight / 2.;
        for (int h = -hsize; h <= hsize; ++h)
            for (int w = -wsize; w <= wsize; ++w)
            {
                int nh = std::max(0, std::min(h + center.y, height-1));
                int nw = std::max(0, std::min(w + center.x, width-1));
                int hh = std::min(h + hsize, mheight-1);
                int ww = std::min(w + wsize, mwidth-1);

                uchar val = amask.at<uchar>(hh, ww);
                if (val >= 128) continue;

                float aalpha = 0;//val / 255.;
                Vec3b dcolor = dst.at<Vec3b>(nh, nw);
                dst.at<Vec3b>(nh, nw) = bcolor * (1-aalpha) * halpha + dcolor * (1-halpha);
            }
    }
}

void StupidPaint::paint(Mat &dst, const Mat &src)
{
    vector<Brush> strokes;
    printf("generating strokes...\n");
    generate_strokes(strokes, src, 8);

    printf("placing strokes...\n");
    strokes_placement(dst, src, strokes);
}