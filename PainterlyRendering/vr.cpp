#include "PainterlyService.h"
#include "SplineStroke.h"

#include <cv.h>
#include <highgui.h>
using namespace cv;

#include <iostream>
#include <list>
using namespace std;

#define DEBUG

bool imreadf(Mat &m, const char *filename);
void propagate(list<SplineStroke> *strokes_queue, int nlayer,
               const Mat &u, const Mat &v);
void removeStrokes(list<SplineStroke> *strokes_queue, int nlayer, const Mat &src);
void addStrokes(list<SplineStroke> *strokes_queue, const list<SplineStroke> *pre_strokes_queue,
                int nlayer, Mat &src);
void getCoverage(Mat &coverage, const list<SplineStroke> &strokes_queue);

int main(int argc, char *argv[])
{
    if (argc != 6){
        cout << "usage: ./vr image flowDir outDir start end" << endl;
        return 1;
    }

    char im[256], flow[256], out[256], buf[256];
    int frameStart, frameEnd;

    memset(im, 0, sizeof(im));
    memset(flow, 0, sizeof(flow));
    memset(out, 0, sizeof(out));
    strcpy(im, argv[1]);
    strcpy(flow, argv[2]);
    strcat(flow, "%s%03d.yml");
    strcpy(out, argv[3]);
    strcat(out, "%s%03d.jpg");
    frameStart = atoi(argv[4]);
    frameEnd = atoi(argv[5]);

    Mat canvas, src, preRender, u, v;
    PainterlyService ps;
    const size_t nlayer = ps.nLayers();
    const int preRenderInterval = 10;
    int preload, frame;
    list<SplineStroke> strokes_queue[nlayer], pre_strokes_queue[nlayer];
    
    preload = std::min(frameStart + preRenderInterval, frameEnd);
    sprintf(buf, im, preload);
    cout << "pre-rendering image " << buf << endl;

    src = imread(buf);
    ps.setSourceImage(src);
    ps.generate_strokes(pre_strokes_queue, nlayer);
    src.copyTo(preRender);
    ps.paint_layer(preRender, pre_strokes_queue, nlayer);
    ps.fixEdges(preRender);
    sprintf(buf, out, "pre", preload);
    imwrite(buf, preRender);

    sprintf(buf, im, frameStart);
    cout << "loading image " << buf << endl;

    src = imread(buf);
    ps.setSourceImage(src);
    ps.generate_strokes(strokes_queue, nlayer);
    preRender.copyTo(canvas);
    ps.paint_layer(canvas, strokes_queue, nlayer);
    ps.fixEdges(canvas);
    sprintf(buf, out, "out", frameStart);
    imwrite(buf, canvas);

    for (int count = 1; count <= (frameEnd - frameStart); ++count)
    {
        frame = count + frameStart;
        sprintf(buf, im, frame);
        cout << "loading image " << buf << endl;
        src = imread(buf);
        ps.setSourceImage(src);
        sprintf(buf, flow, "u", frame-1);
        imreadf(u, buf);
        sprintf(buf, flow, "v", frame-1);
        imreadf(v, buf);

#ifdef DEBUG
            for (size_t l = 0; l < nlayer; ++l)
            {
                preRender.copyTo(canvas);
                ps.paint_layer(canvas, strokes_queue[l], l);
                sprintf(buf, out, "propagate", frame*100+l);
                imwrite(buf, canvas);
            }
#endif
        propagate(strokes_queue, nlayer, u, v);

#ifdef DEBUG
            for (size_t l = 0; l < nlayer; ++l)
            {
                preRender.copyTo(canvas);
                ps.paint_layer(canvas, strokes_queue[l], l);
                sprintf(buf, out, "remove", frame*100+l);
                imwrite(buf, canvas);
            }
#endif
        removeStrokes(strokes_queue, nlayer, src);

        if (count % preRenderInterval == 0)
        {
#ifdef DEBUG
            for (size_t l = 0; l < nlayer; ++l)
            {
                preRender.copyTo(canvas);
                ps.paint_layer(canvas, strokes_queue[l], l);
                sprintf(buf, out, "before-add", frame*100+l);
                imwrite(buf, canvas);
            }
#endif
//            addStrokes(strokes_queue, pre_strokes_queue, nlayer, src);
#ifdef DEBUG
            for (size_t l = 0; l < nlayer; ++l)
            {
                preRender.copyTo(canvas);
                ps.paint_layer(canvas, strokes_queue[l], l);
                sprintf(buf, out, "after-add", frame*100+l);
                imwrite(buf, canvas);
            }
#endif
            
            preload = std::min(frameStart + count + preRenderInterval, frameEnd);
            sprintf(buf, im, preload);
            cout << "pre-rendering image " << buf << endl;

            src = imread(buf);
            ps.setSourceImage(src);
            ps.generate_strokes(pre_strokes_queue, nlayer);
            src.copyTo(preRender);
            ps.paint_layer(preRender, pre_strokes_queue, nlayer);
            ps.fixEdges(preRender);
            sprintf(buf, out, "pre", preload);
            imwrite(buf, preRender);
        }

        ps.setSourceImage(src);
        preRender.copyTo(canvas);
        ps.paint_layer(canvas, strokes_queue, nlayer);
        ps.fixEdges(canvas);
        sprintf(buf, out, "out", frame);
        imwrite(buf, canvas);
    }
    
	return 0;
}

bool imreadf(Mat &m, const char *filename)
{
    printf("reading file %s... ", filename);
    
    FileStorage fs(filename, cv::FileStorage::READ);

    if (!fs.isOpened())
    {
        fs.release();
        printf("failed\n");
        return false;
    }
    
    fs["matrix"] >> m;
    fs.release();
    printf("done\n");
    return true;
}

void propagate(list<SplineStroke> *strokes_queue, int nlayer,
               const Mat &u, const Mat &v)
{
    assert(u.type() == CV_64F);

    cout << "propagating strokes...";

    int offset, step = u.step1(), width = u.cols, height = u.rows;
    double *pu = (double*)u.data, *pv = (double*)v.data;
    
    for (int i = 0; i < nlayer; ++i)
    {
        if (strokes_queue[i].empty()) continue;

        for (list<SplineStroke>::iterator iter = strokes_queue[i].begin(); iter != strokes_queue[i].end(); ++iter)
        {
            for (int k = 0; k < iter->nPoints(); ++k)
            {
                Point p = iter->get(k);
                if (p.y >= height || p.y < 0 || p.x >= width || p.x < 0) continue;

                offset = p.y * step + p.x;
                p.x += cvRound(pu[offset]);
                p.y += cvRound(pv[offset]);
                iter->setControlPoint(k, p);
            }
        }
    }

    cout << "done" << endl;
}

void removeStrokes(list<SplineStroke> *strokes_queue, int nlayer, const Mat &src)
{
    assert(src.type() == CV_8UC3);
    cout << "removing strokes...";
    
    const int threshold = 20 * 3;
    const double fadeStep = -0.1;
    int width = src.cols, height = src.rows;

    for (int i = 0; i < nlayer; ++i)
    {
        if (strokes_queue[i].empty()) continue;

        list<SplineStroke>::iterator iter = strokes_queue[i].begin();
        int R = iter->getRadius();
        Vec3i bcolor = iter->getColor();
        while (iter != strokes_queue[i].end())
        {
            int pointDiff = 0;
            for (int k = 0; k < iter->nPoints(); ++k)
            {
                // color diff of brush covered areas
                Point p = iter->get(k);
                int count = 0, diff = 0;
                for (int hh = -R; hh <= R; ++hh)
                    for (int ww = -R; ww <= R; ++ww)
                    {
                        int hhh = p.y + hh;
                        int www = p.x + ww;
                        if (hhh >= height || www >= width || hhh < 0 || www < 0) continue;

                        Vec3i scolor = src.at<Vec3b>(hhh, www);
                        Vec3i diffcolor = scolor - bcolor;
                        diff += norm(diffcolor, NORM_L1);
                        count++;
                    }

                if (count <= 0) count = 1;
                diff /= count;
                if (diff >= threshold) ++pointDiff;
            }

            // more than half of control points are diff to underneath color
            if (pointDiff >= 0.5 * iter->nPoints()) iter->changeOpacity(fadeStep);
            if (iter->isTransparent()) iter = strokes_queue[i].erase(iter);
            else iter++;
        }
    }

    cout << "done" << endl;
}

void addStrokes(list<SplineStroke> *strokes_queue, const list<SplineStroke> *pre_strokes_queue,
                int nlayer, Mat &src)
{
    cout << "adding strokes..." << endl;

    int width = src.cols, height = src.rows;
    Mat coverage(height, width, CV_32S);
    for (int i = 0; i < nlayer; ++i)
    {
        if (pre_strokes_queue[i].empty()) continue;

        coverage.setTo(0);
        getCoverage(coverage, strokes_queue[i]);
        for (list<SplineStroke>::const_iterator iter = pre_strokes_queue[i].begin(); iter != pre_strokes_queue[i].end(); ++iter)
        {
            // 获取笔刷区域
            int R = iter->getRadius();
            Point center = iter->get(0);
            int brushLen = R * iter->nPoints();
            int maxBrushLen = 2 * brushLen + 1;
            Mat brush = Mat::zeros(maxBrushLen, maxBrushLen, CV_8U);
            for (int k = 0; k < iter->nPoints(); ++k)
            {
                Point point = iter->get(k);
                for (int h = -R; h <= R; ++h)
                    for (int w = -R; w <= R; ++w)
                    {
                        int hh = std::min(std::max(point.y-center.y+brushLen+h, 0), maxBrushLen-1);
                        int ww = std::min(std::max(point.x-center.x+brushLen+w, 0), maxBrushLen-1);
                        brush.at<uchar>(hh, ww) = 1;
                    }
            }

            // 获取笔刷覆盖区域
            int area = 0, covered = 0;
            for (int h = -brushLen; h <= brushLen; ++h)
                for (int w = -brushLen; w <= brushLen; ++w)
                {
                    if (brush.at<uchar>(h+brushLen, w+brushLen) == 0) continue;
                    ++area;

                    int hh = std::min(std::max(center.y+h, 0), height-1);
                    int ww = std::min(std::max(center.x+w, 0), width-1);
                    if (coverage.at<int>(hh, ww) > 0) ++covered;
                }

            // 如果笔刷覆盖区域有超过一半没有渲染，则添加笔刷
            if (covered < 0.5 * area)
                strokes_queue[i].push_front(*iter);
        }
    }

    cout << "done" << endl;
}

void getCoverage(Mat &coverage, const list<SplineStroke> &strokes_queue)
{
    assert(coverage.data != NULL && coverage.type() == CV_32S);

    coverage.setTo(0);
    int width = coverage.cols, height = coverage.rows;
    for (list<SplineStroke>::const_iterator iter = strokes_queue.begin(); iter != strokes_queue.end(); iter++)
    {
        int R = iter->getRadius();
        for (int i = 0; i < iter->nPoints(); ++i)
        {
            Point point = iter->get(i);
            for(int hh = -R; hh <= R; hh++)
                for(int ww = -R; ww <= R; ww++){
                    int bhh = point.y + hh, bww = point.x + ww;
                    if (bhh < 0 || bhh >= height || bww >= width || bww < 0) continue;

                    coverage.at<int>(bhh, bww) += 1;
                }
        }
    }
}
