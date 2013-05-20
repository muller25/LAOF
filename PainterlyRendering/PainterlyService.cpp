#include "PainterlyService.h"
#include "ColorModels.h"
#include "RBF.h"

#include <cstdio>
#include <cmath>
#include <ctime>

#include <highgui.h>
using namespace cv;

#define PA 0.4

#define DEBUG

PainterlyService::PainterlyService()
{
    m_init = false;
    m_width = m_height = 0;

    m_brush_radius.push_back(16);
    m_brush_radius.push_back(8);
    m_brush_radius.push_back(4);
    m_brush_radius.push_back(2);
    
    // 设置渲染参数, Impressionism style
    const int max_stroke_length = 8;
    const int min_stroke_length = 2;
    const int threshold = 20 * 3;
    const float grid_size = 2.0;

    currentStyle.blur_factor = 0.5;
    currentStyle.max_stroke_length = max_stroke_length;
	currentStyle.min_stroke_length = min_stroke_length;
	currentStyle.curvature_filter = 1.0;
	currentStyle.threshold = threshold;
	currentStyle.alpha = 1.0;
    currentStyle.grid_size = grid_size;
	currentStyle.jitter_r = 0.0;
	currentStyle.jitter_g = 0.0;
	currentStyle.jitter_b = 0.0;
	currentStyle.jitter_hue = 0.0;
	currentStyle.jitter_sat = 0.0;
	currentStyle.jitter_val = 0.0;

//    setTexture("brush_map.bmp");
    setTexture();
}

void PainterlyService::setSourceImage(const Mat &src)
{
    assert(src.type() == CV_8UC3);

    if (m_width != src.cols || m_height != src.rows) {
        m_width = src.cols, m_height = src.rows;
    }

    src.copyTo(m_src);
    m_init = true;
    bilateralFilter(m_src, m_reference, 9, 300, 300);
    m_orient = Mat::zeros(m_height, m_width, CV_PERCISION);
    getEdgeMap();
    getStrokeOrientation();

}

void PainterlyService::setTexture(const char *path)
{
    if (NULL == path) m_texture = Mat::ones(100, 100, CV_8U) * 255;
    else m_texture = imread(path, CV_LOAD_IMAGE_GRAYSCALE);
    assert(m_texture.data != NULL);
}

// 提取图像边缘，修改m_edge_map
void PainterlyService::getEdgeMap()
{
    assert(m_init);
    
    Mat smooth, gray;
    cvtColor(m_reference, smooth, CV_BGR2GRAY);
    bitwise_not(smooth, gray);
    
    // 对灰度图像进行边缘检测
    const float edge_dect_thres = 100;

    m_edge_map = Mat::zeros(m_height, m_width, CV_8UC3);
    Canny(gray, smooth,  edge_dect_thres,  edge_dect_thres * 3, 3);
    m_src.copyTo(m_edge_map, smooth);
}

void PainterlyService::getStrokeOrientation()
{
    assert(m_init);

    Mat rbf;
    RBF::rbf_interpolate(rbf, m_src);
    
    // gradient
    Mat gray, gx, gy, abs_gx, abs_gy, gm;
    cvtColor(m_src, gray, CV_BGR2GRAY);
    Sobel(gray, gx, CV_PERCISION, 1, 0);
    Sobel(gray, gy, CV_PERCISION, 0, 1);
    convertScaleAbs(gx, abs_gx);
    convertScaleAbs(gy, abs_gy);
    addWeighted(abs_gx, 0.5, abs_gy, 0.5, 0, gm, CV_PERCISION);
    
    // weighted orientation
    double maxVal, minVal, range;
    minMaxLoc(gm, &minVal, &maxVal);
    range = maxVal - minVal;
    if (fabs(range) <= 1e-6) range = 1;
    
    for (int h = 0; h < m_height; ++h)
        for (int w = 0; w < m_width; ++w)
        {
            PERCISION dx = gx.at<PERCISION>(h, w);
            PERCISION dy = gy.at<PERCISION>(h, w);
            PERCISION lo = atan2(dy, dx) + CV_PI / 2.;
            PERCISION go = rbf.at<PERCISION>(h, w);
            PERCISION weight = (gm.at<PERCISION>(h, w) - minVal) / range;
//            m_orient.at<PERCISION>(h, w) = weight * lo + (1 - weight) * go;
            m_orient.at<PERCISION>(h, w) = lo;
        }
}

void PainterlyService::make_spline_stroke(SplineStroke &spline_stroke, int x0, int y0, int R,
                                          const Mat &canvas, const Mat &coverage)
{
    assert(m_init && coverage.type() == CV_32S);
    
    Vec3i bcolor = m_reference.at<Vec3b>(y0, x0);
    unsigned int r = bcolor[2];
    unsigned int g = bcolor[1];
    unsigned int b = bcolor[0];

    double jr = 0.0;        // jittered red
    double jg = 0.0;        // jittered green       
    double jb = 0.0;        // jittered blue
    double jh = 0.0;        // jittered hue
    double js = 0.0;        // jittered saturation
    double jv = 0.0;        // jittered value

    // convert rgb to hsv
    ColorModels:: rgb_to_hsv((double) r / 256.0, (double) g / 256.0, (double) b / 256.0, &jh, &js, &jv);

    // jitter hsv
    srand(time(NULL));
    jh += currentStyle.jitter_hue * (rand() % 1000 - 500) / 1000.0 * 360.0 ;
    js += currentStyle.jitter_sat * (rand() % 1000 - 500) / 1000.0;
    jv += currentStyle.jitter_val * (rand() % 1000 - 500) / 1000.0;

    // convert hsv back to rgb
    ColorModels::hsv_to_rgb(&jr, &jg, &jb, jh, js, jv);

    // update jittered color
    r = (int) (jr * 256);
    g = (int) (jg * 256);
    b = (int) (jb * 256);

    // jitter rgb
    jr = currentStyle.jitter_r * (rand() % 1000 - 500) / 1000.0;
    jg = currentStyle.jitter_g * (rand() % 1000 - 500) / 1000.0;
    jb = currentStyle.jitter_b * (rand() % 1000 - 500) / 1000.0;

    // update jittered color
    r += (int) (jr * 256);
    g += (int) (jg * 256);
    b += (int) (jb * 256);

    // init spline stroke
    const double angle_threshold = CV_PI / 3.;
    double btheta = m_orient.at<PERCISION>(y0, x0);
    spline_stroke.set(R, btheta, bcolor);

    // add control points to spline stroke
    int x = x0, y = y0;
    for (int k = 0; k < currentStyle.max_stroke_length; k++)
    {
        if (x < 0 || x >= m_width || y < 0 || y >= m_height) break;

        Vec3i rcolor = m_reference.at<Vec3b>(y, x);
        Vec3i ccolor = canvas.at<Vec3b>(y, x);
        Vec3i drc = rcolor - ccolor;
        Vec3i drb = rcolor - bcolor;
        double d1 = norm(drc, NORM_L1);
        double d2 = norm(drb, NORM_L1);
        double theta = m_orient.at<PERCISION>(y, x);
        
        if (k > currentStyle.min_stroke_length && (d1 < d2 || fabs(btheta-theta) > angle_threshold || coverage.at<int>(y, x) >= 1))
            break;

        x += R * cos(theta);
        y += R * sin(theta);
        spline_stroke.add(x, y);
    }
}

/*
void PainterlyService::make_spline_stroke(SplineStroke &spline_stroke,
                                          int x0, int y0, int R,
                                          const Mat &dst, const Mat &ref)
{
    assert(m_init);
    
    uchar *pref = ref.data, *pdst = dst.data;
    int step = ref.step1();
    int color_index = y0 * step + x0 * 3;     
    unsigned int r = pref[color_index + 2];
    unsigned int g = pref[color_index + 1];
    unsigned int b = pref[color_index + 0];

    int x = x0;             // current x coord
    int y = y0;             // current y coord
    double jr = 0.0;        // jittered red
    double jg = 0.0;        // jittered green       
    double jb = 0.0;        // jittered blue
    double jh = 0.0;        // jittered hue
    double js = 0.0;        // jittered saturation
    double jv = 0.0;        // jittered value

    // convert rgb to hsv
    ColorModels:: rgb_to_hsv((double) r / 256.0, (double) g / 256.0, (double) b / 256.0, &jh, &js, &jv);

    // jitter hsv
    srand(time(NULL));
    jh += currentStyle.jitter_hue * (rand() % 1000 - 500) / 1000.0 * 360.0 ;
    js += currentStyle.jitter_sat * (rand() % 1000 - 500) / 1000.0;
    jv += currentStyle.jitter_val * (rand() % 1000 - 500) / 1000.0;

    // convert hsv back to rgb
    ColorModels::hsv_to_rgb(&jr, &jg, &jb, jh, js, jv);

    // update jittered color
    r = (int) (jr * 256);
    g = (int) (jg * 256);
    b = (int) (jb * 256);

    // jitter rgb
    jr = currentStyle.jitter_r * (rand() % 1000 - 500) / 1000.0;
    jg = currentStyle.jitter_g * (rand() % 1000 - 500) / 1000.0;
    jb = currentStyle.jitter_b * (rand() % 1000 - 500) / 1000.0;

    // update jittered color
    r += (int) (jr * 256);
    g += (int) (jg * 256);
    b += (int) (jb * 256);

    // init spline stroke
    spline_stroke.set(R, r, g, b, x0, y0);

    // add control point to spline stroke
    spline_stroke.add(x0, y0);

    // add control points to spline stroke
    int offset, tmp, dist, grid_step = currentStyle.grid_size * R;
    for (int k = 1; k <= currentStyle.max_stroke_length; k++)
    {
        x = x+R*cos(m_grad_orient[y*m_width+x]+PI/2);
        y = y+R*sin(m_grad_orient[y*m_width+x]+PI/2);

        if (x < 0 || x >= m_width || y < 0 || y >= m_height) break;

        // 计算 area_error
        int area_error = 0;
        for (int ii = -(grid_step / 2); ii <= grid_step / 2; ii++)
        {
            for (int jj = -(grid_step / 2); jj <= grid_step / 2; jj++)
            {
                int iii = x + ii;
                int jjj = y + jj;

                //边界控制
                if (iii < 0) iii = 0;
                if (iii >= m_width) iii = m_width - 1;
                if (jjj < 0) jjj = 0;
                if (jjj >= m_height) jjj = m_height - 1;

                //计算(x,y)区域最大差值点及差值
                dist = 0;
                offset = jjj * step + iii * 3;
                for (int k = 0; k < 3; ++k)
                {
                    tmp = pdst[offset+k] - pref[offset+k];
                    dist = tmp * tmp;
                }
                area_error += sqrt((float)dist);
            }
        }

        // 当区域平均差值大于阈值时，绘制笔刷
        area_error /= (grid_step * grid_step);
        if (k >= currentStyle.min_stroke_length && area_error <= currentStyle.threshold) break;

        spline_stroke.add(x, y);
    }
}
*/
// generate brush strokes
// strokes_queue - stroke queue (out)
// ref - smoothed src image
// R - brush stroke radius
void PainterlyService::generate_strokes(list<SplineStroke> &strokes_queue, int R, const Mat &canvas)
{
    assert(m_init && canvas.type() == CV_8UC3);

    strokes_queue.clear();
    int grid_step = currentStyle.grid_size * R;
    int area = grid_step * grid_step;
    Mat coverage = Mat::zeros(m_height, m_width, CV_32S);
    clock_t start = clock();
    for (int j = m_height - 1; j >= 0; j -= grid_step)
    {
        for (int i = 0; i < m_width; i += grid_step)
        {
            // 计算 area_error
            double area_error = 0, max_error = -1;
            int index_i = i, index_j = j;
            for (int ii = -(grid_step / 2); ii <= grid_step / 2; ii++)
                for (int jj = -(grid_step / 2); jj <= grid_step / 2; jj++)
                {
                    // 计算(x,y)区域最大差值点及差值
                    int iii = std::min(std::max(i + ii, 0), m_width-1);
                    int jjj = std::min(std::max(j + jj, 0), m_height-1);
                    Vec3i ccolor = canvas.at<Vec3b>(jjj, iii);
                    Vec3i rcolor = m_reference.at<Vec3b>(jjj, iii);
                    Vec3i diff = ccolor - rcolor;
                    double current_error = norm(diff);
                    
                    if (max_error < current_error)
                    {
                        max_error = current_error;
                        index_i = iii;
                        index_j = jjj;
                    }

                    area_error += current_error;
                }

            // 当区域平均差值大于阈值时，绘制笔刷
            area_error /= area;
            if (area_error <= currentStyle.threshold) continue;

            // 生成CubicBSpline曲线
            SplineStroke stroke, cbstroke;
            make_spline_stroke(stroke, index_i, index_j, R, canvas, coverage);

            Vec3b color = stroke.getColor();
            double angle = stroke.getAngle();
            cbstroke.set(R, angle, color);

            const double animate_time = 1.0;
            Point2d newPoint, prePoint(-1, -1);
            double time_step = animate_time / (stroke.nPoints() * R);
            for (double t = 0; t <= animate_time; t += time_step)
            { 
                stroke.cubic_b_spline(newPoint, t);
                newPoint.x = std::max(std::min(newPoint.x, (double)(m_width-1)), 0.);
                newPoint.y = std::max(std::min(newPoint.y, (double)(m_height-1)), 0.);
                if (newPoint == prePoint) continue;

                // coverage map               
                int nx = cvRound(newPoint.x), ny = cvRound(newPoint.y);
                for (int h = -R; h <= R; ++h)
                    for (int w = -R; w <= R; ++w)
                    {
                        int hh = std::max(std::min(ny+h, m_height-1), 0);
                        int ww = std::max(std::min(nx+w, m_width-1), 0);
                        coverage.at<int>(hh, ww) += 1;
                    }

                cbstroke.add(nx, ny);
                newPoint = prePoint;
            }
            
            strokes_queue.push_back(cbstroke);
        }
    }

#ifdef DEBUG
    imshow("coverage", coverage * 255);
    waitKey(0);
#endif
    
    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("time of running stroke construction steps: %.2f\n", duration);
}
    
// Find areas that we need to paint and draw strokes in them
// ref - smoothed src image
// R - brush stroke radius
// strokes_queue - stroke queue
void PainterlyService::paint_layer(Mat &canvas, const list<SplineStroke> &strokes_queue, int layerId)
{
    assert(m_init && canvas.type() == CV_8UC3);
    if (canvas.empty()) canvas = Mat::zeros(m_height, m_width, m_src.type());

    int R = m_brush_radius[layerId];
    Mat brush_tex, bg;
    resize(m_texture, brush_tex, Size(R*2+1, R*2+1), 0, 0, INTER_CUBIC);
    canvas.copyTo(bg);
    
    // CubicBSpline版本, 得到笔刷的control points
    printf("strokes size: %d\n", strokes_queue.size());

    clock_t start = clock();
    for (list<SplineStroke>::const_iterator iter = strokes_queue.begin(); iter != strokes_queue.end(); iter++)
    {
        if(iter->nPoints() <= 0) continue;

        double alpha = iter->getAlpha();
        Point point;
        Vec3b bcolor = iter->getColor();
        for (int i = 0; i < iter->nPoints(); ++i)
        {
            point = iter->get(i);
            Vec3b bgcolor = bg.at<Vec3b>(point.y, point.x);
            Vec3b ccolor = canvas.at<Vec3b>(point.y, point.x);
            Vec3b color = bgcolor * (1-alpha) + (bcolor * (1-PA) + ccolor * PA) * alpha;

            circle(canvas, point, R, Scalar(color[0], color[1], color[2]), -1);

#ifdef DEBUG
//            imshow("painter", canvas);
//            waitKey(10);
#endif
            // for(int hh = -R; hh <= R; hh++)
            //     for(int ww= -R; ww <= R; ww++){
            //         bii = std::min(std::max(point.y + hh, 0), m_height-1);
            //         bjj = std::min(std::max(point.x + ww, 0), m_width-1);
            //         offset = bii * step + bjj * 3;
            //         alpha = iter->nAlpha();//(1 - PA) * iter->nAlpha();
            //         r = pref[offset+2]*(1-alpha) + iter->ColorR()*alpha;
            //         g = pref[offset+1]*(1-alpha) + iter->ColorG()*alpha;
            //         b = pref[offset]*(1-alpha) + iter->ColorB()*alpha;

            //         p_xx = newPoint.x + ww;
            //         p_yy = newPoint.y + hh;
            //         if(!(p_xx >= 0 && p_yy >=0 && p_xx < m_width && p_yy < m_height)) continue;

            //         int ii = hh + R, jj = ww + R;
            //         int I = brush_tex.at<uchar>(ii, jj);
            //         int Index = p_yy * step + p_xx * 3;

            //         // 附加纹理
            //         pdst[Index+2] = r*(I/255.0) + pdst[Index+2]*(1- I/255.0);
            //         pdst[Index+1] = g*(I/255.0) + pdst[Index+1]*(1- I/255.0);
            //         pdst[Index+0] = b*(I/255.0) + pdst[Index+0]*(1- I/255.0);
            //         // pdst[Index+2] = r;
            //         // pdst[Index+1] = g;
            //         // pdst[Index] = b;
                    
            //         // 计算高度图
            //         int height_map_index = p_yy * m_width + p_xx;
            //         m_count_pass[height_map_index]++;
            //         m_sum_pass[height_map_index] += I;
            //         m_height_maps.at<uchar>(p_yy, p_xx) = (int)((m_sum_pass[height_map_index])/(m_count_pass[height_map_index]));
            // }
        }
    }

    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("time of running CubicBSpline drawing steps: %.2f\n", duration);
}

void PainterlyService::paint_layer(Mat &canvas, const list<SplineStroke> *strokes_queue, int nlayer)
{
    assert(m_init && (nlayer <= m_brush_radius.size()));

    for(int i = 0; i < nlayer; ++i)
    {
        if(m_brush_radius[i] == 0) continue;
        paint_layer(canvas, strokes_queue[i], i);
    }
}

// image rendering
void PainterlyService::render(Mat &canvas)
{
    assert(m_init);

    int nlayer = m_brush_radius.size();
    list<SplineStroke> strokes_queue[nlayer];
    render(canvas, strokes_queue, nlayer);

#ifdef DEBUG
    imwrite("edge_map.jpg", m_edge_map);
    imwrite("reference.jpg", m_reference);
    Mat show;
    RBF::plot(show, m_orient, m_src, 8);
    imwrite("orient.jpg", show);
#endif
}

// image rendering
void PainterlyService::render(Mat &canvas, list<SplineStroke> *strokes_queue, int nlayer)
{
    assert(m_init && (nlayer <= m_brush_radius.size()));
    if (canvas.empty()) canvas = Mat::zeros(m_height, m_width, m_src.type());
    
    clock_t begin, finish, start, end;
    double duration;
    Mat bg = Mat::zeros(m_height, m_width, m_src.type());
    int R;
    
    begin = clock();
    for(int i = 0; i < nlayer; ++i)
    {
        R = m_brush_radius[i];
        if(R == 0) continue;
        
        printf("====== Painting layer with stroke size: %d ======\n", R);
        start = clock();
        render(canvas, strokes_queue[i], i);
        finish = clock();
        duration = (double)(finish - start) / CLOCKS_PER_SEC;
        printf("time of running paint_layer(): %.2f\n", duration);
    }
    
    end = clock();
    duration = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time %.2f !!!!!\n", duration);
}

void PainterlyService::render(Mat &canvas, list<SplineStroke> &strokes_queue, int layerId)
{
    generate_strokes(strokes_queue, m_brush_radius[layerId], canvas);
    paint_layer(canvas, strokes_queue, layerId);
}
