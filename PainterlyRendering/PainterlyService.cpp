#include "PainterlyService.h"

#include "CubicBSpline.h"
#include "ColorModels.h"

#include <iostream>
#include <cmath>
#include <ctime>

#define PA 0.4

PainterlyService::PainterlyService()
{
    m_init = false;
    m_count_pass = m_sum_pass = m_grad_orient = NULL;
     
    // 设置笔刷大小，不绘制 R=32
    int radius = 1;
    memset(m_brush_radius, 0, sizeof(int) * m_nlayer);
    for (int i = m_nlayer-1; i >= 1; i--){
        m_brush_radius[i] = radius;
        radius <<= 1;
    }
    
    // 设置渲染参数, Impressionism style
    const int max_stroke_length = 16;
    const int min_stroke_length = 2;
    const int threshold = 46;
    const float grid_size = 2.0;

    currentStyle.num_layers = m_nlayer;
    currentStyle.brush_radius = m_brush_radius;
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

    setTexture();
}

void PainterlyService::clear()
{
    if (m_count_pass != NULL) delete []m_count_pass;
    if (m_sum_pass != NULL) delete []m_sum_pass;
    if (m_grad_orient != NULL) delete []m_grad_orient;

    m_count_pass = m_sum_pass = m_grad_orient = NULL;
    m_init = false;
}

void PainterlyService::setSourceImage(const Mat &src)
{
    assert(src.channels() == 3);

    clear();
    m_width = src.cols, m_height = src.rows;
    src.copyTo(m_src);
    m_dst.create(m_height, m_width, CV_8UC3);
    m_height_maps.create(m_height, m_width, CV_8U);
    
    int size = m_width * m_height;
    m_count_pass = new int[size];
    m_sum_pass = new int[size];
    m_grad_orient = new double[size];

    memset(count_pass, 0, size * sizeof(int));
    memset(sum_pass, 0, size * sizeof(int));

    getEdgeMap();    
    getStrokeOrientation();
    m_init = true;
}

void PainterlyService::setTexture(const char *path)
{
    if (path = NULL) m_texture.create(100, 100 CV_8U, Scalar::all(255));
    else m_texture = imread(path, CV_LOAD_IMAGE_GRAYSCALE);
}

// 提取图像边缘，修改m_edge_map
void PainterlyService::getEdgeMap()
{
    Mat gray, smooth;

    cvtColor(m_src, gray, CV_BGR2GRAY);    
    blur(gray, smooth, Size(3, 3));
    bitwise_not(smooth, gray);
    
    // 对灰度图像进行边缘检测
    const float edge_dect_thres = 100;

    m_edge_map = Mat::zeros(m_height, m_width, CV_8UC3);
    Canny(gray, smooth,  edge_dect_thres,  edge_dect_thres * 3, 3);
    m_src.copyTo(m_edge_map, smooth);
}

void PainterlyService::getStrokeOrientation()
{
    assert(m_grad_orient != NULL);

    Mat gray, sobel_x, sobel_y, sobel8_x, sobel8_y;
    cvtColor(m_src, gray, CV_BGR2GRAY);

	Sobel(gray, sobel_x, CV_16S, 1, 0, 3);//x方向求导
	Sobel(gray, sobel_y, CV_16S, 0, 1, 3);//y方向求导
	convertScaleAbs(sobel_x, sobel8_x);
	convertScaleAbs(sobel_y, sobel8_y);

    for(int i=0; i< m_height; i++){
        for(int j=0; j < m_width; j++){
            double gy = (double)sobel8_y.at<uchar>(i, j);
            double gx = (double)sobel8_x.at<uchar>(i, j);
            m_grad_orient[i*m_width+j] = atan2(gy, gx);
        }
    }
}

void PainterlyService::make_spline_stroke(SplineStroke &spline_stroke,
                                          int x0, int y0, int R, const Mat &ref)
{
    assert(m_init);
    
    uchar *pref = ref.data, *pdst = dst.data;
    int step = ref.step1();
    int color_index = y0 * step + x0 * 3;     
    unsigned int r = (uchar)pref[color_index + 2];
    unsigned int g = (uchar)pref[color_index + 1];
    unsigned int b = (uchar)pref[color_index + 0];

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
    for (int k = 1; k <= currentStyle.max_stroke_length; k++)
    {
        if (x < 0 || x >= m_width || y < 0 || y >= m_height) break;
        
        int dif_index = (y * step + x * 3) ;
        int d1 = 0, d2 = -r - g - b;
        for (int c = 0; c < 3; ++c)
        {
            d1 += pref[dif_index + c] - pdst[dif_index + c];
            d2 += pref[dif_index + c];
        }
        
        if (k > currentStyle.min_stroke_length && (abs(d1) < abs(d2))) break;

        x = x+R*cos(m_grad_orient[y*m_width+x]+PI/2);
        y = y+R*sin(m_grad_orient[y*m_width+x]+PI/2);

        spline_stroke.add(x, y);
    }
}

// Difference Image - Calculates the difference between
// diff - strores result (out)
// ref - smoothed src image
// R - brush stroke radius
void PainterlyService::difference_image(int *diff, const Mat &ref, int R)
{
    assert(m_dst.type() == ref.type() && diff != NULL && m_init);
    
    uchar *pdst = m_dst.data, *pref = m_ref.data;
    int step = ref.step1();
    int dif_index, img_index, dif_temp;

    // 计算单个点差异
    for (int j = 0; j < m_height; j++)
    {
        for (int i = 0; i < m_width; i++)
        {
            dif_index = j * m_width + i;
            img_index = j * step + i * 3;
            diff[dif_index] = 0;
            for (int k = 0; k < m_src.channels(); ++k){
                dif_temp = pdst[img_index + k] - pref[img_index + k];
                diff[dif_index] +=  dif_temp * dif_temp;
            }

            diff[dif_index] = sqrt((float)diff[dif_index]);
        }
    }
}

// generate brush strokes
// strokes_queue - stroke queue (out)
// ref - smoothed src image
// R - brush stroke radius
void PainterlyService::generate_strokes(vector<SplineStroke> &strokes_queue,
                                        const Mat &ref, int R)
{
    assert(m_init);
    
    int grid_step = currentStyle.grid_size * R; 
    clock_t start, finish;
    double duration;
    
    start = clock();

    int *diff = new int[m_width * m_height];
    difference_image(diff, ref, R);

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    cout << "time of running difference_image(): " << duration <<endl;

    strokes_queue.clear();
    start = clock();
    for (int j = height - 1; j >= 0; j -= grid_step)
    {
        for (int i = 0; i < width; i += grid_step)
        {
            int area_error = 0.0;
            int current_error = -1;
            int index_i = i, index_j = j;

            // 计算 area_error
            for (int ii = -(grid_step / 2); ii <= grid_step / 2; ii++)
            {
                for (int jj = -(grid_step / 2); jj <= grid_step / 2; jj++)
                {
                    int iii = i + ii;
                    int jjj = j + jj;

                    //边界控制
                    if (iii < 0) iii = 0;
                    if (iii >= m_width) iii = m_width - 1;
                    if (jjj < 0) jjj = 0;
                    if (jjj >= m_height) jjj = m_height - 1;

                    //计算(x,y)区域最大差值点及差值
                    if (current_error <= diff[jjj * width + iii])
                    {
                        current_error = diff[jjj * width + iii];
                        index_i = iii;
                        index_j = jjj;
                    }

                    area_error += diff[jjj * width + iii];
                }
            }

            // 当区域平均差值大于阈值时，绘制笔刷
            area_error /= (grid_step * grid_step);
            if (area_error > currentStyle.threshold)
            {
                SplineStroke stroke;
                make_spline_stroke(stroke, index_i, index_j, R, ref);
                strokes_queue.push_back(stroke);
            }
        }
    }

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    cout << "time of running stroke construction steps: " << duration <<endl;

    delete []diff;
}

// Find areas that we need to paint and draw strokes in them
// ref - smoothed src image
// R - brush stroke radius
// strokes_queue - stroke queue
void PainterlyService::paint_layer(const Mat &ref, int R,
                                   const vector<SplineStroke> &strokes_queue)
{
    assert(m_init);
     
    CubicBSpline CB;
    int step = ref.step1();
    uchar *pdst = m_dst.data;
    Mat brush_tex(R*2+1, R*2+1, m_texture.type());
    
    resize(m_texture, brush_tex, brush_tex.size(), 0, 0, INTER_CUBIC);
    
    // CubicBSpline版本, 得到笔刷的control points
    cout << "strokes size: " << strokes_queue.size() << endl;

    clock_t start = clock();
    for(int j = 0; j < strokes_queue.size(); ++j)
    {	 	
        const SplineStroke stroke = strokes_queue[j];
        if(stroke.nPoints() <= 0) continue;
        
        // 制定brush中每bristle确定颜色
        int bii, bjj, cbii, cbjj;
        const Point startPoint = stroke.getStartPoint();

        // 绘制CubicBSpline曲线
        Point2d newPoint;
        uchar r, g, b;
        int p_xx, p_yy;
        
        // t的细分程度=系数*控制点个数*笔刷半径
        const int animate_time = 1.0;
        double time_step = animate_time / (1.8 * stroke.nPoints() * stroke.nRadius());
        for (double t = 0; t <= animate_time; t += time_step)
        { 
            stroke.cubic_b_spline(nPoint, t);

            for(int hh = -R; hh <= R; hh++)
                for(int ww= -R; ww <= R; ww++){
                    bii = startPoint.y + hh;
                    bjj = startPoint.x + ww;
                    cbii = 50 + hh;
                    cbjj = 50 + ww;
                    if (bii < 0) bii = 0;
                    if (bii >= height) bii = height - 1;
                    if (bjj < 0) bjj = 0;
                    if (bjj >= width) bjj = width -1;

                    r = ref.at<uchar>(bii, bjj, 2)*PA + (uchar)stroke.getColorR()*(1-PA);
                    g = ref.at<uchar>(bii, bjj, 1)*PA + (uchar)stroke.getColorG()*(1-PA);
                    b = ref.at<uchar>(bii, bjj, 0)*PA + (uchar)stroke.getColorB()*(1-PA);

                    p_xx = newPoint.x + ww;
                    p_yy = newPoint.y + hh;
                    if(!(p_xx >= 0 && p_yy >=0 && p_xx < width && p_yy < height)) continue;

                    int ii = hh + R, jj = ww + R;
                    int I = brush_tex.at<uchar>(ii, jj);
                    int Index = (int)p_yy * step + (int)p_xx * 3;

                    // 附加纹理
                    pdst[Index+2] = r*(I/255.0) + pdst[Index+2]*(1- I/255.0);
                    pdst[Index+1] = g*(I/255.0) + pdst[Index+1]*(1- I/255.0);
                    pdst[Index+0] = b*(I/255.0) + pdst[Index+0]*(1- I/255.0);

                    // 计算高度图
                    int height_map_index = (int)p_yy*m_width + (int)p_xx;
                    count_pass[height_map_index]++;
                    sum_pass[height_map_index] += I;
                    height_maps.at<uchar>(p_yy, p_xx) = (int)((sum_pass[height_map_index])/(count_pass[height_map_index]));
                }

		}
    }

    clock_t finish = clock();
    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    cout<<"time of running CubicBSpline drawing steps: "<< duration <<endl;
}

// image rendering
void PainterlyService::render()
{
    assert(m_init);
    
    vector<SplineStroke> strokes_queue[6];
    clock_t begin, finish, start, end;
    Mat ref;
    
    begin = clock();
    for(int i = 0; i != m_nlayer; ++i)
    {
        if(m_brush_radius[i] == 0) continue;
        
        cout << "====== Painting Layer:" << i << ", stroke size: " << m_brush_radius[i] << " ======" << endl;
        start = clock();

        int ksize = (m_brush_radius[i] + 0.5) * 2 + 1; // kernel size        
        GaussianBlur(m_src, ref, Size(ksize, ksize), m_brush_radius[i]);

        generate_strokes(strokes_queue[i], ref, m_brush_radius[i]);
        paint_layer(ref, m_brush_radius[i], strokes_queue[i]);

        finish = clock();
        duration = (double)(finish - start) / CLOCKS_PER_SEC;
        cout<<"time of running paint_layer(): "<< duration <<endl;
    }

    end = clock();
    cout << "Total time" << (double)(end - begin) / CLOCKS_PER_SEC << "!!!!!\n" << endl;
}
