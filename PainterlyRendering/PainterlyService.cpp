#include "PainterlyService.h"
#include "CubicBSpline.h"
#include"SplineStroke.h"
#include <cv.h>
#include<iostream>
#include"time.h"
#include<vector>

#define PI 3.1415926
#define MAX(a,b)     (a>b?a:b)
#define MAX(a,b,c)   MAX(MAX(a,b),c)
#define MIN(a,b)     (a<b?a:b)
#define MIN(a,b,c)   MIN(MIN(a,b),c)
#define PA 0.4
using namespace std;

int PainterlyService::NUM_STROKES_PER_FRAME = 500;     // number of strokes per frame
int PainterlyService::ANIMATE_STROKE = 1;             // animation flag
float PainterlyService::ANIMATE_SPEED = 0.1;          // animation speed
int PainterlyService::layer_done = 0;                 // current layer done
int PainterlyService::paint_done = 0;                 // painting done
float PainterlyService::loading_bar = 0.0;            // green loading bar flag
float PainterlyService::loading_bar_draw = 0.0;       // red loading bar flag
float PainterlyService::animate_time = 0.0;           // current animation frame
float PainterlyService::load_step = 0.0;              // loading bar progression step
int PainterlyService::is_locked_for_animation = 0;    // queue locked status
int PainterlyService::first_run = 1;                  // first time entering draw_scene
time_t PainterlyService::current_time = 0;            // random number generator seed
int PainterlyService::output_written = 0;             // output file written       


IplImage* PainterlyService::src_image;            // source image
IplImage*  PainterlyService::dst_image;            // destination image
IplImage* PainterlyService::ref_image;            // reference image
//IplImage* PainterlyService::buf_image;            // buffer image
IplImage* PainterlyService::tex_image;
IplImage* PainterlyService::msk_image;            //mask image;
IplImage* PainterlyService::map_image;
IplImage* PainterlyService::brush_tex;           //纹理笔刷
IplImage* PainterlyService::brush_test;          //测试笔刷形状；
IplImage* PainterlyService::edge_image;          //边缘图；

bool PainterlyService::useTexture;
int* PainterlyService::dif_image;                        // difference image

IplImage* PainterlyService::gradx_maps;                 //x方向梯度
IplImage* PainterlyService::grady_maps;                 //y方向梯度
double* PainterlyService::grad_orient;

IplImage* PainterlyService::height_maps;                 //新高度图；
IplImage* PainterlyService::new_msk_maps;                //新透明图；
int* PainterlyService::count_pass;                       //通过像素点数次
int* PainterlyService::sum_pass;                         //通过像素点像素值之和

IplImage* PainterlyService::color_map;

bool is_showControlLines = true;
int divs = 50; //控制精细度

float* PainterlyService::depth_image;                    // depth image

Queue* PainterlyService::strokes_queue;
map<string, PainterlyStyle> PainterlyService::styles;
PainterlyStyle PainterlyService::currentStyle;

int color_count = 0;
// ofstream fout("f://output.txt");

int PainterlyService::descend(const void * a, const void * b)
{
    int * aa = (int *) a;
    int * bb = (int *) b;

    if (*aa < *bb)
        return 1;
    return 0;
}

/******************************************************************************/
/*** Make Spline Stroke                                                     ***/
/***                                                                        ***/
/*** Input: x0, y0                                                          ***/
/***        brush size R                                                    ***/
/***        ref_image                                                       ***/
/*** Modifies: brushes_queue                                                ***/
/******************************************************************************/
SplineStroke* PainterlyService::make_spline_stroke(int x0, int y0, int R, const IplImage * ref_image)
{
    int color_index = y0 * ref_image->widthStep + x0*3;     // array index //widthStep
    unsigned int r = (uchar)ref_image->imageData[color_index + 2];                   // red
    unsigned int g = (uchar)ref_image->imageData[color_index + 1];               // green
    unsigned int b = (uchar)ref_image->imageData[color_index + 0];               // blue

    int x = x0;             // current x coord
    int y = y0;             // current y coord

//    float last_dx = 0;      // last x direction
//    float last_dy = 0;      // last y direction

    double jr = 0.0;        // jittered red
    double jg = 0.0;        // jittered green       
    double jb = 0.0;        // jittered blue
    double jh = 0.0;        // jittered hue
    double js = 0.0;        // jittered saturation
    double jv = 0.0;        // jittered value

    int k;            // for loop counters

    SplineStroke * spline_stroke;
/*
    int Gx[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}};

    int Gy[3][3] = {
        {-1, -2, -1},
        { 0,  0,  0},
        { 1,  2,  1}};
*/
    // convert rgb to hue  
    ColorModels:: rgb_to_hsv((double) r / 256.0, (double) g / 256.0, (double) b / 256.0, &jh, &js, &jv);

    // jitter hsv
    jh += currentStyle.jitter_hue * (rand() % 1000 - 500) / 1000.0 * 360.0 ;
    js += currentStyle.jitter_sat * (rand() % 1000 - 500) / 1000.0;
    jv += currentStyle.jitter_val * (rand() % 1000 - 500) / 1000.0;

    // convert hue back to rgb
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

    // create empty spline stroke
    spline_stroke = SplineStrokeService::spline_stroke_create(R, r, g, b,x0,y0);

    // add control point to spline stroke
    SplineStrokeService::spline_stroke_add(spline_stroke, x0, y0);

    // add control points to spline stroke
	 
    for (k = 1; k <= currentStyle.max_stroke_length; k++)
    {
//        float gx = 0;
//        float gy = 0;
//        float dx, dy;
//        float sum = 0;
        if (x < 0 || x >= ref_image->width || y < 0 || y >= ref_image->height)
            return spline_stroke;
        int dif_index = (y * ref_image->widthStep + x*3) ;//widthStep
        int d1 = (uchar)ref_image->imageData[dif_index] + (uchar)ref_image->imageData[dif_index + 1] + (uchar)ref_image->imageData[dif_index + 2] - (uchar)dst_image->imageData[dif_index] - (uchar)dst_image->imageData[dif_index + 1] -(uchar) dst_image->imageData[dif_index + 2];
        int d2 = (uchar)ref_image->imageData[dif_index] + (uchar)ref_image->imageData[dif_index + 1] + (uchar)ref_image->imageData[dif_index + 2] - r - g - b;

        if (k > currentStyle.min_stroke_length && (abs(d1) < abs(d2)))
            return spline_stroke;

		//for (i = -1; i <=1; i++)
        /*for(i = 0;i<3;i++)//以前是这样的
          {
          //for (j = -1; j <=1; j++)
          for(j = 0 ;j<3;j++)//以前是这样的
          {
          int ii = x - i - 1;//以前是这样
          int jj = y - j - 1;
          //int ii = x - i;
          //int jj = y - j;

          int sobel_index;

          if (ii < 0 || ii >= ref_image->width)
          continue;
          if (jj < 0 || jj >= ref_image->height)
          continue;

          sobel_index = jj * ref_image->widthStep + ii * 3;//widthStep

          gx += Gx[j][i] * (((uchar)ref_image->imageData[sobel_index] + (uchar)ref_image->imageData[sobel_index + 1] + (uchar)ref_image->imageData[sobel_index + 2] )/ 3);
          gy += Gy[j][i] * (((uchar)ref_image->imageData[sobel_index] + (uchar)ref_image->imageData[sobel_index + 1] + (uchar)ref_image->imageData[sobel_index + 2]) / 3);
          // gx += Gx[j][i]*((uchar)ref_image->imageData[sobel_index]*0.11+(uchar)ref_image->imageData[sobel_index+1]*0.59+(uchar)ref_image->imageData[sobel_index+2]*0.30);
          //gy += Gy[j][i]*((uchar)ref_image->imageData[sobel_index]*0.11+(uchar)ref_image->imageData[sobel_index+1]*0.59+(uchar)ref_image->imageData[sobel_index+2]*0.30);
          }

          }*/
        // gx = gradx_maps->imageData[y*ref_image->widthStep+x*3];
        //gy = grady_maps->imageData[y*ref_image->widthStep+x*3];

		/* sum = sqrt(gx * gx + gy * gy);
           spline_stroke->sumGradient = spline_stroke->sumGradient+sum;//计算笔划梯度总和

           if (gx * gx + gy * gy <= 0.0001)
           return spline_stroke;

           gx /= sum;
           gy /= sum;
           dx = -gy;
           dy = gx;

           if (last_dx * dx + last_dy * dy < 0)
           {
           dx = -dx;
           dy = -dy;
           }

           dx = currentStyle.curvature_filter * dx + (1.0 - currentStyle.curvature_filter) * last_dx;
           dy = currentStyle.curvature_filter * dy + (1.0 - currentStyle.curvature_filter) * last_dy;
           sum = sqrt(dx * dx + dy * dy);
		 
           if (sum < 0.001)
           return spline_stroke;

           dx = dx / sum;
           dy = dy / sum;
		 
           x = x + R * dx;
           y = y + R * dy;*/
        if (x < 0 || x >= ref_image->width || y < 0 || y >= ref_image->height)
			return spline_stroke;
		else{
            x = x+R*cos(grad_orient[y*src_image->width+x]+PI/2);
            y = y+R*sin(grad_orient[y*src_image->width+x]+PI/2);
		}

		/*last_dx = dx;
          last_dy =  dy;

          //判断是不是边缘
          int edge_index = y*edge_image->widthStep+x*3;
          unsigned int  rr= (uchar)edge_image->imageData[edge_index + 2];
          unsigned int  gg= (uchar)edge_image->imageData[edge_index + 1];
          unsigned int  bb= (uchar)edge_image->imageData[edge_index + 0];

          if (x < 0 || x >= ref_image->width || y < 0 || y >= ref_image->height||(rr>10&&gg>10&&bb>10))
          return spline_stroke;*/

        SplineStrokeService::spline_stroke_add(spline_stroke, x, y);
    }

    return spline_stroke;
}

/******************************************************************************/
/*** Difference Image                                                       ***/
/*** - Calculates the difference between src_image and dst_image and        ***/
/***   stores it to dif_image                                               ***/
/***                                                                        ***/
/*** Input: src_image                                                       ***/
/***        dst_image                                                       ***/
/***        dif_image                                                       ***/
/***                                                                        ***/
/*** Modifies: dif_image                                                    ***/
/******************************************************************************/
void PainterlyService::difference_image(int * dif_image, const IplImage * ref_image, 
                                        const IplImage * dst_image, int R)
{
    int i, j;
    int dif_index,img_index;
    int dif_temp;

    //计算单个点差异
// #pragma omp parallel for private(j,dif_index,img_index,dif_temp)//这里并行修改
    for (i = 0; i < src_image->width; i++)
    {
        //#pragma omp parallel for
        for (j = 0; j < src_image->height; j++)
        {
            dif_index = j * ref_image->width + i;
            img_index =j * ref_image->widthStep + i*3;//widthStep

            dif_temp = (uchar)dst_image->imageData[img_index] - (uchar)ref_image->imageData[img_index];
            dif_image[dif_index] =  dif_temp * dif_temp;
            dif_temp = (uchar)dst_image->imageData[img_index + 1] - (uchar)ref_image->imageData[img_index + 1];
            dif_image[dif_index] += dif_temp * dif_temp;
            dif_temp = (uchar)dst_image->imageData[img_index + 2] - (uchar)ref_image->imageData[img_index + 2];
            dif_image[dif_index] += dif_temp * dif_temp;
            dif_image[dif_index] = sqrt((float)dif_image[dif_index]);
        }

    }

    //差异滤波

}

/******************************************************************************/
/*** Paint Layer                                                            ***/	
/*** - Find areas that we need to paint and draw strokes in them            ***/
/***                                                                        ***/
/*** Input: ref_image                                                       ***/
/***        dst_image                                                       ***/
/***        Brush Size R                                                    ***/
/*** Modifies: dst_image                                                    ***/
/******************************************************************************/
void PainterlyService::paint_layer(IplImage * dst_image, IplImage * ref_image, int R)
{
    cout << "image paint_layer()" << endl;
    
    //atexit(Exit);
    int grid_step = currentStyle.grid_size * R; 
    int i, j;
    clock_t   start,   finish;   
    double     duration;

    load_step = (float) src_image->width  / ((src_image->height / grid_step) * (src_image->width / grid_step));
    loading_bar = 0;

    start = clock();
    cout << "before difference_image()" << endl;
    
    difference_image(dif_image, ref_image, dst_image,R);
    //difference_image(dif_image, ref_image, src_image, R);//parellel
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    cout<<"time of running difference_image(): "<< duration <<endl;

    start = clock();
// #pragma omp parallel for private(i)	//这里并行修改
    for (j = src_image->height - 1; j >= 0; j -= grid_step)
    {
        for (i = 0; i < src_image->width; i += grid_step)
        {
            int area_error = 0.0;
            int current_error = -1;
            int ii, jj;
            int index_i = i, index_j = j;

            //#pragma omp parallel for private(jj) reduction(+:area_error)
            for (ii = -(grid_step / 2); ii <= grid_step / 2; ii++)
            {
                for (jj = -(grid_step / 2); jj <= grid_step / 2; jj++)
                {
                    int iii = i + ii;
                    int jjj = j + jj;
                    //边界控制
                    //排除图像像素点小于0区域					
                    if (iii < 0)
                        iii = 0;
                    //排除图像像素点大于sizeX区域
                    if (iii >= src_image->width)
                        iii = src_image->width - 1;
                    //排除图像像素点小于0区域
                    if (jjj < 0)
                        jjj = 0;
                    //排除图像像素点大于sizeX区域
                    if (jjj >= src_image->height)
                        jjj = src_image->height - 1;
                    //计算(x,y)区域最大差值点及差值
                    //#pragma omp critical(current_error)
                    if (current_error <= dif_image[jjj * src_image->width + iii])
                    {
                        current_error = dif_image[jjj * src_image->width + iii];
                        index_i = iii;
                        index_j = jjj;
                    }

                    area_error += dif_image[jjj * src_image->width + iii];
                }
            }

            // 当区域平均差值大于阈值时，绘制笔刷
            area_error /= (grid_step * grid_step);
            if (area_error > currentStyle.threshold)
            {
                SplineStroke * stroke = make_spline_stroke(index_i, index_j, R, ref_image);

// #pragma omp critical(strokes_queue)//这里并行修改
                {
                    QueueService::queue_add(strokes_queue, stroke);
                }
                loading_bar += load_step;
            }
        }
    }
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    cout<<"time of running stroke construction steps: "<< duration <<endl;
    start = clock();

// CubicBSpline版本
//spline_strokes->a queue of control points->rgb, x,y 
    //得到笔刷的control points
    cout<<"strokes size:"<<strokes_queue->size<<endl;
    i=strokes_queue->size;
    SplineStroke *stroke;
    //#pragma omp parallel for
    //for(i =strokes_queue->size ;i >0 ; ++i)
    //while(strokes_queue->size>0)

// #pragma omp parallel for
    for(j=0;j<i;++j)
    {	 	
// #pragma omp critical(strokes_queue)//这里并行修改
		{
            stroke= (SplineStroke *)QueueService::queue_remove(strokes_queue);
        }//比上面高效
		// fout<<stroke->sumGradient/stroke->num_points<<" ";

        /*
          制定brush中每bristle确定颜色
        */
        std::vector<Brush_color> brushcolor;
        Brush_color bcolor;
        int point_x,point_y;
        int bi,bj,bii,bjj;
        int cbii,cbjj;
        point_x = stroke->start_point_x;
        point_y = stroke->start_point_y;
        //#pragma omp parallel for private(bj,bcolor,bii,bjj,point_x,point_y)	//这里并行修改
        //#pragma omp parallel for
        for(bi=-R;bi<=R;bi++){
            for(bj=-R;bj<=R;bj++){
                bii = point_y+bi;
                bjj = point_x+bj;
                cbii = 50+bi;
                cbjj = 50+bj;
                if(bii<0) bii =0;
                if(bii>=ref_image->height) bii = ref_image->height-1;
                if(bjj<0) bjj =0;
                if(bjj>=ref_image->width) bjj = ref_image->width -1;
                CvScalar color1 = cvGet2D(ref_image,bii,bjj);
                if(color_count == 0){
                    cvSet2D(color_map,cbii,cbjj,color1);
                }
                bcolor.r = (uchar)color1.val[2]*PA + (uchar)stroke->r*(1-PA);
                bcolor.g = (uchar)color1.val[1]*PA + (uchar)stroke->g*(1-PA);
                bcolor.b = (uchar)color1.val[0]*PA + (uchar)stroke->b*(1-PA);
                brushcolor.push_back(bcolor);
            }
        }
        if(color_count==0){
            cvSaveImage("color.jpg",color_map);
            color_count++;
        }
		 
        // CvScalar color = CV_RGB((double)stroke->r,(double)stroke->g,(double)stroke->b);
        CvPoint pt_now;
       
        //绘制CubicBSpline曲线
        if(stroke->num_points> 0)
        {
            double t;   // time
            double x, y;   // x-y coord
            // uchar r,g,b;
            CubicBSpline CB;
            // int dstIndex=0,texIndex=0,mapIndex=0;
            animate_time = 1.0;

            for (t = 0; t <= animate_time; t += animate_time / (1.8*stroke->num_points * stroke->s))//t的细分程度=系数*控制点个数*笔刷半径
            { 
                CB.cubic_b_spline(stroke, t, &x, &y);
                pt_now.x=x;
                pt_now.y=y;

                /*if(useTexture)
                  {
                  double x_next,y_next,t_next;
                  t_next=t + animate_time / (1.8*stroke->num_points * stroke->s);

                  if(t_next>1) 
                  t_next=t - animate_time / (1.8*stroke->num_points * stroke->s);

						
                  CB.cubic_b_spline(stroke, t_next, &x_next, &y_next);
                  double theta=atan2(y_next-y,x_next-x);


                  if(theta<0) theta+=2*PI;
                  if(theta>2*PI) theta-=2*PI;
                  double xo=x+R*cos(theta-PI/2);
                  double yo=y+R*sin(theta-PI/2);

                  for(n=0.0;n<=1.0;n+=1.0/(1.3*2*R))
                  {
                  x0=xo+n*2*R*cos(theta+PI/2);
                  y0=yo+n*2*R*sin(theta+PI/2);
                  if(x0<0||y0<0||x0>=dst_image->width||y0>=dst_image->height)
                  {
                  continue;
                  }

                  dstIndex=(int)y0*map_image->widthStep+(int)x0*3;
                  texIndex=(int)(tex_image->height*n)*tex_image->widthStep+(int)(tex_image->width*t)*3;
                  //mapIndex=(int)(tex_image->height*n)*tex_image->widthStep+(int)(tex_image->width*t);

                  uchar c=tex_image->imageData[texIndex+1];
								
                  //纯灰度高度笔刷
                  map_image->imageData[dstIndex+0]=(uchar)map_image->imageData[dstIndex+0]*((255.0-(uchar)msk_image->imageData[texIndex])/255.0)+c*(((uchar)msk_image->imageData[texIndex])/255.0);
                  map_image->imageData[dstIndex+1]=(uchar)map_image->imageData[dstIndex+1]*((255.0-(uchar)msk_image->imageData[texIndex])/255.0)+c*(((uchar)msk_image->imageData[texIndex])/255.0);
                  map_image->imageData[dstIndex+2]=(uchar)map_image->imageData[dstIndex+2]*((255.0-(uchar)msk_image->imageData[texIndex])/255.0)+c*(((uchar)msk_image->imageData[texIndex])/255.0);

                  }

                  }*/
                //else
                //{
                //cvCircle(dst_image, pt_now, R, color, -1);
                //} 
                //cvCircle(dst_image, pt_now, R, color, -1);
                int brushR = R;
                int p_xx,p_yy;
                // 计算笔刷的纹理和高度图
                for(int hh = -brushR;hh<= brushR;hh++)
                    for(int ww= -brushR;ww<= brushR;ww++){
                        p_xx = pt_now.x+ww;
                        p_yy = pt_now.y+hh;
                        int ii = hh+brushR;
                        int jj = ww+brushR;
                        //int I = (uchar)*(PainterlyService::brush_tex->imageData+PainterlyService::brush_tex->widthStep*ii+jj*3);
                        CvScalar cc = cvGet2D(PainterlyService::brush_tex,ii,jj);
                        int I = cc.val[0];
                        int index_tex = ii*brushR+jj;
                        if(p_xx>=0&&p_yy>=0&&p_xx<dst_image->width&&p_yy<dst_image->height){
                            int Index=(int)p_yy*dst_image->widthStep+(int)p_xx*3;
                            // int index_tex = jj*brush_tex->width+ii;
                            //dst_image->imageData[Index+2] = (uchar)color.val[2]*(I/255.0)+(uchar)dst_image->imageData[Index+2]*(1- I/255.0);
                            //dst_image->imageData[Index+1] = (uchar)color.val[1]*(I/255.0)+(uchar)dst_image->imageData[Index+1]*(1- I/255.0);
                            //dst_image->imageData[Index+0] = (uchar)color.val[0]*(I/255.0)+(uchar)dst_image->imageData[Index+0]*(1- I/255.0);
                            // 附加纹理
                            dst_image->imageData[Index+2] = (uchar)brushcolor[index_tex].r*(I/255.0)+(uchar)dst_image->imageData[Index+2]*(1- I/255.0);
                            dst_image->imageData[Index+1] = (uchar)brushcolor[index_tex].g*(I/255.0)+(uchar)dst_image->imageData[Index+1]*(1- I/255.0);
                            dst_image->imageData[Index+0] = (uchar)brushcolor[index_tex].b*(I/255.0)+(uchar)dst_image->imageData[Index+0]*(1- I/255.0);

                            //计算高度图
                            int height_map_index = (int)p_yy*height_maps->width+(int)p_xx;
                            count_pass[height_map_index]++;
                            sum_pass[height_map_index] = sum_pass[height_map_index]+I;
                            // int aver = (int)((sum_pass[height_map_index])/(count_pass[height_map_index]));
                            // int msk_index = ii*new_msk_maps->widthStep+jj*3;
                            height_maps->imageData[height_maps->widthStep*p_yy+p_xx] = (int)((sum_pass[height_map_index])/(count_pass[height_map_index]));
                            //height_maps->imageData[height_maps->widthStep*p_yy+p_xx]=(uchar)height_maps->imageData[height_maps->widthStep*p_yy+p_xx]*((255.0-(uchar)new_msk_maps->imageData[msk_index])/255.0)+(uchar)(aver*(((uchar)new_msk_maps->imageData[msk_index])/255.0));
                            //if(j ==2||j==10||j==40||j==60||j==100||j==120||j==200||j==400||j==800||j==1000){
                            if(j<0)	{
                                brush_test->imageData[Index+2] = (uchar)brushcolor[index_tex].r*(I/255.0)+(uchar)brush_test->imageData[Index+2]*(1- I/255.0);
                                brush_test->imageData[Index+1] = (uchar)brushcolor[index_tex].g*(I/255.0)+(uchar)brush_test->imageData[Index+1]*(1- I/255.0);
                                brush_test->imageData[Index+0] = (uchar)brushcolor[index_tex].b*(I/255.0)+(uchar)brush_test->imageData[Index+0]*(1- I/255.0);
                            }
                        }
                    }
            }

//整线位置
#ifdef SHOWMODE
            cvShowImage("show",map_image);
            cvWaitKey(1);
#endif
		}
    }

    // fout.close();

    //for(i = 0; i< dst_image->width; ++i)
    //	for(j = 0;j< dst_image->height; ++j)
    // 	{
    //		int color_index2 = (j * dst_image->width + i) * 3;
    //             dst_image->imageData[color_index2 + 2] = ((img->imageData + j*img->widthStep))[i*img->nChannels + 2]; // B
    //             dst_image->imageData[color_index2 + 1] = ((img->imageData + j*img->widthStep))[i*img->nChannels + 1]; // G
    //	   dst_image->imageData[color_index2] = ((img->imageData + j*img->widthStep))[i*img->nChannels + 0]; // R
    //	}
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    cout<<"time of running CubicBSpline drawing steps: "<< duration <<endl;

    //  cvReleaseImage(&img);
}

void PainterlyService::loadStyles(string name, PainterlyStyle style)
{
    styles.insert(std::make_pair(name, style));
}

