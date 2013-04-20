#include "RenderService.h"
#include <fstream>
using namespace std;

#define PA 0.4
float RenderService::animate_time = 0.0;           // current animation frame
IplImage * RenderService::src_image;            // source image
IplImage * RenderService::dst_image;            // destination image
IplImage * RenderService::ref_image;            // reference image
int * RenderService::dif_image;                        // difference image

IplImage * RenderService::brush;
//IplImage * RenderService::brush_tex;//纹理笔刷
IplImage * RenderService::brush_tex1;//纹理笔刷
IplImage * RenderService::brush_tex2;//纹理笔刷
IplImage * RenderService::brush_tex3;//纹理笔刷
IplImage * RenderService::brush_tex4;//纹理笔刷
IplImage * RenderService::brush_tex5;//纹理笔刷
IplImage * RenderService::brush_tex6;//纹理笔刷

IplImage * RenderService::edge_image; //边缘图
IplImage * RenderService::height_maps;//新高度图
IplImage * RenderService::new_msk_maps;//新透明图

int * RenderService::count_pass;//通过像素点数次
int * RenderService::sum_pass;//通过像素点像素值之和
float * RenderService::depth_image;// depth image 

Queue * RenderService::strokes_queue;//笔划队列
map<string, PainterlyStyle> RenderService::styles;
PainterlyStyle RenderService::currentStyle;


int RenderService::descend(const void * a, const void * b){
    int * aa = (int *) a;
    int * bb = (int *) b;

    if (*aa < *bb)
        return 1;
    return 0;
}


/*****************************************************************************************
			笔划控制点计算
			x0,y0笔划起点
			R笔划大小
			ref_image参考图像
******************************************************************************************/
SplineStroke* RenderService::make_spline_stroke(int x0, int y0, int R, const IplImage * ref_image){

    int color_index = y0 * ref_image->widthStep + x0*3;     // array index //widthStep
    unsigned int r = (uchar)ref_image->imageData[color_index + 2];                   // red
    unsigned int g = (uchar)ref_image->imageData[color_index + 1];               // green
    unsigned int b = (uchar)ref_image->imageData[color_index + 0];               // blue

    int x = x0;             // current x coord
    int y = y0;             // current y coord

    float last_dx = 0;      // last x direction
    float last_dy = 0;      // last y direction

    double jr = 0.0;        // jittered red
    double jg = 0.0;        // jittered green       
    double jb = 0.0;        // jittered blue
    double jh = 0.0;        // jittered hue
    double js = 0.0;        // jittered saturation
    double jv = 0.0;        // jittered value

    int i, j, k;            // for loop counters

    SplineStroke * spline_stroke;

    int Gx[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}};

    int Gy[3][3] = {
        {-1, -2, -1},
        { 0,  0,  0},
        { 1,  2,  1}};

    ColorModels:: rgb_to_hsv((double) r / 256.0, (double) g / 256.0, (double) b / 256.0, &jh, &js, &jv);


    jh += currentStyle.jitter_hue * (rand() % 1000 - 500) / 1000.0 * 360.0 ;
    js += currentStyle.jitter_sat * (rand() % 1000 - 500) / 1000.0;
    jv += currentStyle.jitter_val * (rand() % 1000 - 500) / 1000.0;

    ColorModels::hsv_to_rgb(&jr, &jg, &jb, jh, js, jv);

    r = (int) (jr * 256);
    g = (int) (jg * 256);
    b = (int) (jb * 256);

    jr = currentStyle.jitter_r * (rand() % 1000 - 500) / 1000.0;
    jg = currentStyle.jitter_g * (rand() % 1000 - 500) / 1000.0;
    jb = currentStyle.jitter_b * (rand() % 1000 - 500) / 1000.0;

    r += (int) (jr * 256);
    g += (int) (jg * 256);
    b += (int) (jb * 256);

    spline_stroke = SplineStrokeService::spline_stroke_create(R, r, g, b,x0,y0);

    SplineStrokeService::spline_stroke_add(spline_stroke, x0, y0);

	 
    for (k = 1; k <= currentStyle.max_stroke_length; k++){

        float gx = 0;
        float gy = 0;
        float dx;
        float dy;
        float sum = 0;
        int dif_index = (y * ref_image->widthStep + x*3) ;//widthStep

        int d1 = (uchar)ref_image->imageData[dif_index] + (uchar)ref_image->imageData[dif_index + 1] + (uchar)ref_image->imageData[dif_index + 2] - (uchar)dst_image->imageData[dif_index] - (uchar)dst_image->imageData[dif_index + 1] -(uchar) dst_image->imageData[dif_index + 2];
        int d2 = (uchar)ref_image->imageData[dif_index] + (uchar)ref_image->imageData[dif_index + 1] + (uchar)ref_image->imageData[dif_index + 2] - r - g - b;

        if (k > currentStyle.min_stroke_length && (abs(d1) < abs(d2)))
            return spline_stroke;

        for(i = 0;i<3;i++){
            for(j = 0 ;j<3;j++){
                int ii = x - i - 1;
                int jj = y - j - 1;

                int sobel_index;

                if (ii < 0 || ii >= ref_image->width)
                    continue;
                if (jj < 0 || jj >= ref_image->height)
                    continue;

                sobel_index = jj * ref_image->widthStep + ii * 3;//widthStep

                gx += Gx[j][i] * (((uchar)ref_image->imageData[sobel_index] + (uchar)ref_image->imageData[sobel_index + 1] + (uchar)ref_image->imageData[sobel_index + 2] )/ 3);
                gy += Gy[j][i] * (((uchar)ref_image->imageData[sobel_index] + (uchar)ref_image->imageData[sobel_index + 1] + (uchar)ref_image->imageData[sobel_index + 2]) / 3);
            }
		}

        if (gx * gx + gy * gy <= 0.0001)
            return spline_stroke;

        sum = sqrt(gx * gx + gy * gy);
        spline_stroke->sumGradient = spline_stroke->sumGradient+sum;

        gx /= sum;
        gy /= sum;
        dx = -gy;
        dy = gx;

        if (last_dx * dx + last_dy * dy < 0){
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
        y = y + R * dy;

        last_dx = dx;
        last_dy =  dy;

        //判断是不是边缘
        int edge_index = y*edge_image->widthStep+x*3;
        unsigned int  rr= (uchar)edge_image->imageData[edge_index + 2];
        unsigned int  gg= (uchar)edge_image->imageData[edge_index + 1];
        unsigned int  bb= (uchar)edge_image->imageData[edge_index + 0];

        if (x < 0 || x >= ref_image->width || y < 0 || y >= ref_image->height||(rr>10&&gg>10&&bb>10))
            return spline_stroke;

        SplineStrokeService::spline_stroke_add(spline_stroke, x, y);
    }

    return spline_stroke;
}

/******************************************************************************
   计算两张图差异
   dif_image 差异结果
   ref_image 参考图
   dst_image 渲染结果图
   R 当前笔划大小
******************************************************************************/
 
void RenderService::difference_image(int * dif_image, const IplImage * ref_image, const IplImage * dst_image, int R){
    int i, j;
    int dif_index,img_index;
    int dif_temp;

#pragma omp parallel for private(j,dif_index,img_index,dif_temp)//这里并行修改
    for (i = 0; i < src_image->width; i++)
    {
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

}

ofstream fcc("f:\\out.txt");

void RenderService::paint_layer(IplImage * dst_image, IplImage * ref_image, int R){

    int grid_step = currentStyle.grid_size * R; 
    int i, j;
    int k = 0;

    difference_image(dif_image, ref_image, dst_image,R);//计算画板与参考图差异
	
#pragma omp parallel for private(i)	//这里并行修改
    for (j = src_image->height - 1; j >= 0; j -= grid_step)
    {
        for (i = 0; i < src_image->width; i += grid_step)
        {
            int area_error = 0.0;
            int current_error = -1;
            int ii, jj;
            int index_i = i, index_j = j;

            for (ii = -(grid_step / 2); ii <= grid_step / 2; ii++)
            {
                for (jj = -(grid_step / 2); jj <= grid_step / 2; jj++)
                {
                    int iii = i + ii;
                    int jjj = j + jj;
				
                    if (iii < 0)
                        iii = 0;

                    if (iii >= src_image->width)
                        iii = src_image->width - 1;

                    if (jjj < 0)
                        jjj = 0;

                    if (jjj >= src_image->height)
                        jjj = src_image->height - 1;

                    if (current_error <= dif_image[jjj * src_image->width + iii])
                    {
                        current_error = dif_image[jjj * src_image->width + iii];
                        index_i = iii;
                        index_j = jjj;
                    }

                    area_error += dif_image[jjj * src_image->width + iii];
                }
            }
            area_error /= (grid_step * grid_step);
			 
            if (area_error > currentStyle.threshold)
            {
                SplineStroke * stroke = make_spline_stroke(index_i, index_j, R, ref_image);

#pragma omp critical(strokes_queue)
                {
                    QueueService::queue_add(strokes_queue, stroke);
                }
            }
        }
    }


    i=strokes_queue->size;
    SplineStroke *stroke;

#pragma omp parallel for
    for(j=0;j<i;++j)
    {	 	

#pragma omp critical(strokes_queue)//这里并行修改
		{
            stroke= (SplineStroke *)QueueService::queue_remove(strokes_queue);
        }//比上面高效
        int Rs;
        float aver_value = (stroke->sumGradient/stroke->num_points);
        getStrokeSize(aver_value,R,&Rs);
        fcc<<Rs<<" ";
        IplImage* brush_tex;
		if(Rs == 32)
            brush_tex = brush_tex6;
        else if(Rs == 16)
            brush_tex = brush_tex5;
        else if(Rs == 8)
            brush_tex = brush_tex4;
        else if(Rs == 4)
            brush_tex = brush_tex3;
        else if(Rs == 2)
            brush_tex = brush_tex2;
        else
            brush_tex = brush_tex1;

        std::vector<Brush_color_2> brushcolor;
        Brush_color_2 bcolor;
        int point_x,point_y;
        int bi,bj,bii,bjj;
        point_x = stroke->start_point_x;
        point_y = stroke->start_point_y;
		 
        for(bi=-Rs;bi<=Rs;bi++){
            for(bj=-Rs;bj<=Rs;bj++){
                bii = point_y+bi;
                bjj = point_x+bj;
                if(bii<0) bii =0;
                if(bii>=ref_image->height) bii = ref_image->height-1;
                if(bjj<0) bjj =0;
                if(bjj>=ref_image->width) bjj = ref_image->width -1;
                CvScalar color1 = cvGet2D(ref_image,bii,bjj);
                bcolor.r = (uchar)color1.val[2]*PA + (uchar)stroke->r*(1-PA);
                bcolor.g = (uchar)color1.val[1]*PA + (uchar)stroke->g*(1-PA);
                bcolor.b = (uchar)color1.val[0]*PA + (uchar)stroke->b*(1-PA);
                brushcolor.push_back(bcolor);
            }
        }
        //
        CvScalar color = CV_RGB((double)stroke->r,(double)stroke->g,(double)stroke->b);
        CvPoint pt_now,pt_next;
       
        //绘制CubicBSpline曲线
        if(stroke->num_points> 0)
        {   double t,n;   // time
            double x,x0;   // x coord
            double y,y0;   // y coord
            uchar r,g,b;
            CubicBSpline CB;
            int dstIndex=0,texIndex=0,mapIndex=0;
            animate_time = 1.0;

            for (t = 0; t <= animate_time; t += animate_time / (1.8*stroke->num_points * stroke->s))//t的细分程度=系数*控制点个数*笔刷半径
            { 
                CB.cubic_b_spline(stroke, t, &x, &y);
                pt_now.x=x;
                pt_now.y=y;

                int brushR = Rs;
                int p_xx,p_yy;
                for(int hh = -brushR;hh<= brushR;hh++)
                    for(int ww= -brushR;ww<= brushR;ww++){
                        p_xx = pt_now.x+ww;
                        p_yy = pt_now.y+hh;
                        int ii = hh+brushR;
                        int jj = ww+brushR;
                        //int I = (uchar)*(PainterlyService::brush_tex->imageData+PainterlyService::brush_tex->widthStep*ii+jj*3);
                        //CvScalar cc = cvGet2D(RenderService::brush_tex,ii,jj);
                        CvScalar cc = cvGet2D(brush_tex,ii,jj);
                        int I = cc.val[0];
                        int index_tex = ii*brushR+jj;
                        if(p_xx>=0&&p_yy>=0&&p_xx<dst_image->width&&p_yy<dst_image->height){
                            int Index=(int)p_yy*dst_image->widthStep+(int)p_xx*3;
                            int index_tex = jj*brush_tex->width+ii;
                            dst_image->imageData[Index+2] = (uchar)brushcolor[index_tex].r*(I/255.0)+(uchar)dst_image->imageData[Index+2]*(1- I/255.0);
                            dst_image->imageData[Index+1] = (uchar)brushcolor[index_tex].g*(I/255.0)+(uchar)dst_image->imageData[Index+1]*(1- I/255.0);
                            dst_image->imageData[Index+0] = (uchar)brushcolor[index_tex].b*(I/255.0)+(uchar)dst_image->imageData[Index+0]*(1- I/255.0);
                            //计算高度图

                            int height_map_index = (int)p_yy*height_maps->width+(int)p_xx;
                            count_pass[height_map_index]++;
                            sum_pass[height_map_index] = sum_pass[height_map_index]+I;
                            int aver = (int)((sum_pass[height_map_index])/(count_pass[height_map_index]));
                            int msk_index = ii*new_msk_maps->widthStep+jj*3;
                            height_maps->imageData[height_maps->widthStep*p_yy+p_xx] = (int)((sum_pass[height_map_index])/(count_pass[height_map_index]));
								
                        }
                        else{
                            continue;
                        }
                    }

            }

		}

    }
			
}

void RenderService::loadStyles(string name, PainterlyStyle style)
{
    styles.insert(std::make_pair(name, style));
}

void RenderService::getStrokeSize(float aver_gradient,int R,int * Rs){
	if(aver_gradient<5.0)
		* Rs = R*2;
	else if(aver_gradient>=5.0&&aver_gradient<10.0)
		* Rs = R*2;
	else if(aver_gradient>=10.0&&aver_gradient<60.0)
		* Rs = R*1;
	else if(aver_gradient>=60.0&&aver_gradient<100.0)
		* Rs = R*0.5;
	else if(aver_gradient>=100.0&&aver_gradient<200.0)
		* Rs = R*0.5;
	else
		* Rs = 1;
	if(*Rs<1)
		*Rs =1;
	if(*Rs>32)
		*Rs = 32;

}
