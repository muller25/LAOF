#include "RenderingImage.h"
#include <iostream>
#include <ctime> 

#include "GaussianBlur.h"
#include "PainterlyService.h"
#include "NormalMap.h"
#include "Lighting.h"
#include "MakeBrushMap.h"
#include "Lighting_Effect.h"
#include "Edge_Optimization.h"

using namespace std;

int RenderingImage::max_stroke_length = 0;
int RenderingImage::min_stroke_length = 0;
int RenderingImage::threshold = 0;
float RenderingImage::grid_size = 0.0;
int RenderingImage::R[nlayer];
char RenderingImage::brushMapPath[256];

RenderingImage::RenderingImage(){}
RenderingImage::~RenderingImage(){}

// 需要预先设置 stroke_length, threshold, grid_size, R, brushMapPath
IplImage* RenderingImage::Processing(IplImage* srcImage, IplImage* edgeImage){
    assert(srcImage != NULL);

    int srcWidth = srcImage->width, srcHeight = srcImage->height;
    int srcDepth = srcImage->depth, srcChannels = srcImage->nChannels;
    CvSize srcSize = cvGetSize(srcImage);
    clock_t start, finish, begin, end;
    double duration;

	PainterlyService::src_image = srcImage;
    PainterlyService::dst_image = cvCreateImage(srcSize,8,3);
    PainterlyService::ref_image = cvCreateImage(srcSize, 8, 3);
    PainterlyService::brush_test = cvCreateImage(srcSize, srcDepth, srcChannels);
    PainterlyService::edge_image = edgeImage;
		
    // Impressionism style
    PainterlyService::currentStyle.num_layers = nlayer;
    PainterlyService::currentStyle.brush_radius = R;
    PainterlyService::currentStyle.blur_factor = 0.5;
    PainterlyService::currentStyle.max_stroke_length = max_stroke_length;
	PainterlyService::currentStyle.min_stroke_length = min_stroke_length;
	PainterlyService::currentStyle.curvature_filter = 1.0;
	PainterlyService::currentStyle.threshold = threshold;
	PainterlyService::currentStyle.alpha = 1.0;
    PainterlyService::currentStyle.grid_size = grid_size;
	PainterlyService::currentStyle.jitter_r = 0.0;
	PainterlyService::currentStyle.jitter_g = 0.0;
	PainterlyService::currentStyle.jitter_b = 0.0;
	PainterlyService::currentStyle.jitter_hue = 0.0;
	PainterlyService::currentStyle.jitter_sat = 0.0;
	PainterlyService::currentStyle.jitter_val = 0.0;

    PainterlyService::strokes_queue = QueueService::queue_create();

    int src_size = srcWidth * srcHeight;
    PainterlyService::dif_image = new int[src_size];
    PainterlyService::count_pass = new int[src_size];
    PainterlyService::sum_pass = new int[src_size];
    PainterlyService::grad_orient = new double[src_size];
    
    memset(PainterlyService::count_pass,0,src_size*4);
    memset(PainterlyService::sum_pass,0,src_size*4);

    getStrokeOrientation(PainterlyService::src_image,PainterlyService::grad_orient);
    IplImage *brush = cvLoadImage(brushMapPath);

    PainterlyService::useTexture=true;
    PainterlyService::color_map = cvCreateImage(cvSize(100,100),8,3);
    /*if(PainterlyService::useTexture)//载入纹理
      {
      PainterlyService::tex_image=cvLoadImage("F:\\纹理\\brush.bmp");
      PainterlyService::msk_image=cvLoadImage("F:\\纹理\\mask.bmp");
      if(PainterlyService::tex_image==NULL||PainterlyService::msk_image==NULL)
      {
      cout<<"纹理载入失败！"<<endl;
      return 0;
      }
      }*/
    PainterlyService::map_image=cvCreateImage(srcSize,8,3);
    PainterlyService::height_maps = cvCreateImage(srcSize,8,1);
    CvScalar msk_color = CV_RGB(255.0,255.0,255.0);
    CvPoint msk_p;

    begin=clock();
    for(int i = 0; i != nlayer; ++i)
    {
        if(R[i] != 0){
            cout << "====== Painting Layer:" << i << ", stroke size: " << R[i] << " ======" << endl;
            start = clock();
            //PainterlyService::brush_tex = MakeBrushMap::NewBrushMap(R[i]);
            msk_p.x = R[i], msk_p.y = R[i];
            PainterlyService::brush_tex = cvCreateImage(cvSize(2*R[i]+1,2*R[i]+1),brush->depth,brush->nChannels);
            PainterlyService::new_msk_maps = cvCreateImage(cvSize(2*R[i]+1,2*R[i]+1),brush->depth,brush->nChannels);
            cvCircle(PainterlyService::new_msk_maps, msk_p, R[i], msk_color, -1);
            cvResize(brush,PainterlyService::brush_tex,CV_INTER_CUBIC);//设置透明图大小

            // int ksize = (R[i] + 0.5) * 2 + 1; // kernel size
            // cvSmooth(PainterlyService::src_image, PainterlyService::ref_image, ksize, ksize, R[i]);
            GaussianBlur::gaussian_blur(PainterlyService::ref_image, PainterlyService::src_image, R[i]);
            // cvSaveImage("F://实验图//guass.jpg",PainterlyService::ref_image);
            PainterlyService::paint_layer(PainterlyService::dst_image, PainterlyService::ref_image, R[i]);
            finish = clock();
            duration = (double)(finish - start) / CLOCKS_PER_SEC;
            cout<<"time of running paint_layer(): "<< duration <<endl;
        }
    }
    //cvSaveImage("F://实验图//col.jpg",PainterlyService::brush_test);

    if(PainterlyService::useTexture)
    {
        cout << "using texture" << endl;
        IplImage* graymap=cvCreateImage(srcSize,8,1);
        IplImage* ehgraymap=cvCreateImage(srcSize,8,1);

        cvCvtColor(PainterlyService::map_image,graymap,CV_RGB2GRAY);
        cvEqualizeHist(graymap,ehgraymap);
        IplImage* nmap=NormalMap::computeNormalMap(ehgraymap);
        IplImage* res=Lighting::phong(PainterlyService::dst_image,nmap);
        cvReleaseImage(&(PainterlyService::dst_image));
        PainterlyService::dst_image = res;
    }

    end=clock();
    cout << "Total time" <<(double)(end - begin) / CLOCKS_PER_SEC<<"!!!!!\n"<< endl;

    memset(PainterlyService::count_pass,0,sizeof(PainterlyService::count_pass));
    memset(PainterlyService::sum_pass,0,sizeof(PainterlyService::sum_pass));
    //delete PainterlyService::strokes_queue;
    //delete PainterlyService::dif_image;
    //delete PainterlyService::count_pass;
    //delete PainterlyService::sum_pass;

    return PainterlyService::dst_image;
}

IplImage* RenderingImage::operateLight(IplImage *dst_map, float para){
    assert(dst_map != NULL);
    
	IplImage* light_map = cvCreateImage(cvSize(PainterlyService::src_image->width,PainterlyService::src_image->height),8,3);
	IplImage* ehmap=cvCreateImage(cvSize(PainterlyService::src_image->width,PainterlyService::src_image->height),8,1);
	cvEqualizeHist(PainterlyService::height_maps,ehmap);

	light_map = Lighting_Effect::Add_Light(PainterlyService::src_image,ehmap,dst_map,para);
	return light_map;
}

IplImage* RenderingImage::opereateEdge(IplImage *dst_map, int para){
    assert(dst_map != NULL);
    
	IplImage* edge_map = cvCreateImage(cvSize(PainterlyService::src_image->width,PainterlyService::src_image->height),8,3);
	edge_map = Edge_Optimization::Optimization_edge(PainterlyService::src_image,dst_map,para);
	return edge_map;
}

void RenderingImage::getStrokeOrientation(IplImage* src_image, double* orientation){
    assert(src_image != NULL && orientation != NULL);
    
	//std::vector<array<double, 2> > positions;
    //std::vector<array<double, 1> > values;

	IplImage* gray,*sobel_x,*sobel_y;
	gray = cvCreateImage(cvGetSize(src_image),src_image->depth,1);
	sobel_x = cvCreateImage(cvGetSize(src_image),IPL_DEPTH_16S,1);
	sobel_y = cvCreateImage(cvGetSize(src_image),IPL_DEPTH_16S,1);
	cvCvtColor(src_image,gray,CV_BGR2GRAY);

	cvSobel(gray,sobel_x,1,0,3);//x方向求导
	cvSobel(gray,sobel_y,0,1,3);//y方向求导
		
	IplImage* sobel8_x=cvCreateImage(cvGetSize(sobel_x),IPL_DEPTH_8U,1);
	IplImage* sobel8_y=cvCreateImage(cvGetSize(sobel_y),IPL_DEPTH_8U,1);

	cvConvertScaleAbs(sobel_x,sobel8_x,1,0);
	cvConvertScaleAbs(sobel_y,sobel8_y,1,0);

	/*for(int i=0;i<sobel8_y->height;i+=10){
      for(int j=0;j<sobel8_y->width;j+=10){
      double max_grad = -1;int index_i;int index_j;
      for(int ii=i;ii<i+10&&ii<sobel8_y->height;ii++){
      for(int jj=j;jj<j+10&&jj<sobel8_y->width;jj++){
      int index = ii*sobel8_y->widthStep+jj;
      double gy = (double)sobel8_y->imageData[index];
      double gx = (double)sobel8_x->imageData[index];
      if(sqrt(gy*gy+gx*gx)>max_grad){
      index_i = ii;
      index_j = jj;
      max_grad = sqrt(gy*gy+gx*gx);
      }
      }
      }
      if(max_grad>150){
      array<double, 2> curPos;
      curPos[0] = index_i; curPos[1] = index_j;
      positions.push_back(curPos);
      array<double,1> curVal;
      int ind = index_i*sobel8_x->widthStep+index_j;
      double gy = sobel8_y->imageData[ind];
      double gx = sobel8_x->imageData[ind];
      curVal[0] = atan2(gy,gx);
      values.push_back(curVal);
      }
      }
			
      }*/
    // ThinPlateSpline<2,1> spline(positions, values);
    for(int i=0;i<src_image->height;i++){
        for(int j=0;j<src_image->width;j++){
            int index = i*sobel8_y->widthStep+j;
            double gy = (double)sobel8_y->imageData[index];
            double gx = (double)sobel8_x->imageData[index];
            //if(sqrt(gy*gy+gx*gx)<40){
            // array<double, 2> curPos;
            // curPos[0] = i;
            // curPos[1] = j;
            // PainterlyService::grad_orient[i*src_image->width+j] = spline.interpolate(curPos)[0];
            //  }
            // else{
            PainterlyService::grad_orient[i*src_image->width+j] = atan2(gy,gx);
            // }
        }
    }
}
