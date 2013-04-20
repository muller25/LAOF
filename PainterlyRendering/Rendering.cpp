#include "Rendering.h"


int Rendering::max_stroke_length = 0;
int Rendering::min_stroke_length = 0;
int Rendering::threshold = 0;
float Rendering::grid_size = 0.0;
int Rendering::R[6];


Rendering::Rendering(){


}
Rendering::~Rendering(){

}

IplImage* Rendering::Processing(const QImage* src,QImage* res,IplImage* edgeImage){
 
	if(src == NULL){
		QMessageBox::warning(0,"加载缺失","请先加载图片",QMessageBox::Ok,QMessageBox::Cancel);
		IplImage* dst = cvCreateImage(cvSize(100,100),8,3);
		return dst;
	}
	RenderService::src_image = ImageOperate::QImageToIplImage(*src);

	int size = 6;

	RenderService::dst_image = cvCreateImage(cvSize(RenderService::src_image->width, RenderService::src_image->height),8,3);
	if(res !=NULL)
		RenderService::dst_image = ImageOperate::QImageToIplImage(*res);
	RenderService::ref_image = cvCreateImage(cvSize(RenderService::src_image->width, RenderService::src_image->height),8,3);
	RenderService::edge_image =  edgeImage;
		
		//Impressionism style
	RenderService::currentStyle.initParameter(size,R,0.5,max_stroke_length,min_stroke_length,1.0,threshold,1.0,grid_size,0.0,0.0,0.0,0.0,0.0,0.0); 
	RenderService::strokes_queue = QueueService::queue_create();
	RenderService::dif_image = new int[RenderService::src_image->width*RenderService::src_image->height];
	RenderService::count_pass = new int[RenderService::src_image->width*RenderService::src_image->height];
	RenderService::sum_pass = new int[RenderService::src_image->width*RenderService::src_image->height];

	int src_size = RenderService::src_image->width*RenderService::src_image->height;
	memset(RenderService::count_pass,0,src_size*4);
	memset(RenderService::sum_pass,0,src_size*4);


	IplImage * brush = cvLoadImage("F:\\纹理\\brush_map.bmp");//加载笔刷纹理
	RenderService::height_maps = cvCreateImage(cvSize(RenderService::src_image->width,RenderService::src_image->height),8,1);
	CvScalar msk_color = CV_RGB(255.0,255.0,255.0);
    CvPoint msk_p;

	R[0] = 16;R[1] = 8;R[2]= 4;R[3] = 4;R[4] =2;R[5] = 2;

    
	RenderService::brush_tex1 = cvCreateImage(cvSize(2*1+1,2*1+1),brush->depth,brush->nChannels);
	RenderService::brush_tex2 = cvCreateImage(cvSize(2*2+1,2*2+1),brush->depth,brush->nChannels);
	RenderService::brush_tex3 = cvCreateImage(cvSize(2*4+1,2*4+1),brush->depth,brush->nChannels);
	RenderService::brush_tex4 = cvCreateImage(cvSize(2*8+1,2*8+1),brush->depth,brush->nChannels);
	RenderService::brush_tex5 = cvCreateImage(cvSize(2*16+1,2*16+1),brush->depth,brush->nChannels);
	RenderService::brush_tex6 = cvCreateImage(cvSize(2*32+1,2*32+1),brush->depth,brush->nChannels);

	cvResize(brush,RenderService::brush_tex1,CV_INTER_CUBIC);//设置透明图大小
	cvResize(brush,RenderService::brush_tex2,CV_INTER_CUBIC);//设置透明图大小
    cvResize(brush,RenderService::brush_tex3,CV_INTER_CUBIC);//设置透明图大小
	cvResize(brush,RenderService::brush_tex4,CV_INTER_CUBIC);//设置透明图大小
	cvResize(brush,RenderService::brush_tex5,CV_INTER_CUBIC);//设置透明图大小
	cvResize(brush,RenderService::brush_tex6,CV_INTER_CUBIC);//设置透明图大小

	
    for( int i = 0; i !=size; ++i){
		//if(R[i]!=0){
	
			//PainterlyService::brush_tex = cvCreateImage(cvSize(2*R[i]+1,2*R[i]+1),brush->depth,brush->nChannels);

			//cvResize(brush,PainterlyService::brush_tex,CV_INTER_CUBIC);//设置透明图大小

			GaussianBlur::gaussian_blur(RenderService::ref_image, RenderService::src_image, R[i]);
			RenderService::paint_layer(RenderService::dst_image, RenderService::ref_image, R[i]);
				
	  // }
    }

	memset(RenderService::count_pass,0,sizeof(RenderService::count_pass));
	memset(RenderService::sum_pass,0,sizeof(RenderService::sum_pass));

	 return RenderService::dst_image;
}

IplImage* Rendering::operateLight( QImage* res ){
	
	IplImage* dst_map = ImageOperate::QImageToIplImage(*res);
	IplImage* light_map = cvCreateImage(cvSize(RenderService::src_image->width,RenderService::src_image->height),8,3);
	IplImage* ehmap=cvCreateImage(cvSize(RenderService::src_image->width,RenderService::src_image->height),8,1);
	cvEqualizeHist(RenderService::height_maps,ehmap);

	//light_map = Lighting_Effect::Add_Light(RenderService::src_image,ehmap,dst_map);
	return light_map;

}

IplImage* Rendering::opereateEdge( QImage* res ){

	IplImage* dst_map = ImageOperate::QImageToIplImage(*res);
	IplImage* edge_map = cvCreateImage(cvSize(RenderService::src_image->width,RenderService::src_image->height),8,3);
	//edge_map = Edge_Optimization::Optimization_edge(RenderService::src_image,dst_map);
	return edge_map;


}