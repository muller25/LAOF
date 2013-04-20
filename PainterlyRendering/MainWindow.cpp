#include "MainWindow.h"
#include"qpainter.h"
#include<QDebug>

MainWindow::MainWindow(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	ui.setupUi(this);

	/*进度条数据*/
	m_angle = 0;
    m_timerId = -1;
    m_delay = 40;
    m_color = Qt::red;

	//pi = new QProgressIndicator();
	//l = new load();
	op_id = 0;
	initData();
	createConnects();

}

MainWindow::~MainWindow()
{

	if(src!=NULL){
		delete src;
		src = NULL;
	}
	/*if(tmp!=NULL){
		delete tmp;
		tmp = NULL;
	}*/
	if(res!=NULL){
		delete res;
		res = NULL;
	}
	if(imageView!=NULL){
		delete imageView;
		imageView = NULL;
	}
	if(edgeView!=NULL){
		delete edgeView;
		edgeView = NULL;
	}
	cvReleaseImage(&ss);
//	cvReleaseImage(&cedge);
//	cvReleaseImage(&gray);
//  cvReleaseImage(&sedge);
//	cvReleaseImage(&image);

}

void MainWindow::initData(){

	src = NULL;
	res = NULL; 
	//tmp = NULL;
	ss = 0 ;
	imageView = new PWidget;
	edgeView = new EdgeWidget;

	cedge = NULL, gray = NULL, sedge = NULL,image = NULL;

	/*图片显示滚动条有效办法*/
	//ui.scrollArea->styleSheet("QScrollArea{}");
	QSizePolicy policy = imageView->sizePolicy();
	policy.setHorizontalStretch(1);
	policy.setVerticalStretch(1);
	policy.setHorizontalPolicy(QSizePolicy::Minimum);
	policy.setVerticalPolicy(QSizePolicy::Minimum);
	imageView->setSizePolicy(policy);
	ui.scrollArea->setWidget(imageView);
	ui.scrollArea->setStyleSheet("background-image: url(../icon/pps.png);");

	/*边缘图显示滚动条有效办法*/
	QSizePolicy edgePolicy = edgeView->sizePolicy();
	edgePolicy.setHorizontalStretch(1);
	edgePolicy.setVerticalStretch(1);
	edgePolicy.setHorizontalPolicy(QSizePolicy::Minimum);
	edgePolicy.setVerticalPolicy(QSizePolicy::Minimum);
	edgeView->setSizePolicy(edgePolicy);
	ui.edgeScrollArea->setWidget(edgeView);

	
	ui.maxSpinBox->setRange(8,32);   //设置范围大小
	ui.maxSlider->setRange(8,32);
	ui.maxSpinBox->setValue(16);
	ui.maxSlider->setValue(16);

	ui.minSpinBox->setRange(1,4);
	ui.minSlider->setRange(1,4);
	ui.minSpinBox->setValue(2);
	ui.minSlider->setValue(2);

	ui.gridSpinBox->setRange(1,4);
	ui.gridSlider->setRange(1,4);
	ui.gridSpinBox->setValue(2);
	ui.gridSlider->setValue(2);

	ui.thresholdSpinBox->setRange(20,100);
	ui.thresholdSlider->setRange(20,100);
	ui.thresholdSpinBox->setValue(46);
	ui.thresholdSlider->setValue(46);

	ui.repairEdgeSpinBox->setRange(50,250);
	ui.repairEdgeSlider->setRange(50,250);
	ui.repairEdgeSpinBox->setValue(100);
	ui.repairEdgeSlider->setValue(100);

	ui.lightSpinBox->setRange(10,40);
	ui.lightHorizontalSlider_2->setRange(10,40);
	ui.lightSpinBox->setValue(25);
	ui.lightHorizontalSlider_2->setValue(25);

	ui.edgeSpinBox->setRange(0,300);
	ui.edgeHorizontalSlider->setRange(0,300);
	ui.edgeSpinBox->setValue(100);
	ui.edgeHorizontalSlider->setValue(100);

	/*图片特效处理部分*/
	ui.sHslider->setRange(0,512);
	ui.hHslider->setRange(0,512);
	ui.bHslider->setRange(0,512);
	ui.cHslider->setRange(0,50);

	ui.sHslider->setValue(234);
	ui.hHslider->setValue(250);
	ui.bHslider->setValue(230);
	ui.cHslider->setValue(12);

}

void MainWindow::createConnects(){

	connect(ui.maxSpinBox,SIGNAL(valueChanged(int)),ui.maxSlider,SLOT(setValue(int)));
	connect(ui.maxSlider,SIGNAL(valueChanged(int)),ui.maxSpinBox,SLOT(setValue(int)));

	connect(ui.minSpinBox,SIGNAL(valueChanged(int)),ui.minSlider,SLOT(setValue(int)));
	connect(ui.minSlider,SIGNAL(valueChanged(int)),ui.minSpinBox,SLOT(setValue(int)));

	connect(ui.gridSpinBox,SIGNAL(valueChanged(int)),ui.gridSlider,SLOT(setValue(int)));
	connect(ui.gridSlider,SIGNAL(valueChanged(int)),ui.gridSpinBox,SLOT(setValue(int)));

	connect(ui.thresholdSpinBox,SIGNAL(valueChanged(int)),ui.thresholdSlider,SLOT(setValue(int)));
	connect(ui.thresholdSlider,SIGNAL(valueChanged(int)),ui.thresholdSpinBox,SLOT(setValue(int)));

	connect(ui.repairEdgeSpinBox,SIGNAL(valueChanged(int)),ui.repairEdgeSlider,SLOT(setValue(int)));
	connect(ui.repairEdgeSlider,SIGNAL(valueChanged(int)),ui.repairEdgeSpinBox,SLOT(setValue(int)));

	connect(ui.lightSpinBox,SIGNAL(valueChanged(int)),ui.lightHorizontalSlider_2,SLOT(setValue(int)));
	connect(ui.lightHorizontalSlider_2,SIGNAL(valueChanged(int)),ui.lightSpinBox,SLOT(setValue(int)));

	connect(ui.edgeSpinBox,SIGNAL(valueChanged(int)),ui.edgeHorizontalSlider,SLOT(setValue(int)));
	connect(ui.edgeHorizontalSlider,SIGNAL(valueChanged(int)),ui.edgeSpinBox,SLOT(setValue(int)));
	connect(ui.edgeHorizontalSlider,SIGNAL(valueChanged(int)),this,SLOT(gainEdge()));

	connect(ui.actionOpen,SIGNAL(triggered()),this,SLOT(onMenuLoad()));
	connect(ui.actionRendering,SIGNAL(triggered()),this,SLOT(onMenuRender()));
	//connect(ui.actionRendering,SIGNAL(triggered()),pi,SLOT(startAnimation()));

	connect(ui.actionSave,SIGNAL(triggered()),this,SLOT(onMenuSave()));
	connect(ui.actionRedo,SIGNAL(triggered()),this,SLOT(onMenuRedo()));
	connect(ui.actionEdge,SIGNAL(triggered()),this,SLOT(onMenuEdge()));
	connect(ui.actionDefines,SIGNAL(triggered()),this,SLOT(onOK()));
	connect(ui.actionPrevious,SIGNAL(triggered()),this,SLOT(onMenuOptiRedo()));

	connect(ui.actionOptimizedLight,SIGNAL(triggered()),this,SLOT(onMenuLight()));
	/*图像特效处理部分*/
	connect(ui.sHslider,SIGNAL(valueChanged(int)),this,SLOT(colorsAdjust()));
	connect(ui.hHslider,SIGNAL(valueChanged(int)),this,SLOT(colorsAdjust()));
	connect(ui.bHslider,SIGNAL(valueChanged(int)),this,SLOT(colorsAdjust()));
	connect(ui.cHslider,SIGNAL(valueChanged(int)),this,SLOT(colorsAdjust()));

	ui.loadImage->addAction(ui.actionOpen);
	ui.saveImage->addAction(ui.actionSave);
	ui.RenderImage->addAction(ui.actionRendering);
	ui.RenderImage->addAction(ui.actionRedo);
	ui.gainEdge->addAction(ui.actionEdge);
	ui.gainEdge->addAction(ui.actionOptimizedLight);
	ui.gainEdge->addAction(ui.actionPrevious);
	ui.menu_2->addAction(ui.actionDefines);

}

int MainWindow::getMaxStrokeLenght(){
	return ui.maxSlider->value();
}
void MainWindow::onOK(){
	QImage* e = ImageOperate::IplImageToQImage(ss);
	src = e;
}
void MainWindow::colorsAdjust(){
	 if(src == NULL)
		QMessageBox::warning(this,"加载缺失","请先加载图片",QMessageBox::Ok,QMessageBox::Cancel);
	 else{
		 ss = ImageOperate::QImageToIplImage(*src);
		 ss = ColorAdjust::hsv_ad(ss,(int)ui.sHslider->value(),(int)ui.hHslider->value(),(int)ui.bHslider->value(),(float)((int)ui.cHslider->value()/10.0));
		 QImage* e = ImageOperate::IplImageToQImage(ss);
		 imageView->setImage(e->copy());
	 }
}


void MainWindow::onMenuLoad()
{
	//res = NULL;
	QString path = QFileDialog::getOpenFileName(0,"Load",".","All files (*.*)");
	if (path.isEmpty())
	{
		return;
	}
	if (src != NULL)
	{
		delete src;
	}
	src = new QImage(path);//通过路径获取图像;

	op_id = 0;//这个参数标记是哪方面优化

	//tmp = src;
	imageView->setImage(src->copy());
	initEdge();
	gainEdge();
}
void MainWindow::onMenuRender(){

	RenderingImage::max_stroke_length = (int)ui.maxSpinBox->value();
	RenderingImage::min_stroke_length = (int)ui.minSpinBox->value();
	RenderingImage::threshold = (int)ui.thresholdSpinBox->value();
	RenderingImage::grid_size = (float)ui.gridSpinBox->value();
	/*Rendering::max_stroke_length = (int)ui.maxSpinBox->value();
	Rendering::min_stroke_length = (int)ui.minSpinBox->value();
	Rendering::threshold = (int)ui.thresholdSpinBox->value();
	Rendering::grid_size = (float)ui.gridSpinBox->value();*/
	setRSize();
	//tmp = res;
	//loadThread* load = new loadThread();
	//load->run();
	//ui.progress->startAnimation();
	//ui.progress->show();

	double startTime,endTime;
	startTime = clock();
	IplImage* p = RenderingImage::Processing(src,cedge);
	//IplImage* p = Rendering::Processing(src,res,cedge);
	res = ImageOperate::IplImageToQImage(p);
	imageView->setImage(res->copy());
	endTime = clock();
	ui.textEdit->setText(QString::number((double)(endTime - startTime) / CLOCKS_PER_SEC)+"s");

}

void MainWindow::onMenuRedo(){

	if(src ==NULL){
		QMessageBox::warning(this,"返回原图","没有加载图片",QMessageBox::Ok,QMessageBox::Cancel);
	}
	else{
		//res = tmp;
		//if(res == NULL)
			imageView->setImage(src->copy());
			op_id = 0;
		//else
			//imageView->setImage(res->copy());
	}

}
void MainWindow::onMenuSave(){
	if(res ==NULL){
		QMessageBox::warning(this,"保存图像","当前没有图像可以保存",QMessageBox::Ok,QMessageBox::Cancel);
	}
	else{
		QString file =QFileDialog::getSaveFileName(0,tr("choose a filename to save under"), QDir::currentPath(),"Images(*.jpg)");
		res->save(file,"jpg");
	}

}


void MainWindow::initEdge(){
	
	image = ImageOperate::QImageToIplImage(*src);
	cedge = cvCreateImage(cvSize(image->width,image->height), IPL_DEPTH_8U, 3);

    // 将彩色图像转换为灰度图像
    gray = cvCreateImage(cvSize(image->width,image->height), IPL_DEPTH_8U, 1);
    sedge = cvCreateImage(cvSize(image->width,image->height), IPL_DEPTH_8U, 1);
    cvCvtColor(image, gray, CV_BGR2GRAY);
}

void MainWindow::gainEdge(){
	if(gray == NULL){
		QMessageBox::warning(this,"图像边缘","请先加载图像",QMessageBox::Ok,QMessageBox::Cancel);
	}
	else{
		cvSmooth( gray, sedge, CV_BLUR, 3, 3, 0 );
		cvNot( gray, sedge );

		// 对灰度图像进行边缘检测
		cvCanny(gray, sedge, (float)ui.edgeHorizontalSlider->value(), (float)ui.edgeHorizontalSlider->value()*3, 3);
		cvZero( cedge );
		// copy edge points
		cvCopy( image, cedge, sedge );
		QImage* e = ImageOperate::IplImageToQImage(cedge);
		edgeView->setImage(*e);
	}
}

void MainWindow::setRSize(){
	if(ui.checkBox_32->isChecked())
		RenderingImage::R[0] = 32;
	else
		RenderingImage::R[0] = 0;
	if(ui.checkBox_16->isChecked())
		RenderingImage::R[1] = 16;
	else
		RenderingImage::R[1] = 0;
	if(ui.checkBox_8->isChecked())
		RenderingImage::R[2] = 8;
	else
		RenderingImage::R[2] = 0;
	if(ui.checkBox_4->isChecked())
		RenderingImage::R[3] = 4;
	else
		RenderingImage::R[3] = 0;
	if(ui.checkBox_2->isChecked())
		RenderingImage::R[4] = 2;
	else
		RenderingImage::R[4] = 0;
	if(ui.checkBox_1->isChecked())
		RenderingImage::R[5] = 1;
	else
		RenderingImage::R[5] = 0;
}

void MainWindow::onMenuEdge(){

	int value = (int)ui.repairEdgeSpinBox->value();
	if(res == NULL){
		QMessageBox::warning(this,"边缘修复","没有渲染图片",QMessageBox::Ok,QMessageBox::Cancel);
	}
	else{
		if(op_id ==0){
			op_id = 1;
			re = res;
		}
		else{
			op_id = 3;
			reLight = res;
		}
		IplImage* p = RenderingImage::opereateEdge(res,value);
		//IplImage* p = Rendering::opereateEdge(res);
		res = ImageOperate::IplImageToQImage(p);
		imageView->setImage(res->copy());

	}

}

void MainWindow::onMenuLight(){
	int value = (int)ui.lightSpinBox->value();
	float para = (float)value/(100.00);
	if(res == NULL){
		QMessageBox::warning(this,"光照处理","没有渲染图片",QMessageBox::Ok,QMessageBox::Cancel);
	}
	else{
		if(op_id ==0){
			op_id = 2;
			re = res;
		}
		else{
			op_id = 4;
			reEdge = res;
		}
		IplImage* p = RenderingImage::operateLight(res,para);
		//IplImage* p = Rendering::operateLight(res);
		res = ImageOperate::IplImageToQImage(p);
		imageView->setImage(res->copy());
	}

}

void MainWindow::onMenuOptiRedo(){
	if(op_id == 0){

		QMessageBox::warning(this,"返回上一步","还没有进行优化处理",QMessageBox::Ok,QMessageBox::Cancel);
	}
	else if(op_id == 1||op_id == 2){
		res = re;
	}
	else if(op_id == 3){
		res = reLight;
		op_id = op_id -2;
	}
	else{
		res = reEdge;
		op_id = op_id -2;
	}
	imageView->setImage(res->copy());
}

