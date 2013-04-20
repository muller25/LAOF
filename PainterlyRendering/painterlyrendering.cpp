#include "painterlyrendering.h"

PainterlyRendering::PainterlyRendering(QWidget *parent, Qt::WFlags flags)
	: QDialog(parent, flags)
{
	createItems();
	createLayouts();
	createConnections();
}

PainterlyRendering::~PainterlyRendering()
{	
	if( loadImage!=NULL ){
		delete loadImage;
		loadImage = NULL;
	}
	if( saveImage!=NULL ){
		delete saveImage;
		saveImage = NULL;
	}
	if(imageView!=NULL){
		delete imageView;
		imageView = NULL;
	}
	if( src!=NULL ){
		delete src;
		src = NULL;
	}
	if( res!=NULL ){
		delete res;
		res = NULL;
	}
}

void PainterlyRendering::createItems(){
	src = NULL;
	res = NULL;
	imageView = new PWidget;
	paintArea = new QScrollArea(this);
	paintArea->setWidget(imageView);

	max_stroke_length_spinbox = new QSpinBox;
    max_stroke_length_slider = new QSlider(Qt::Horizontal);
    max_stroke_length_spinbox->setRange(8,32);   //设置范围大小
    max_stroke_length_slider->setRange(8,32);
	max_stroke_length_spinbox->setValue(16);

	loadImage = new QPushButton("Load");
	renderImage = new QPushButton("Rendering");
	saveImage = new QPushButton("Save");

}

void PainterlyRendering::createLayouts(){

	QVBoxLayout* functionLayout = new QVBoxLayout;
	functionLayout->addStretch();
	//QHBoxLayout* sp = new QHBoxLayout;
	//sp->addStretch();
	//sp->addWidget(max_stroke_length_spinbox);
	//functionLayout->addWidget(max_stroke_length_slider);
	//functionLayout->addLayout(sp);

	functionLayout->addStretch();	
	functionLayout->addWidget(loadImage);
	functionLayout->addWidget(renderImage);
	functionLayout->addWidget(saveImage);
	functionLayout->addStretch();	

	QHBoxLayout* mainLayout = new QHBoxLayout;
	mainLayout->addWidget(paintArea);
	mainLayout->addLayout(functionLayout);
	setLayout(mainLayout);
}

void PainterlyRendering::createConnections(){

	connect(loadImage,SIGNAL(clicked()),this,SLOT(onPushLoad()));
	connect(saveImage,SIGNAL(clicked()),this,SLOT(onPushSave()));
	connect(renderImage,SIGNAL(clicked()),this,SLOT(onPushRender()));

	connect(max_stroke_length_spinbox,SIGNAL(valueChanged(int)),max_stroke_length_slider,SLOT(setValue(int)));
    connect(max_stroke_length_slider,SIGNAL(valueChanged(int)),max_stroke_length_spinbox,SLOT(setValue(int)));

}

void PainterlyRendering::resize(){

	int w = (src->width()) > 800 ? 800 : (src->width());
	int h = (src->height()) > 600 ? 600 : (src->height());
	this->setFixedSize(w + loadImage->width() + 32, h + 25);
	
}

void PainterlyRendering::onPushLoad()
{

	QString path = QFileDialog::getOpenFileName(0,"Load",".","All files (*.*)");
	if (path.isEmpty())
	{
		return;
	}
	if (src != NULL)
	{
		delete src;
	}
	src = new QImage(path);//通过路径获取图像；
	imageView->setImage(src->copy());
	resize();
}

/*void PainterlyRendering::onPushRender(){

	IplImage* p = RenderingImage::Processing(src);
	res = ImageOperate::IplImageToQImage(p);
	imageView->setImage(res->copy());
	resize();
}*/

void PainterlyRendering::onPushSave(){

	QString file =QFileDialog::getSaveFileName(0,tr("choose a filename to save under"), QDir::currentPath(),"Images(*.jpg)");
	res->save(file,"jpg");
}