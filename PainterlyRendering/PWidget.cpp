#include "PWidget.h"
#include"qpainter.h"
#include "cv.h"
#include<QTextStream>
#include<iostream>
#include<loadThread.h>
using namespace std;

PWidget::PWidget() 
{		
		setWindowTitle(tr("Paint Demo")); 
		displayImage	= new QImage(600, 700, QImage::Format_RGB32);
		displayImage->fill(qRgb(230, 230, 230));

}
PWidget::~PWidget(){

	if(displayImage!= NULL){
		delete displayImage;
		displayImage = NULL;
	}
}
void PWidget::setImage(const QImage &image)
{
	*displayImage = image.convertToFormat(QImage::Format_RGB32);
	adjustSize();
	update();
	updateGeometry();
}

QSize PWidget::sizeHint() const
{
	return displayImage->size();
}

void PWidget::paintEvent(QPaintEvent *event) 
{		
	   //loadThread* load = new loadThread();
	   //load->start();
        QPainter painter(this); 
		//painter.drawImage(QPoint(0, 0), *displayImage);

		if(displayImage->width()<559 && displayImage->height()<519){
			
			int x = 559/2 - (displayImage->width())/2;
			int y = 519/2 - (displayImage->height())/2;
			painter.drawPixmap(QPoint(x,y),QPixmap::fromImage(*displayImage));
		}
		else{
			painter.drawImage(QPoint(0,0),*displayImage);
		}
		
}