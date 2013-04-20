#ifndef PAINTERLYRENDERING_H
#define PAINTERLYRENDERING_H

#include <QtGui/QDialog>
#include<QpushButton>
#include <QLayout>
#include<QScrollArea>
#include<QFileDialog>
#include<QLabel>
#include <QSlider>
#include <QSpinBox>
#include"PWidget.h"
#include "ImageOperate.h"
#include "RenderingImage.h"

class PainterlyRendering : public QDialog
{
	Q_OBJECT

public:
	PainterlyRendering(QWidget *parent = 0, Qt::WFlags flags = 0);
	~PainterlyRendering();
	void createItems();
	void createLayouts();
	void createConnections();
	void resize();
public slots:

	void onPushLoad();
	void onPushSave(); 
	void onPushRender();

private:

	QScrollArea* paintArea;
	QPushButton* loadImage;
	QPushButton* renderImage;
	QPushButton* saveImage;
	QLabel* label;
	QSpinBox * max_stroke_length_spinbox;
	QSlider *max_stroke_length_slider;
	PWidget* imageView;


	QImage* src;
	QImage* res;
};

#endif // PAINTERLYRENDERING_H
