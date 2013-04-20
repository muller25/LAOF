#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui/QMainWindow>
#include "ui_MainWindow.h"
#include "PWidget.h"
#include "EdgeWidget.h"
#include <QFileDialog>
#include <QLayout>
#include <QScrollArea>
#include <QRect>
#include<QMessageBox>
#include "ImageOperate.h"
#include "RenderingImage.h"
#include "Rendering.h"
#include "ColorAdjust.h"
#include "loadThread.h"
#include"load.h"
#include "QProgressIndicator.h"

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = 0, Qt::WFlags flags = 0);
	~MainWindow();
	void initData();
	void createConnects();
	void initEdge();
	int getMaxStrokeLenght();
	void setRSize();


public slots:
	void onMenuLoad();
    void onMenuRender();
	void onMenuSave();
	void onMenuRedo();
	void gainEdge();
	void onMenuEdge();
	void onMenuLight();
	void colorsAdjust();
	void onOK();
	void onMenuOptiRedo();

private:
	Ui::MainWindow ui;

	//QScrollArea* scrollArea;
	QImage* src;
	//QImage* tmp;
	QImage* res;
	IplImage* ss;
	PWidget* imageView;
	EdgeWidget* edgeView;
	IplImage *cedge, *gray, *sedge,*image;//±ﬂ‘µ–≈œ¢

	int m_angle;
    int m_timerId;
    int m_delay;
	QColor m_color;
	int op_id;
	QImage* re,*reLight,*reEdge;
	//QProgressIndicator* pi;
	//load* l;
};

#endif 