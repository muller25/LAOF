#ifndef PWIDGET_H
#define PWIDGET_H

#include<QWidget>

class PWidget : public QWidget 
{ 
public: 
        PWidget();
		~PWidget();
		void setImage(const QImage &image);
		QSize			sizeHint() const;
 
protected: 
        void paintEvent(QPaintEvent *event); 
private:

	QImage* displayImage;
};
#endif