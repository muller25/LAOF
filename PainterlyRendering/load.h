#ifndef LOAD_H
#define LOAD_H

#include <QtGui/QWidget>
#include "ui_load.h"

class load : public QWidget
{
	Q_OBJECT

public:
	load(QWidget *parent = 0, Qt::WFlags flags = 0);
	~load();
	 bool isAnimated () const;
public slots:
	void startAnimation();
	void stopAnimation();

protected:
    virtual void timerEvent(QTimerEvent * event); 
    virtual void paintEvent(QPaintEvent * event);

private:
	Ui::Form ui;

	int m_angle;
    int m_timerId;
    int m_delay;
	bool m_displayedWhenStopped;
    QColor m_color;
};

#endif // QPROBAR_H