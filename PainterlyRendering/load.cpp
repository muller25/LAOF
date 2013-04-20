#include "load.h"
#include<qpainter.h>

load::load(QWidget *parent, Qt::WFlags flags)
	: QWidget(parent,flags)
{
	ui.setupUi(this);
	m_angle = 0;
    m_timerId = -1;
    m_delay = 40;
    m_color = Qt::red;
	startAnimation();
	//connect(ui.pushButton,SIGNAL(clicked()),this,SLOT(stopAnimation()));
}

bool load::isAnimated () const
{
    return (m_timerId != -1);
}

void load::startAnimation()
{
    m_angle = 0;

    if (m_timerId == -1)
        m_timerId = startTimer(m_delay);
	
}

void load::stopAnimation()
{
    if (m_timerId != -1)
        killTimer(m_timerId);

    m_timerId = -1;

    update();
}

load::~load()
{

}
void load::timerEvent(QTimerEvent * /*event*/)
{
    m_angle = (m_angle+30)%360;

    update();
}

void load::paintEvent(QPaintEvent * /*event*/)
{
	if (!isAnimated())
        return;
    int width = qMin(this->width(), this->height());
    
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);
    
    int outerRadius = (width-1)*0.2;
    int innerRadius = (width-1)*0.2*0.38;

    int capsuleHeight = outerRadius - innerRadius;
    int capsuleWidth  = (width > 32 ) ? capsuleHeight *.23 : capsuleHeight *.35;
    int capsuleRadius = capsuleWidth/2;

    for (int i=0; i<12; i++)
    {
        QColor color = m_color;
        color.setAlphaF(1.0f - (i/12.0f));
        p.setPen(Qt::NoPen);
        p.setBrush(color);       
        p.save();
        p.translate(rect().center());
        p.rotate(m_angle - i*30.0f);
        p.drawRoundedRect(-capsuleWidth*0.5, -(innerRadius+capsuleHeight), capsuleWidth, capsuleHeight, capsuleRadius, capsuleRadius);
        p.restore();
    }
}
