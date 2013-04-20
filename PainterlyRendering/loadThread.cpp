#include "loadThread.h"
loadThread::loadThread(QObject *parent) : QThread(parent) 
{ 

	progressBar = new load();

}
loadThread::~loadThread(){
	
	delete progressBar;
	progressBar = NULL;

}

void loadThread::run(){
    
	//progressBar->setWindowFlags(Qt::FramelessWindowHint);
	//progressBar->setAttribute(Qt::WA_TranslucentBackground, true);
	progressBar->show();

}