#include<QThread>
#include "load.h"
class loadThread : public QThread 
{ 
    Q_OBJECT 
public: 
    loadThread(QObject *parent = 0); 
    ~loadThread(); 
    virtual void run(); 
private: 
	load* progressBar;
};