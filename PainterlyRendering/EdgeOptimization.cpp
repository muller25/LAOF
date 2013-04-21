#include"EdgeOptimization.h"
#include "ImportantMap.h"

EdgeOptimization::EdgeOptimization(void){}
EdgeOptimization::~EdgeOptimization(void){}

IplImage * EdgeOptimization::OptimizationEdge(IplImage * src,IplImage * ref,int para){
	ImportantMap::important_energy = new double[src->width*src->height];
	ImportantMap::compute_important_map(ImportantMap::important_energy,src);
	for(int i=0;i<src->height;i++)
		for(int j=0;j<src->width;j++){
			if(ImportantMap::important_energy[i*src->width+j]>para){
				CvScalar ss = cvGet2D(src,i,j);
				cvSet2D(ref,i,j,ss);
			}
		}
	return ref;
}
