#include "GCoptimization.h"
#include "Image.h"
#include "ImageIO.h"

#include <vector>

int width;

double smoothFn(int p1, int p2, int l1, int l2, void *pData)
{
    const double penalty = 0.3;
    double *ptr = (double *)pData;

    return fabs(l1-l2);
}

int main(int argc, char *argv[])
{
    const char *inIm = "/home/iaml/Projects/exp/lena/in/lena.avi.%03d.bmp";
    const char *inFlow = "/home/iaml/Projects/exp/lena/out/%s%03d.yml";
    char buf[256];
    DImage im, u, v, extra, result;
    int height, size, labels = 5, idx = 0;
    std::vector<DImage> vec;
    
    sprintf(buf, inIm, idx);
    imread(im, buf);
    sprintf(buf, inFlow, "u", idx);
    imreadf(u, buf);
    sprintf(buf, inFlow, "v", idx);
    imreadf(v, buf);

    vec.push_back(im);
    vec.push_back(u);
    vec.push_back(v);
    mergec(extra, vec);
    vec.clear();
    assert(extra.nChannels() == 5);
    
    size = im.nSize(), width = im.nWidth(), height = im.nHeight();
    double *data = new double [size * labels];

    memset(data, 0, sizeof(double)*size*labels);
    result.create(width, height);

    try{
		GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(width, height, labels);
		gc->setDataCost(data);
        gc->setSmoothCost(smoothFn, extra.ptr());
        
		printf("Before optimization energy is %.6f\n",gc->compute_energy());
		gc->swap(2);
		printf("After optimization energy is %.6f\n",gc->compute_energy());

		for ( int  i = 0; i < size; i++ )
			result[i] = gc->whatLabel(i);

		delete gc;
	}
	catch (GCException e){
		e.Report();
	}

    imshow("result", result);
    imwait(0);
    
    delete []data;
    
    return 0;
}
