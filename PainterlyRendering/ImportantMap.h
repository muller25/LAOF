#ifndef _ImportantMap_H
#define _ImportantMap_H

#include <cv.h>

class ImportantMap{
	public:
		ImportantMap(void);
		~ImportantMap(void);
		static double * important_energy;
		static void compute_important_map(double * important_energy,IplImage * src_image);
};

#endif
