/*
  Copyright (C) 2006 Pedro Felzenszwalb

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#include <cstdio>
#include "SegmentImage.h"
#include "Flow2Color.h"
#include "ImageIO.h"

int main(int argc, char **argv)
{
    if (argc != 5)
    {
        fprintf(stderr, "usage: %s k min input output\n", argv[0]);
        return 1;
    }
  
    double k = atof(argv[1]);
    int min_size = atoi(argv[2]);
	
    printf("loading input image.\n");
    DImage im, seg;
    UCImage out;
    imread(im, argv[3]);
	
    printf("processing\n");
    int num_ccs = segment_image(seg, im, k, min_size); 

    printf("seg %.6f .. %.6f\n", seg.min(), seg.max());

    substract(seg, (seg.max()-seg.min())/2);
    flow2color(out, seg, seg);
    imwrite(argv[4], out);

    printf("got %d components\n", num_ccs);
    printf("done! uff...thats hard work.\n");

    return 0;
}
