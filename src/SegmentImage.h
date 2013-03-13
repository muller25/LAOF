#ifndef _SegmentImage_H
#define _SegmentImage_H

#include "Image.h"
#include "Maths.h"
#include "SegmentGraph.h"
#include <vector>
#include <map>

// dissimilarity measure between pixels
template<class T>
inline T diff(Image<T> &features, int x1, int y1, int x2, int y2)
{
    assert(features.ptr() != NULL);

    int width = features.nWidth(), channels = features.nChannels();
    int offset1, offset2;
    T dist = 0;
    T *pf = features.ptr();

    offset1 = (y1 * width + x1) * channels;
    offset2 = (y2 * width + x2) * channels;
    for (int i = 0; i < channels; ++i)
        dist += (pf[offset1] - pf[offset2]) * (pf[offset1] - pf[offset2]);
    
    return dist;
}

/*
 * Segment an image
 *
 * Returns a color image representing the segmentation.
 *
 * im: image to segment.
 * sigma: to smooth the image.
 * c: constant for treshold function.
 * min_size: minimum component size (enforced by post-processing stage).
 * num_ccs: number of connected components in the segmentation.
 */
template <class T, class T1 >
int segment_image(Image<T> &seg, const Image<T1> &im, double c, int min_size)
{
    const int r = 2;
    const double eps = 1e-2;
    
    int width = im.nWidth(), height = im.nHeight(), channels = im.nChannels();
    int num_ccs;
    std::vector< Image<T1> > ims, sIms;
    Image<T1> tmp, features;
    
    split(ims, im);
    for (int i = 0; i < channels; ++i)
    {
        GuidedFilter(tmp, ims[i], ims[i], r, eps);
        sIms.push_back(tmp);
    }

    mergec(features, sIms);
    
    // build graph
    edge *edges = new edge[width*height*4];
    int num = 0;
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            if (x < width-1)
            {
                edges[num].a = y * width + x;
                edges[num].b = y * width + (x+1);
                edges[num].w = diff(features, x, y, x+1, y);
                num++;
            }
            if (y < height-1)
            {
                edges[num].a = y * width + x;
                edges[num].b = (y+1) * width + x;
                edges[num].w = diff(features, x, y, x, y+1);
                num++;
            }
            if ((x < width-1) && (y < height-1))
            {
                edges[num].a = y * width + x;
                edges[num].b = (y+1) * width + (x+1);
                edges[num].w = diff(features, x, y, x+1, y+1);
                num++;
            }
            if ((x < width-1) && (y > 0))
            {
                edges[num].a = y * width + x;
                edges[num].b = (y-1) * width + (x+1);
                edges[num].w = diff(features, x, y, x+1, y-1);
                num++;
            }
        }
    }

    // segment
    universe *u = segment_graph(width*height, num, edges, c);
  
    // post process small components
    for (int i = 0; i < num; i++)
    {
        int a = u->find(edges[i].a);
        int b = u->find(edges[i].b);
        if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
            u->join(a, b);
    }
    delete [] edges;
    num_ccs = u->num_sets();

    int idx, label, count = 0;
    std::map<int, int> hash;
    std::map<int, int>::iterator iter;
    
    seg.create(width, height);
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            idx = y * width + x;
            label = u->find(idx);

            iter = hash.find(label);
            if (iter == hash.end())
            {
                hash.insert(std::pair<int, int>(label, count));
                label = count++;
            }
            else
                label = iter->second;
            
            seg[idx] = label;
        }
    }
    
    delete u;
    return num_ccs;
}

#endif
