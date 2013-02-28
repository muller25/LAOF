#ifndef _UTILS_H
#define _UTILS_H

#include <cv.h>
using namespace cv;

#include <iostream>
using namespace std;

template<class T>
bool array_match(T *exp, T *act, int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        if (fabs(exp[i] - act[i]) <= 10E-6) continue;

        cout << "expected[" << i << "(" << exp[i]
             << ") != actual[" << i << "(" << act[i] << ")\n";
        return false;
    }

    return true;
}

template<class T>
void setData(Mat &m, const T *data)
{
    int rows = m.rows;
    int cols = m.cols;
    int channels = m.channels();
    int r, c, k, offset;
    int step = m.step / sizeof(T);
    T *ptr = (T *)m.data;
    
    for (r = 0; r < rows; r++)
    {
        for (c = 0; c < cols; c++)
        {
            offset = r * step + c * channels;
            for (k = 0; k < channels; k++)
                ptr[offset+k] = data[offset+k];
        }
    }
}

template<class T>
bool matrix_match(const Mat &m1, const Mat &m2)
{
    assert(matchAll(m1, m2));

    int rows = m1.rows, cols = m1.cols, channels = m1.channels();
    int step = m1.step / sizeof(T);
    T *p1 = (T *)m1.data, *p2 = (T *)m2.data;
    int offset;
    
    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            offset = r * step + c * channels;
            for (int k = 0; k < channels; k++)
            {
                if (fabs(p1[offset+k] - p2[offset+k]) > ESP)
                {
                    std::cout << "expected[" << r << "][" << c << "][" << k << "](" << p1[offset+k]
                              << ") != actual[" << r << "][" << c << "][" << k << "](" << p2[offset+k] << ")\n";
                   
                    return false;
                }
            }
            
        }
    }

    return true;
}

#endif
