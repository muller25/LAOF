#ifndef _Image_H
#define _Image_H

#include <cv.h>

#ifndef ESP
#define ESP 1e-6
#endif

#define R 2
#define G 1
#define B 0

template <class T>
class Image
{
public:
    Image()
    {
        pData = NULL;
        width = height = channels = elements = 0;
    }

    Image(int _width, int _height, int _channels=1, T val=0);

    template <typename T1>
    Image(const Image<T1> &im):pData(NULL){copyFrom(im);};

    Image(const Image<T> &im):pData(NULL){copyFrom(im);}
    
    virtual ~Image(){release();}

    inline void create(int _width, int _height, int _channels=1, T val=0);
    
    inline void release()
    {
        if (pData != NULL)
            delete []pData;
        
        pData = NULL;
        width = height = channels = elements = 0;
    }

    void convertTo(cv::Mat &im) const;
    void convertFrom(const cv::Mat &im);
    
    template <typename T1>
    inline void set(T1 val)
    {
        for (int i = 0; i < elements; ++i)
            pData[i] = (T)val;
    }

    template <typename T1>
    inline void copyTo(Image<T1> &im) const;

    template <typename T1>
    inline void copyFrom(const Image<T1> &im);
    
    void threshold(T minVal, T maxVal);
    
    template <typename T1>
    inline bool match2D(const Image<T1> &im) const;
 
    template <typename T1>
    inline bool match3D(const Image<T1> &im) const;

    inline bool isFloat() const
    {
        return (typeid(T) == typeid(double) ||
                typeid(T) == typeid(float) ||
                typeid(T) == typeid(long double));
    }

    inline T& operator[](int idx){return pData[idx];}
    inline T& operator[](int idx) const{return pData[idx];}

    template <typename T1>
    inline Image<T>& operator=(const Image<T1> &m);
    inline Image<T>& operator=(const Image<T> &m);

    template <typename T1>
    void invTo(Image<T1> &m) const;
    
    inline bool isEmpty() const {return (pData==NULL || width==0 || height==0 || channels==0);}
    inline int nonZeros() const;
    inline int nSize() const {return width*height;}
    inline int nWidth() const{return width;}
    inline int nHeight() const{return height;}
    inline int nChannels() const{return channels;}
    inline int nElements() const{return elements;}
    inline T* ptr(){return pData;}
    inline T* ptr() const {return pData;}
    inline T min() const;
    inline T max() const;
    inline void normalize();

private:
    T * pData;
    int width, height, channels, elements;
};

template <class T>
Image<T>::Image(int _width, int _height, int _channels, T val)
{
    pData = NULL;
    create(_width, _height, _channels, val);
}

template <class T>
void Image<T>::convertTo(cv::Mat &im) const
{
    if (typeid(T) == typeid(uchar))
        im.create(height, width, CV_8UC(channels));
    else if (typeid(T) == typeid(int))
        im.create(height, width ,CV_32SC(channels));
    else if (typeid(T) == typeid(short))
        im.create(height, width ,CV_16SC(channels));
    else if (typeid(T) == typeid(unsigned short))
        im.create(height, width ,CV_16UC(channels));
    else if (typeid(T) == typeid(float))
        im.create(height, width, CV_32FC(channels));
    else if (typeid(T) == typeid(double))
        im.create(height, width, CV_64FC(channels));
    else {
        printf("unknown data type, can't convert Image<T> to cv::Mat\n");
        im.create(height, width, CV_64FC(channels));
    }
    
    int step = im.step1(), ioffset, offset;
    T *pi = (T *)im.data;
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = (h * width + w) * channels;
            ioffset = h * step + w * channels;
            for (int k = 0; k < channels; ++k)
                pi[ioffset + k] = pData[offset + k];
        }
    }
}

template <class T>
void Image<T>::convertFrom(const cv::Mat &im)
{
    int height = im.rows, width = im.cols, channels = im.channels();
    create(width, height, channels);

    int step = im.step1(), ioffset, offset;
    T *pi = (T *)im.data;
    for (int h = 0; h < height; ++h)
    {
        for (int w = 0; w < width; ++w)
        {
            offset = (h * width + w) * channels;
            ioffset = h * step + w * channels;
            for (int k = 0; k < channels; ++k)
                pData[offset + k] = pi[ioffset + k];
        }
    }
}

template <class T>
template <typename T1>
void Image<T>::copyTo(Image<T1> &im) const
{
    im.create(width, height, channels);

    for (int i = 0; i < elements; i++)
        im[i] = pData[i];
}

template <class T>
template <typename T1>
void Image<T>::copyFrom(const Image<T1> &im)
{
    create(im.nWidth(), im.nHeight(), im.nChannels());

    for (int i = 0; i < elements; ++i)
        pData[i] = im[i];
}

template <class T>
void Image<T>::create(int _width, int _height, int _channels, T val)
{
    if (width != _width || height != _height || channels != _channels ||
        pData == NULL)
    {
        release();
        width = _width;
        height = _height;
        channels = _channels;
        elements = width * height * channels;
        pData = new T[elements];
        assert(pData != NULL);
    }

    set(val);
}

template <class T>
void Image<T>::threshold(T minVal, T maxVal)
{
    for (int i = 0; i < elements; ++i)
    {
        if (pData[i] < minVal) pData[i] = minVal;
        if (pData[i] > maxVal) pData[i] = maxVal;
    }
}

template <class T>
int Image<T>::nonZeros() const
{
    int count = 0;
    for (int i = 0; i < elements; ++i)
        if (fabs(pData[i]) < ESP) ++count;

    return count;
}

template <class T>
T Image<T>::min() const
{
    T val = pData[0];
    for (int i = 1; i < elements; ++i)
        if (pData[i] < val) val = pData[i];

    return val;
}

template <class T>
T Image<T>::max() const
{
    T val = pData[0];
    for (int i = 1; i < elements; ++i)
        if (pData[i] > val) val = pData[i];

    return val;
}

template <class T>
void Image<T>::normalize()
{
    T maxm = max(), minm = min();
    for (int i = 0; i < elements; ++i)
        pData[i] = (pData[i] - minm) / (maxm - minm);
}

template <class T>
template <typename T1>
bool Image<T>::match2D(const Image<T1> &im) const
{
    return (width==im.nWidth() && height==im.nHeight());
}

template <class T>
template <typename T1>
bool Image<T>::match3D(const Image<T1> &im) const
{
    return (match2D(im) && channels == im.nChannels());
}

template <class T>
template <typename T1>
Image<T>& Image<T>::operator=(const Image<T1> &m)
{
    copyFrom(m);
    return *this;
}

template <class T>
Image<T>& Image<T>::operator=(const Image<T> &m)
{
    copyFrom(m);
    return *this;
}

template <class T>
template <typename T1>
void Image<T>::invTo(Image<T1> &m) const
{
    cv::Mat mat;
    convertTo(mat);

    // invert only accept float point matrix
    m.convertFrom(mat.inv(cv::DECOMP_SVD));
}

typedef Image<double> DImage;
typedef Image<uchar> UCImage;
typedef Image<int> IImage;

#endif
