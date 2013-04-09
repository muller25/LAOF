#include <opencv2/imgproc/imgproc.hpp>
#include <highgui.h>
#include <cv.h>

void LabHist(cv::Mat imat, cv::Mat maskMat=cv::Mat());

int main(int argc, char *argv[])
{
    cv::Mat src = cv::imread(argv[1]);
    LabHist(src);
    
    return 0;
}

// normalized lab histogram
void LabHist(cv::Mat imat, cv::Mat maskMat)
{
    cv::Mat lab, hist;
    
    // convert to lab color space
    // cv::Mat imat, lab, maskMat, tmp;

    // if (typeid(T) != typeid(float))
    // {
    //     im.covertTo(tmp);
    //     tmp.convertTo(imat, CV_32FC3);
    // }
    // else
    //     im.convertTo(imat);
        
    // cv::cvtColor(imat, lab, CV_BGR2Lab);

    // // convert mask
    // if (typeid(M) != typeid(uchar))
    // {
    //     mask.convertTo(tmp);
    //     tmp.convertTo(maskMat, CV_8U);
    // }
    // else
    //     mask.convertTo(maskMat);

    // const int lbins = 10;
    const int abins = 12;
    const int bbins = 12;
    // const int histSize = {lbins, abins, bbins};
    const int histSize[] = {abins, bbins};
    // const int channels = {0, 1, 2};
    const int channels[] = {1, 2};
        
    // l varies from [0, 100]
    // const float lranges[] = {0, 101};
    // a varies from [-127, 127]
    const float aranges[] = {-127, 128};
    // b varies from [-127, 127]
    const float branges[] = {-127, 128};

    // const float *ranges[] = {lranges, aranges, branges};
    const float *ranges[] = {aranges, branges};
        
    cv::Mat histMat;
    cv::calcHist(&imat, 1, channels, maskMat,
                 histMat, 2, histSize, ranges, true, false);

    double maxVal = 0;
    cv::minMaxLoc(histMat, 0, &maxVal, 0, 0);

    int scale = 10;
    cv::Mat histImg = cv::Mat::zeros(abins*scale, bbins*scale, CV_8UC3);
    for (int a = 0; a < abins; ++a)
    {
        for (int b = 0; b < bbins; ++b)
        {
            histMat.at<float>(a, b) /= maxVal;
            float binVal = histMat.at<float>(a, b);
            int intensity = cvRound(binVal*255);
            cv::rectangle(histImg, cv::Point(a*scale, b*scale),
                          cv::Point((a+1)*scale-1, (b+1)*scale-1),
                          cv::Scalar::all(intensity), CV_FILLED);
        }
    }

    cv::Mat backproj;
    calcBackProject(&imat, 1, channels, histMat, backproj, ranges, 255, true);

    cv::imshow("source", imat);
    cv::imshow("a-b histogram", histImg);
    cv::imshow("backproj", backproj);
    cv::waitKey();
}
