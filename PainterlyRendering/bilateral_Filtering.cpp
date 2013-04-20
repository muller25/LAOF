#include"bilateral_Filtering.h"
IplImage *  bilateral_filter::bilateral_Filters(char *picturePath){
   IplImage *m_image=cvLoadImage(picturePath);  
   int w=m_image->width;  
   int h=m_image->height;  
     
   IplImage *m_ba=cvCreateImage(cvSize(w,h),IPL_DEPTH_8U,3);  

   for(int i=1;i<h-1;i++)  
      for(int j=1;j<w-1;j++)  
      {  
        
         CvScalar cs[9];  
         CvScalar css = cvScalarAll(0);  
           CvScalar cc=cvGet2D(m_image,i,j);  
         int a=0;  
         double sum=0;  
         double c[9],s[9];  
           
           
         for(int k=i-1;k<i+2;k++)  
            for(int q=j-1;q<j+2;q++)  
            {  
               cs[a]=cvGet2D(m_image,k,q);  
               c[a]=exp(-((pow((double)k-i,2)+pow((double)(q-j),2))/200));  
               s[a]=exp(-(pow(((cs[a].val[0]-cc.val[0])),2)+pow(((cs[a].val[1]-cc.val[1])),2)+pow(((cs[a].val[2]-cc.val[2])),2))/200);  
                 
               a++;  

                }  
         for(int y=0;y<9;y++)  
         {    
            sum=sum+c[y]*s[y];  
         }  

           
             for(int x=0;x<9;x++)  
          {  
             css.val[0]+=(c[x]*s[x]*(cs[x].val[0]))/sum;  
             css.val[1]+=(c[x]*s[x]*(cs[x].val[1]))/sum;  
             css.val[2]+=(c[x]*s[x]*(cs[x].val[2]))/sum;  

          }  

          cvSet2D(m_ba,i,j,css);  

             
          } 
	  cvReleaseImage(&m_image);
	  return m_ba;
}

