#include"DivideIntoGrids.h";
#include"ImportantMap.h"
#include<iostream>
using namespace std;
double * DivideIntoGrids::energy;
int DivideIntoGrids::w;
int DivideIntoGrids::h;

std::priority_queue<node> q;
DivideIntoGrids::DivideIntoGrids(void){

}

DivideIntoGrids::~DivideIntoGrids(void){

}
void DivideIntoGrids::setData(IplImage* src){//获取能量值，图像宽与长
	energy = new double[src->width*src->height];
	w = src->width;
	h = src->height;
	ImportantMap::compute_important_map(energy,src);
	cout<<w<<" "<<h<<endl;
}

double DivideIntoGrids::calculate_sum_energy(CvPoint start,CvPoint end){//计算这个矩形能量之和
	double sum_energy = 0.0;
	for(int i = start.y;i<=end.y;i++){
		for(int j=start.x;j<=end.x;j++){
			sum_energy = sum_energy+energy[i*w+j];
		}
	}
	return sum_energy;
}

std::priority_queue<node> DivideIntoGrids::make_grids(CvPoint start,CvPoint end){
	double s = DivideIntoGrids::calculate_sum_energy(start,end);
	//cout<<s<<" ";
	int area = (end.y - start.y) * (end.x - start.x);
    //cout<<area<<" ";
	double E1=0.0,E2 = s,aver_e1 = 0.0,aver_e2 = 0.0;
	int area1 = 0,area2 = area;
	double min_ = 1000000,m = 0.0;
	int ind_x,ind_y;

	if((s<THRESHOLD&&area<R_MAX)||(area<R_MIN)){//如果格子不能再分割了，加入队列中
		node grid;grid.top = start;grid.bottom = end;grid.areas = area;
		//cout<<"grid.top = ("<<grid.top.x<<","<<grid.top.y<<"),grid.bottom = ("<<grid.bottom.x<<","<<grid.bottom.y<<") ";
		cout<<s<<","<<area<<" ";
		q.push(grid);
	}
	else{
		for(int i=start.y+1;i<end.y;i++){//一行一行确定最佳分割线，为什么i<=end.y，不要计算到最后一行，防止area2 = 0
			for(int j=start.x;j<=end.x;j++){
				E1 = E1 + energy[i*w+j];
				E2 = E2 - energy[i*w+j];
			}
			area1 = area1 + (end.x - start.x);
			area2 = area2 - (end.x - start.x);
			aver_e1 = E1/(double)area1;
			aver_e2 = E2/(double)area2;
			double tmp1;int tmp2;
			if(aver_e1 > aver_e2) tmp1 = aver_e1;aver_e1 = aver_e2;aver_e2 = tmp1;
			if(area1 < area2) tmp2 = area1;area1 = area2;area2 = tmp2;
			m = ((aver_e1+δ)/(aver_e2+δ))*((double)area1/(double)area2);
			if(m<min_){
				min_ = m; ind_x = end.x;ind_y = i;
			}
		}
		area1 = 0; area2 = area; E1 = 0.0;E2 = s;
		for(int i = start.x+1 ;i<end.x;i++){//为什么i<=end.x，不要计算到最后一列，防止area2 = 0
			for(int j=start.y ; j<= end.y;j++){
				E1 = E1 + energy[j*w+i];
				E2 = E2 - energy[j*w+i];
			}
			area1 = area1 + (end.y - start.y);
			area2 = area2 - (end.y - start.y);
			aver_e1 = E1/(double)area1;
			aver_e2 = E2/(double)area2;
			double tmp1;int tmp2;
			if(aver_e1 > aver_e2) tmp1 = aver_e1;aver_e1 = aver_e2;aver_e2 = tmp1;
			if(area1 < area2) tmp2 = area1;area1 = area2;area2 = tmp2;
			m = ((aver_e1+δ)/(aver_e2+δ))*((double)area1/(double)area2);
			if(m<min_){
				min_ = m; ind_x = i;ind_y = end.y;
			}
		}
		if(ind_y == end.y){
			CvPoint start1,end1;
			start1.x = start.x;start1.y = start.y; end1.x = ind_x;end1.y = ind_y;
			make_grids(start1,end1);
			CvPoint start2,end2;
			start2.x = ind_x+1;start2.y = start.y; end2.x = end.x;end2.y = end.y;
			make_grids(start2,end2);
		}
		else{
			CvPoint start1,end1;
			start1.x = start.x;start1.y = start.y; end1.x = ind_x;end1.y = ind_y;
			make_grids(start1,end1);
			CvPoint start2,end2;
			start2.x = start.x;start2.y = ind_y+1; end2.x = end.x;end2.y = end.y;
			make_grids(start2,end2);
		}
	}
	return q;
}