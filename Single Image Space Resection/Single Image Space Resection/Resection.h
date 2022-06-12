#pragma once

#include "Eigen/Eigen"
#include "Point.h"
#include "stdio.h"
#define ITER_TIMES 100

using namespace Eigen;

class Resection
{
public:
	std::map<int, Point2f> mImgCoor;
	std::map<int, Point3f> mObjCoor;
	double mEOP[6] = { 0 };//Xs Ys Zs phi omega kappa(mm,rad)
	double mIOP[3] = { 0 };//f x0 y0(pix)
	double mParamK[2] = { 0 };//k1,k2;
	double mParamP[2] = { 0 };//p1,p2;

public:
	/**
	@brief 公有函数
	设置影像外方位元素初值
	
	@param EOP 影像外方位元素
	@param size 数组长度
	*/
	void setEOP(double* EOP, int size);
	/**
	@brief 公有函数
	设置影像内方位元素初值
	
	@param  IOP 影像内方位元素
	@param  size 数组长度
	*/
	void setIOP(double* IOP, int size);
	/**
	@brief 公有函数
	设置影像量测坐标和控制点坐标

	@param  imgCoor 像点坐标
	@param  objCoor 控制点坐标
	*/
	void setCoor(std::map<int, Point2f>& imgCoor, std::map<int, Point3f>& objCoor);
	/**
	@brief 公有函数
	计算方位元素和畸变系数

	@param  str 输出文件名
	*/
	void calculate(std::string str);
};