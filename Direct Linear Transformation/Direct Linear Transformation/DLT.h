#pragma once

#include "Eigen/Eigen"
#include "Point.h"
#include "stdio.h"
#define ITER_TIMES 100
using namespace Eigen;


class DLT
{
public:
	//左像片控制点、待定点、检查点的像点坐标(pix)
	std::map<int, Point2f> mLeftCtrImgCoor;
	std::map<int, Point2f> mLeftUknImgCoor;
	std::map<int, Point2f> mLeftChkImgCoor;
	//右像片控制点、待定点、检查点的像点坐标(pix)
	std::map<int, Point2f> mRightCtrImgCoor;
	std::map<int, Point2f> mRightUknImgCoor;
	std::map<int, Point2f> mRightChkImgCoor;

	//左右像片同名待定点的改正坐标
	std::map<int, Point2f> mLeftCorImgCoor;
	std::map<int, Point2f> mRightCorImgCoor;

	//控制点的物方坐标
	std::map<int, Point3f> mCtrObjCoor;
	//检查点的物方坐标
	std::map<int, Point3f> mChkObjCoor;
	//待定点的物方坐标
	std::map<int, Point3f> mUknObjCoor;


	//L系数矩阵
	Matrix<double, 3, 4> mLeftMatL;
	Matrix<double, 3, 4> mRightMatL;
	//畸变系数
	double mLeftParamK[2] = { 0 };//k1,k2;
	double mLeftParamP[2] = { 0 };//p1,p2;
	double mRightParamK[2] = { 0 };//k1,k2;
	double mRightParamP[2] = { 0 };//p1,p2;
	//内方位元素x0,y0(pix)
	double mLeft_x0 = 0, mLeft_y0 = 0;
	double mRight_x0 = 0, mRight_y0 = 0;
	
	//输出文件名
	std::string mstr;

public:
	/**
	@brief 公有函数
	设置影像控制点像点量测坐标

	@param  leftCtrImgCoor		左片控制点像点坐标
	@param  rightCtrImgCoor		右片控制点像点坐标
	*/
	void setCtrImgCoor(std::map<int, Point2f>& leftCtrImgCoor, std::map<int, Point2f>& rightCtrImgCoor);
	/**
	@brief 公有函数
	设置影像待定点像点量测坐标

	@param  leftUknImgCoor		左片待定点像点坐标
	@param  rightUknImgCoor		右片待定点像点坐标
	*/
	void setUknImgCoor(std::map<int, Point2f>& leftUknImgCoor, std::map<int, Point2f>& rightUknImgCoor);
	/**
	@brief 公有函数
	设置影像检查点像点量测坐标

	@param  leftChkImgCoor		左片检查点像点坐标
	@param  rightChkImgCoor		右片检查点像点坐标
	*/
	void setChkImgCoor(std::map<int, Point2f>& leftChkImgCoor, std::map<int, Point2f>& rightChkImgCoor);
	/**
	@brief 公有函数
	设置控制点的物方坐标

	@param  ctrObjCoor		控制点的物方坐标
	*/
	void setCtrObjCoor(std::map<int, Point3f> ctrObjCoor);
	/**
	@brief 公有函数
	设置检查点的物方坐标

	@param  chkObjCoor		控制点的物方坐标
	*/
	void setChkObjCoor(std::map<int, Point3f> chkObjCoor);
	/**
	@brief 公有函数
	设置输出文件名前缀

	@param  str		输出文件名前缀
	*/
	void setStr(std::string str);
	/**
	@brief 公有函数
	计算初始系数和内方位元素
	*/
	void initLValue();
	/**
	@brief 公有函数
	计算L系数，相当于后方交会
	*/
	void calculateLValue();

	/*
	@brief 公有函数
	改正像点坐标，并筛选出同名点
	*/
	void correctImgCoor();

	/**
	@brief 公有函数
	计算初始待定点物方坐标
	*/
	void initUknObjCoorValue();

	/**
	@brief 公有函数
	计算待定点物方坐标，相当于前方交会
	*/
	void calculateUknObjCoor();
};