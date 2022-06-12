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
	@brief ���к���
	����Ӱ���ⷽλԪ�س�ֵ
	
	@param EOP Ӱ���ⷽλԪ��
	@param size ���鳤��
	*/
	void setEOP(double* EOP, int size);
	/**
	@brief ���к���
	����Ӱ���ڷ�λԪ�س�ֵ
	
	@param  IOP Ӱ���ڷ�λԪ��
	@param  size ���鳤��
	*/
	void setIOP(double* IOP, int size);
	/**
	@brief ���к���
	����Ӱ����������Ϳ��Ƶ�����

	@param  imgCoor �������
	@param  objCoor ���Ƶ�����
	*/
	void setCoor(std::map<int, Point2f>& imgCoor, std::map<int, Point3f>& objCoor);
	/**
	@brief ���к���
	���㷽λԪ�غͻ���ϵ��

	@param  str ����ļ���
	*/
	void calculate(std::string str);
};