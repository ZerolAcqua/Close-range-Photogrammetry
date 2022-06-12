#pragma once

#include "Eigen/Eigen"
#include "Point.h"
#include "stdio.h"
#define ITER_TIMES 100
using namespace Eigen;


class DLT
{
public:
	//����Ƭ���Ƶ㡢�����㡢������������(pix)
	std::map<int, Point2f> mLeftCtrImgCoor;
	std::map<int, Point2f> mLeftUknImgCoor;
	std::map<int, Point2f> mLeftChkImgCoor;
	//����Ƭ���Ƶ㡢�����㡢������������(pix)
	std::map<int, Point2f> mRightCtrImgCoor;
	std::map<int, Point2f> mRightUknImgCoor;
	std::map<int, Point2f> mRightChkImgCoor;

	//������Ƭͬ��������ĸ�������
	std::map<int, Point2f> mLeftCorImgCoor;
	std::map<int, Point2f> mRightCorImgCoor;

	//���Ƶ���﷽����
	std::map<int, Point3f> mCtrObjCoor;
	//������﷽����
	std::map<int, Point3f> mChkObjCoor;
	//��������﷽����
	std::map<int, Point3f> mUknObjCoor;


	//Lϵ������
	Matrix<double, 3, 4> mLeftMatL;
	Matrix<double, 3, 4> mRightMatL;
	//����ϵ��
	double mLeftParamK[2] = { 0 };//k1,k2;
	double mLeftParamP[2] = { 0 };//p1,p2;
	double mRightParamK[2] = { 0 };//k1,k2;
	double mRightParamP[2] = { 0 };//p1,p2;
	//�ڷ�λԪ��x0,y0(pix)
	double mLeft_x0 = 0, mLeft_y0 = 0;
	double mRight_x0 = 0, mRight_y0 = 0;
	
	//����ļ���
	std::string mstr;

public:
	/**
	@brief ���к���
	����Ӱ����Ƶ������������

	@param  leftCtrImgCoor		��Ƭ���Ƶ��������
	@param  rightCtrImgCoor		��Ƭ���Ƶ��������
	*/
	void setCtrImgCoor(std::map<int, Point2f>& leftCtrImgCoor, std::map<int, Point2f>& rightCtrImgCoor);
	/**
	@brief ���к���
	����Ӱ������������������

	@param  leftUknImgCoor		��Ƭ�������������
	@param  rightUknImgCoor		��Ƭ�������������
	*/
	void setUknImgCoor(std::map<int, Point2f>& leftUknImgCoor, std::map<int, Point2f>& rightUknImgCoor);
	/**
	@brief ���к���
	����Ӱ����������������

	@param  leftChkImgCoor		��Ƭ�����������
	@param  rightChkImgCoor		��Ƭ�����������
	*/
	void setChkImgCoor(std::map<int, Point2f>& leftChkImgCoor, std::map<int, Point2f>& rightChkImgCoor);
	/**
	@brief ���к���
	���ÿ��Ƶ���﷽����

	@param  ctrObjCoor		���Ƶ���﷽����
	*/
	void setCtrObjCoor(std::map<int, Point3f> ctrObjCoor);
	/**
	@brief ���к���
	���ü�����﷽����

	@param  chkObjCoor		���Ƶ���﷽����
	*/
	void setChkObjCoor(std::map<int, Point3f> chkObjCoor);
	/**
	@brief ���к���
	��������ļ���ǰ׺

	@param  str		����ļ���ǰ׺
	*/
	void setStr(std::string str);
	/**
	@brief ���к���
	�����ʼϵ�����ڷ�λԪ��
	*/
	void initLValue();
	/**
	@brief ���к���
	����Lϵ�����൱�ں󷽽���
	*/
	void calculateLValue();

	/*
	@brief ���к���
	����������꣬��ɸѡ��ͬ����
	*/
	void correctImgCoor();

	/**
	@brief ���к���
	�����ʼ�������﷽����
	*/
	void initUknObjCoorValue();

	/**
	@brief ���к���
	����������﷽���꣬�൱��ǰ������
	*/
	void calculateUknObjCoor();
};