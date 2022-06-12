#include "DLT.h"

#include "DLT.h"
#include "math.h"
#include <iostream>
#include <fstream> 
#include <istream>


void DLT::setCtrImgCoor(std::map<int, Point2f>& leftCtrImgCoor, std::map<int, Point2f>& rightCtrImgCoor)
{
	this->mLeftCtrImgCoor = leftCtrImgCoor;
	this->mRightCtrImgCoor = rightCtrImgCoor;
}
void DLT::setUknImgCoor(std::map<int, Point2f>& leftUknImgCoor, std::map<int, Point2f>& rightUknImgCoor)
{
	this->mLeftUknImgCoor = leftUknImgCoor;
	this->mRightUknImgCoor = rightUknImgCoor;
}
void DLT::setChkImgCoor(std::map<int, Point2f>& leftChkImgCoor, std::map<int, Point2f>& rightChkImgCoor)
{
	this->mLeftChkImgCoor = leftChkImgCoor;
	this->mRightChkImgCoor = rightChkImgCoor;
}
void DLT::setCtrObjCoor(std::map<int, Point3f> ctrObjCoor)
{
	this->mCtrObjCoor = ctrObjCoor;
}
void DLT::setChkObjCoor(std::map<int, Point3f> chkObjCoor)
{
	this->mChkObjCoor = chkObjCoor;
}
void DLT::setStr(std::string str)
{
	this->mstr = str;
}

void DLT::initLValue()
{
	MatrixXd matA();
	int ID;
	double X = 0, Y = 0, Z = 0;
	double x = 0, y = 0;
	int num = 0;
	int i = 0;//计数变量
	double L1, L2, L3, L4;
	double L5, L6, L7, L8;
	double L9, L10, L11;
	//左像片近似值计算
	{
		//计算L系数
		num = mLeftCtrImgCoor.size();
		MatrixXd matA = MatrixXd::Zero(2 * num, 11);
		MatrixXd matL = MatrixXd::Zero(2 * num, 1);
		i = 0;
		for (auto it = this->mLeftCtrImgCoor.begin(); it != mLeftCtrImgCoor.end(); it++)
		{
			ID = it->first;
			x = it->second.mx;
			y = it->second.my;
			X = this->mCtrObjCoor.at(ID).mx;
			Y = this->mCtrObjCoor.at(ID).my;
			Z = this->mCtrObjCoor.at(ID).mz;

			matA(i * 2, 0) = X;
			matA(i * 2, 1) = Y;
			matA(i * 2, 2) = Z;
			matA(i * 2, 3) = 1;
			matA(i * 2, 4) = 0;
			matA(i * 2, 5) = 0;
			matA(i * 2, 6) = 0;
			matA(i * 2, 7) = 0;
			matA(i * 2, 8) = x * X;
			matA(i * 2, 9) = x * Y;
			matA(i * 2, 10) = x * Z;

			matA(i * 2 + 1, 0) = 0;
			matA(i * 2 + 1, 1) = 0;
			matA(i * 2 + 1, 2) = 0;
			matA(i * 2 + 1, 3) = 0;
			matA(i * 2 + 1, 4) = X;
			matA(i * 2 + 1, 5) = Y;
			matA(i * 2 + 1, 6) = Z;
			matA(i * 2 + 1, 7) = 1;
			matA(i * 2 + 1, 8) = y * X;
			matA(i * 2 + 1, 9) = y * Y;
			matA(i * 2 + 1, 10) = y * Z;

			matL(i * 2, 0) = -x;
			matL(i * 2 + 1, 0) = -y;

			i++;
		}
		MatrixXd temp = (matA.transpose() * matA).inverse() * matA.transpose() * matL;
		L1 = this->mLeftMatL(0, 0) = temp(0, 0); L2 = this->mLeftMatL(0, 1) = temp(1, 0); L3 = this->mLeftMatL(0, 2) = temp(2, 0); L4 = this->mLeftMatL(0, 3) = temp(3, 0);
		L5 = this->mLeftMatL(1, 0) = temp(4, 0); L6 = this->mLeftMatL(1, 1) = temp(5, 0); L7 = this->mLeftMatL(1, 2) = temp(6, 0); L8 = this->mLeftMatL(1, 3) = temp(7, 0);
		L9 = this->mLeftMatL(2, 0) = temp(8, 0); L10 = this->mLeftMatL(2, 1) = temp(9, 0); L11 = this->mLeftMatL(2, 2) = temp(10, 0); this->mLeftMatL(2, 3) = 1;

		//计算内方位元素x0,y0
		this->mLeft_x0 = -(L1 * L9 + L2 * L10 + L3 * L11) / (L9 * L9 + L10 * L10 + L11 * L11);
		this->mLeft_y0 = -(L5 * L9 + L6 * L10 + L7 * L11) / (L9 * L9 + L10 * L10 + L11 * L11);

	}
	//右像片近似值计算
	{
		//计算L系数
		num = mRightCtrImgCoor.size();
		MatrixXd matA = MatrixXd::Zero(2 * num, 11);
		MatrixXd matL = MatrixXd::Zero(2 * num, 1);
		i = 0;
		for (auto it = this->mRightCtrImgCoor.begin(); it != mRightCtrImgCoor.end(); it++)
		{
			ID = it->first;
			x = it->second.mx;
			y = it->second.my;
			X = this->mCtrObjCoor.at(ID).mx;
			Y = this->mCtrObjCoor.at(ID).my;
			Z = this->mCtrObjCoor.at(ID).mz;

			matA(i * 2, 0) = X;
			matA(i * 2, 1) = Y;
			matA(i * 2, 2) = Z;
			matA(i * 2, 3) = 1;
			matA(i * 2, 4) = 0;
			matA(i * 2, 5) = 0;
			matA(i * 2, 6) = 0;
			matA(i * 2, 7) = 0;
			matA(i * 2, 8) = x * X;
			matA(i * 2, 9) = x * Y;
			matA(i * 2, 10) = x * Z;

			matA(i * 2 + 1, 0) = 0;
			matA(i * 2 + 1, 1) = 0;
			matA(i * 2 + 1, 2) = 0;
			matA(i * 2 + 1, 3) = 0;
			matA(i * 2 + 1, 4) = X;
			matA(i * 2 + 1, 5) = Y;
			matA(i * 2 + 1, 6) = Z;
			matA(i * 2 + 1, 7) = 1;
			matA(i * 2 + 1, 8) = y * X;
			matA(i * 2 + 1, 9) = y * Y;
			matA(i * 2 + 1, 10) = y * Z;

			matL(i * 2, 0) = -x;
			matL(i * 2 + 1, 0) = -y;

			i++;
		}
		MatrixXd temp = (matA.transpose() * matA).inverse() * matA.transpose() * matL;
		L1 = this->mRightMatL(0, 0) = temp(0, 0); L2 = this->mRightMatL(0, 1) = temp(1, 0); L3 = this->mRightMatL(0, 2) = temp(2, 0); L4 = this->mRightMatL(0, 3) = temp(3, 0);
		L5 = this->mRightMatL(1, 0) = temp(4, 0); L6 = this->mRightMatL(1, 1) = temp(5, 0); L7 = this->mRightMatL(1, 2) = temp(6, 0); L8 = this->mRightMatL(1, 3) = temp(7, 0);
		L9 = this->mRightMatL(2, 0) = temp(8, 0); L10 = this->mRightMatL(2, 1) = temp(9, 0); L11 = this->mRightMatL(2, 2) = temp(10, 0); this->mRightMatL(2, 3) = 1;

		//计算内方位元素x0,y0
		this->mRight_x0 = -(L1 * L9 + L2 * L10 + L3 * L11) / (L9 * L9 + L10 * L10 + L11 * L11);
		this->mRight_y0 = -(L5 * L9 + L6 * L10 + L7 * L11) / (L9 * L9 + L10 * L10 + L11 * L11);

	}
}
void DLT::calculateLValue()
{
	//计算L系数和x0,y0初始值
	initLValue();

	std::ofstream outfile;   //输出流
	
	//---- 左像片“后方交会” ----
	{
		outfile.open("output/" + this->mstr + "_LeftL.rep", std::ios::trunc);
		if (!outfile.is_open())
		{
			std::cout << "output failed" << std::endl;
			return;
		}

		//控制点个数
		int num = this->mLeftCtrImgCoor.size();
		//计数变量
		int count = 0, i = 0;
		//内方位元素
		double x0 = this->mLeft_x0, y0 = this->mLeft_y0;
		//物方点的物方坐标
		double X, Y, Z;
		//像点坐标
		double x, y;
		//L系数
		double L1, L2, L3, L4;
		double L5, L6, L7, L8;
		double L9, L10, L11;
		//分母A
		double A = 1;
		// r = 径向
		double r;
		//畸变系数
		double k1 = this->mLeftParamK[0], k2 = this->mLeftParamK[1];
		double p1 = this->mLeftParamP[0], p2 = this->mLeftParamP[1];
		//法方程系数阵
		MatrixXd matA = MatrixXd::Zero(2 * num, 15);
		MatrixXd matX = MatrixXd::Zero(15, 1);	//L1 ... L11 k1 k2 p1 p2
		//法方程常数项
		MatrixXd matL = MatrixXd::Zero(2 * num, 1);
		//记录上一次迭代的fx的数值(pix)
		double tmpOld, tmpNew;
		double tmpA;
		double tmpB;
		double tmpC;
		double gamma3sq;

		x0 = this->mLeft_x0, y0 = this->mLeft_y0;
		L1 = this->mLeftMatL(0, 0), L2 = this->mLeftMatL(0, 1), L3 = this->mLeftMatL(0, 2), L4 = this->mLeftMatL(0, 3);
		L5 = this->mLeftMatL(1, 0), L6 = this->mLeftMatL(1, 1), L7 = this->mLeftMatL(1, 2), L8 = this->mLeftMatL(1, 3);
		L9 = this->mLeftMatL(2, 0), L10 = this->mLeftMatL(2, 1), L11 = this->mLeftMatL(2, 2);
		
		outfile << "======== 初值计算结果 ========" << std::endl;
		outfile << "L系数" << std::endl;
		outfile << L1 << " " << L2 << " " << L3 << " " << L4 << " " << std::endl;
		outfile << L5 << " " << L6 << " " << L7 << " " << L8 << " " << std::endl;
		outfile << L9 << " " << L10 << " " << L11 << " " << std::endl;
		outfile << "畸变参数" << std::endl;
		outfile <<"k1 " << k1 << " k2 " << k2 << " p1 " << p1 << " p2 " << p2 << " " << std::endl;
		outfile << "像主点坐标(pix)" << std::endl;
		outfile << "x0 " << x0 << " y0 " << y0 << " " << std::endl;

		gamma3sq = 1 / (L9 * L9 + L10 * L10 + L11 * L11);
		tmpA = gamma3sq * (L1 * L1 + L2 * L2 + L3 * L3) - x0 * x0;
		tmpB = gamma3sq * (L5 * L5 + L6 * L6 + L7 * L7) - y0 * y0;
		tmpC = gamma3sq * (L1 * L5 + L2 * L6 + L3 * L7) - x0 * y0;
		tmpOld = sqrt((tmpA * tmpB - tmpC * tmpC) / tmpB);

		outfile << "======== 迭代计算结果 ========" << std::endl;
		for (i = 0; i < ITER_TIMES; i++)
		{
			count = 0;
			//对每个像点列方程
			for (auto it = this->mLeftCtrImgCoor.begin(); it != mLeftCtrImgCoor.end(); it++)
			{
				x = it->second.mx;
				y = it->second.my;
				// 对每一个像点搜索它的物方坐标
				X = this->mCtrObjCoor.at(it->first).mx;
				Y = this->mCtrObjCoor.at(it->first).my;
				Z = this->mCtrObjCoor.at(it->first).mz;

				r = sqrt(pow(x - x0, 2) + pow(y - y0, 2));
				A = L9 * X + L10 * Y + L11 * Z + 1;


				matL(2 * count, 0) = x / A;
				matL(2 * count + 1, 0) = y / A;

				matA(2 * count, 0) = X / A;
				matA(2 * count, 1) = Y / A;
				matA(2 * count, 2) = Z / A;
				matA(2 * count, 3) = 1 / A;
				matA(2 * count, 4) = 0;
				matA(2 * count, 5) = 0;
				matA(2 * count, 6) = 0;
				matA(2 * count, 7) = 0;
				matA(2 * count, 8) = x * X / A;
				matA(2 * count, 9) = x * Y / A;
				matA(2 * count, 10) = x * Z / A;
				matA(2 * count, 11) = (x - x0) * pow(r, 2);
				matA(2 * count, 12) = (x - x0) * pow(r, 4);
				matA(2 * count, 13) = pow(r, 2) + 2 * pow(x - x0, 2);
				matA(2 * count, 14) = 2 * (x - x0) * (y - y0);

				matA(2 * count + 1, 0) = 0;
				matA(2 * count + 1, 1) = 0;
				matA(2 * count + 1, 2) = 0;
				matA(2 * count + 1, 3) = 0;
				matA(2 * count + 1, 4) = X / A;
				matA(2 * count + 1, 5) = Y / A;
				matA(2 * count + 1, 6) = Z / A;
				matA(2 * count + 1, 7) = 1 / A;
				matA(2 * count + 1, 8) = y * X / A;
				matA(2 * count + 1, 9) = y * Y / A;
				matA(2 * count + 1, 10) = y * Z / A;
				matA(2 * count + 1, 11) = (y - y0) * pow(r, 2);
				matA(2 * count + 1, 12) = (y - y0) * pow(r, 4);
				matA(2 * count + 1, 13) = 2 * (x - x0) * (y - y0);
				matA(2 * count + 1, 14) = pow(r, 2) + 2 * pow(y - y0, 2);

				count++;
			}
			matA = -matA;
			//计算新值
			matX = (matA.transpose() * matA).inverse() * matA.transpose() * matL;

			//更新未知数 L1 ... L11 k1 k2 p1 p2 
			L1 = matX(0, 0);
			L2 = matX(1, 0);
			L3 = matX(2, 0);
			L4 = matX(3, 0);
			L5 = matX(4, 0);
			L6 = matX(5, 0);
			L7 = matX(6, 0);
			L8 = matX(7, 0);
			L9 = matX(8, 0);
			L10 = matX(9, 0);
			L11 = matX(10, 0);
			k1 = matX(11, 0);
			k2 = matX(12, 0);
			p1 = matX(13, 0);
			p2 = matX(14, 0);
			//更新内方位元素
			x0 = -(L1 * L9 + L2 * L10 + L3 * L11) / (L9 * L9 + L10 * L10 + L11 * L11);
			y0 = -(L5 * L9 + L6 * L10 + L7 * L11) / (L9 * L9 + L10 * L10 + L11 * L11);




			gamma3sq = 1 / (L9 * L9 + L10 * L10 + L11 * L11);
			tmpA = gamma3sq * (L1 * L1 + L2 * L2 + L3 * L3) - x0 * x0;
			tmpB = gamma3sq * (L5 * L5 + L6 * L6 + L7 * L7) - y0 * y0;
			tmpC = gamma3sq * (L1 * L5 + L2 * L6 + L3 * L7) - x0 * y0;
			tmpNew = sqrt((tmpA * tmpB - tmpC * tmpC) / tmpB);

			//输出
			outfile << "---- 第" << i + 1 << "次迭代 ----" << std::endl;
			outfile << "L系数" << std::endl;
			outfile << L1 << " " << L2 << " " << L3 << " " << L4 << " " << std::endl;
			outfile << L5 << " " << L6 << " " << L7 << " " << L8 << " " << std::endl;
			outfile << L9 << " " << L10 << " " << L11 << " " << std::endl;
			outfile << "畸变参数" << std::endl;
			outfile << "k1 " << k1 << " k2 " << k2 << " p1 " << p1 << " p2 " << p2 << " " << std::endl;
			outfile << "像主点坐标(pix)" << std::endl;
			outfile << "x0 " << x0 << " y0 " << y0 << " " << std::endl;
			outfile << "dfx(pix) " << tmpNew - tmpOld << std::endl;

			if (fabs(tmpNew - tmpOld) < 1e-4)
			{
				break;
			}
			tmpOld = tmpNew;
		}
		if (i < ITER_TIMES)
		{
			outfile << "======== 迭代成功！" << "共迭代" << i + 1 << "次 ======== " << std::endl;
			outfile << "参与平差的控制点个数:" << count << std::endl;
			MatrixXd matV = matA * matX - matL;
			outfile << "像点残差(pix):" << std::endl;
			int j = 0;
			for (auto it = this->mLeftCtrImgCoor.begin(); it != mLeftCtrImgCoor.end(); it++, j++)
			{
				outfile << it->first << " vx: " << matV(j * 2, 0) << " vy: " << matV(j * 2 + 1, 0) << std::endl;
			}
			outfile << "======== 像点最或是值(pix) ========" << std::endl;
			j = 0;
			for (auto it = this->mLeftCtrImgCoor.begin(); it != mLeftCtrImgCoor.end(); it++, j++)
			{
				outfile << it->first << " " << matV(j * 2, 0) + it->second.mx << " " << matV(j * 2 + 1, 0) + it->second.my << std::endl;
			}

			outfile << "======== 精度统计 ========" << std::endl;
			// 单位权中误差
			int r = num * 2 - 15;
			double sigma = sqrt((matV.transpose() * matV)(0, 0) / r);
			outfile << "单位权中误差：" << sigma << std::endl;
			//协因数阵
			MatrixXd matQxhxh = (matA.transpose() * matA).inverse();
			//未知数中误差
			MatrixXd matMxhxh = matQxhxh.diagonal().cwiseSqrt() * sigma;
			//std::cout << "diag(Qxhxh)" << std::endl << matQxhxh.diagonal() << std::endl;
			outfile << "未知数中误差：" << std::endl;
			outfile << "L系数" << std::endl;
			outfile << matMxhxh(0, 0) << " " << matMxhxh(1, 0) << " " << matMxhxh(2, 0) << " " << matMxhxh(3, 0) << " " << std::endl;
			outfile << matMxhxh(4, 0) << " " << matMxhxh(5, 0) << " " << matMxhxh(6, 0) << " " << matMxhxh(7, 0) << " " << std::endl;
			outfile << matMxhxh(8, 0) << " " << matMxhxh(9, 0) << " " << matMxhxh(10, 0) << " " << std::endl;
			outfile << "畸变参数" << std::endl;
			outfile << "k1 " << matMxhxh(11, 0) << " k2 " << matMxhxh(12, 0) << " p1 " << matMxhxh(13, 0) << " p2 " << matMxhxh(14, 0) << " " << std::endl;


			this->mLeftMatL(0, 0) = L1, this->mLeftMatL(0, 1) = L2, this->mLeftMatL(0, 2) = L3, this->mLeftMatL(0, 3) = L4;
			this->mLeftMatL(1, 0) = L5, this->mLeftMatL(1, 1) = L6, this->mLeftMatL(1, 2) = L7, this->mLeftMatL(1, 3) = L8;
			this->mLeftMatL(2, 0) = L9, this->mLeftMatL(2, 1) = L10, this->mLeftMatL(2, 2) = L11;
			this->mLeft_x0 = x0, this->mLeft_y0 = y0;

			this->mLeftParamK[0] = k1, this->mLeftParamK[1] = k2;
			this->mLeftParamP[0] = p1, this->mLeftParamP[1] = p2;

			outfile << "======== 计算结果 ========" << std::endl;
			outfile << "L系数" << std::endl;
			outfile << L1 << " " << L2 << " " << L3 << " " << L4 << " " << std::endl;
			outfile << L5 << " " << L6 << " " << L7 << " " << L8 << " " << std::endl;
			outfile << L9 << " " << L10 << " " << L11 << " " << std::endl;
			outfile << "畸变参数" << std::endl;
			outfile << "k1 " << k1 << " k2 " << k2 << " p1 " << p1 << " p2 " << p2 << " " << std::endl;
			outfile << "像主点坐标(pix)" << std::endl;
			outfile << "x0 " << x0 << " y0 " << y0 << " " << std::endl;

			{
				outfile << "外方位线元素(mm)" << std::endl;
				MatrixXd tmpXYZ(3, 1);
				Vector3d tmpL(3, 1);
				tmpL(0,0) = -L4; tmpL(1,0) = -L8; tmpL(2,0) = -1;
				Matrix3d tmpA;
				tmpA(0, 0) = L1; tmpA(0, 1) = L2; tmpA(0, 2) = L3;
				tmpA(1, 0) = L5; tmpA(1, 1) = L6; tmpA(1, 2) = L7;
				tmpA(2, 0) = L9; tmpA(2, 1) = L10; tmpA(2, 2) = L11;
				tmpXYZ = tmpA.inverse() * tmpL;
				outfile << "Xs " << tmpXYZ(0,0) << " Ys " << tmpXYZ(1,0) << " Zs " << tmpXYZ(2,0) << std::endl;

	
				double gamma3sq = 1 / (L9 * L9 + L10 * L10 + L11 * L11);
				double A = gamma3sq * (L1 * L1 + L2 * L2 + L3 * L3) - x0 * x0;
				double B = gamma3sq * (L5 * L5 + L6 * L6 + L7 * L7) - y0 * y0;
				double C = gamma3sq * (L1 * L5 + L2 * L6 + L3 * L7) - x0 * y0;
				double fx = sqrt((A * B - C * C) / B);
				outfile << "主距fx(pix) " << fx << std::endl;

				double dbeta = C < 0 ? asin(sqrt(C * C / (A * B))) : asin(-sqrt(C * C / (A * B)));
				double ds = sqrt(A / B) - 1;
				outfile << "比例尺不一系数ds " << ds << " 不正交性角dbeta " << dbeta << std::endl;


				outfile << "外方位角元素(rad)" << std::endl;
				double a3 = sqrt(gamma3sq) * L9;
				double b3 = sqrt(gamma3sq) * L10;
				double c3 = sqrt(gamma3sq) * L11;
				double b2 = sqrt(gamma3sq) * (L6 + L10 * y0) * (1 + ds) * cos(dbeta) / fx;
				double b1 = (L2 * sqrt(gamma3sq) + b3 * x0 + b2 * fx * tan(dbeta)) / fx;

				outfile << "phi " << atan(-a3 / c3) << " omega " << atan(-b3) << " kappa " << atan(b1 / b2) << std::endl;

			}
		}
		outfile.close();
	}

	//---- 右像片“后方交会” ----
	{
		outfile.open("output/" + this->mstr + "_RightL.rep", std::ios::trunc);
		if (!outfile.is_open())
		{
			std::cout << "output failed" << std::endl;
			return;
		}

		//控制点个数
		int num = this->mRightCtrImgCoor.size();
		//计数变量
		int count = 0, i = 0;
		//内方位元素
		double x0 = this->mRight_x0, y0 = this->mRight_y0;
		//物方点的物方坐标
		double X, Y, Z;
		//像点坐标
		double x, y;
		//L系数
		double L1, L2, L3, L4;
		double L5, L6, L7, L8;
		double L9, L10, L11;
		//分母A
		double A = 1;
		// r = 径向
		double r;
		//畸变系数
		double k1 = this->mRightParamK[0], k2 = this->mRightParamK[1];
		double p1 = this->mRightParamP[0], p2 = this->mRightParamP[1];
		//法方程系数阵
		MatrixXd matA = MatrixXd::Zero(2 * num, 15);
		MatrixXd matX = MatrixXd::Zero(15, 1);	//L1 ... L11 k1 k2 p1 p2
		//法方程常数项
		MatrixXd matL = MatrixXd::Zero(2 * num, 1);
		//记录上一次迭代的fx的数值(pix)
		double tmpOld, tmpNew;
		double tmpA;
		double tmpB;
		double tmpC;
		double gamma3sq;

		x0 = this->mRight_x0, y0 = this->mRight_y0;
		L1 = this->mRightMatL(0, 0), L2 = this->mRightMatL(0, 1), L3 = this->mRightMatL(0, 2), L4 = this->mRightMatL(0, 3);
		L5 = this->mRightMatL(1, 0), L6 = this->mRightMatL(1, 1), L7 = this->mRightMatL(1, 2), L8 = this->mRightMatL(1, 3);
		L9 = this->mRightMatL(2, 0), L10 = this->mRightMatL(2, 1), L11 = this->mRightMatL(2, 2);

		outfile << "======== 初值计算结果 ========" << std::endl;
		outfile << "L系数" << std::endl;
		outfile << L1 << " " << L2 << " " << L3 << " " << L4 << " " << std::endl;
		outfile << L5 << " " << L6 << " " << L7 << " " << L8 << " " << std::endl;
		outfile << L9 << " " << L10 << " " << L11 << " " << std::endl;
		outfile << "畸变参数" << std::endl;
		outfile << "k1 " << k1 << " k2 " << k2 << " p1 " << p1 << " p2 " << p2 << " " << std::endl;
		outfile << "像主点坐标(pix)" << std::endl;
		outfile << "x0 " << x0 << " y0 " << y0 << " " << std::endl;

		gamma3sq = 1 / (L9 * L9 + L10 * L10 + L11 * L11);
		tmpA = gamma3sq * (L1 * L1 + L2 * L2 + L3 * L3) - x0 * x0;
		tmpB = gamma3sq * (L5 * L5 + L6 * L6 + L7 * L7) - y0 * y0;
		tmpC = gamma3sq * (L1 * L5 + L2 * L6 + L3 * L7) - x0 * y0;
		tmpOld = sqrt((tmpA * tmpB - tmpC * tmpC) / tmpB);

		outfile << "======== 迭代计算结果 ========" << std::endl;
		for (i = 0; i < ITER_TIMES; i++)
		{
			count = 0;
			//对每个像点列方程
			for (auto it = this->mRightCtrImgCoor.begin(); it != mRightCtrImgCoor.end(); it++)
			{
				x = it->second.mx;
				y = it->second.my;
				// 对每一个像点搜索它的物方坐标
				X = this->mCtrObjCoor.at(it->first).mx;
				Y = this->mCtrObjCoor.at(it->first).my;
				Z = this->mCtrObjCoor.at(it->first).mz;

				r = sqrt(pow(x - x0, 2) + pow(y - y0, 2));
				A = L9 * X + L10 * Y + L11 * Z + 1;


				matL(2 * count, 0) = x / A;
				matL(2 * count + 1, 0) = y / A;

				matA(2 * count, 0) = X / A;
				matA(2 * count, 1) = Y / A;
				matA(2 * count, 2) = Z / A;
				matA(2 * count, 3) = 1 / A;
				matA(2 * count, 4) = 0;
				matA(2 * count, 5) = 0;
				matA(2 * count, 6) = 0;
				matA(2 * count, 7) = 0;
				matA(2 * count, 8) = x * X / A;
				matA(2 * count, 9) = x * Y / A;
				matA(2 * count, 10) = x * Z / A;
				matA(2 * count, 11) = (x - x0) * pow(r, 2);
				matA(2 * count, 12) = (x - x0) * pow(r, 4);
				matA(2 * count, 13) = pow(r, 2) + 2 * pow(x - x0, 2);
				matA(2 * count, 14) = 2 * (x - x0) * (y - y0);

				matA(2 * count + 1, 0) = 0;
				matA(2 * count + 1, 1) = 0;
				matA(2 * count + 1, 2) = 0;
				matA(2 * count + 1, 3) = 0;
				matA(2 * count + 1, 4) = X / A;
				matA(2 * count + 1, 5) = Y / A;
				matA(2 * count + 1, 6) = Z / A;
				matA(2 * count + 1, 7) = 1 / A;
				matA(2 * count + 1, 8) = y * X / A;
				matA(2 * count + 1, 9) = y * Y / A;
				matA(2 * count + 1, 10) = y * Z / A;
				matA(2 * count + 1, 11) = (y - y0) * pow(r, 2);
				matA(2 * count + 1, 12) = (y - y0) * pow(r, 4);
				matA(2 * count + 1, 13) = 2 * (x - x0) * (y - y0);
				matA(2 * count + 1, 14) = pow(r, 2) + 2 * pow(y - y0, 2);

				count++;
			}
			matA = -matA;
			//计算新值
			matX = (matA.transpose() * matA).inverse() * matA.transpose() * matL;

			//更新未知数 L1 ... L11 k1 k2 p1 p2 
			L1 = matX(0, 0);
			L2 = matX(1, 0);
			L3 = matX(2, 0);
			L4 = matX(3, 0);
			L5 = matX(4, 0);
			L6 = matX(5, 0);
			L7 = matX(6, 0);
			L8 = matX(7, 0);
			L9 = matX(8, 0);
			L10 = matX(9, 0);
			L11 = matX(10, 0);
			k1 = matX(11, 0);
			k2 = matX(12, 0);
			p1 = matX(13, 0);
			p2 = matX(14, 0);
			//更新内方位元素
			x0 = -(L1 * L9 + L2 * L10 + L3 * L11) / (L9 * L9 + L10 * L10 + L11 * L11);
			y0 = -(L5 * L9 + L6 * L10 + L7 * L11) / (L9 * L9 + L10 * L10 + L11 * L11);

			gamma3sq = 1 / (L9 * L9 + L10 * L10 + L11 * L11);
			tmpA = gamma3sq * (L1 * L1 + L2 * L2 + L3 * L3) - x0 * x0;
			tmpB = gamma3sq * (L5 * L5 + L6 * L6 + L7 * L7) - y0 * y0;
			tmpC = gamma3sq * (L1 * L5 + L2 * L6 + L3 * L7) - x0 * y0;
			tmpNew = sqrt((tmpA * tmpB - tmpC * tmpC) / tmpB);

			//输出
			outfile << "---- 第" << i + 1 << "次迭代 ----" << std::endl;
			outfile << "L系数" << std::endl;
			outfile << L1 << " " << L2 << " " << L3 << " " << L4 << " " << std::endl;
			outfile << L5 << " " << L6 << " " << L7 << " " << L8 << " " << std::endl;
			outfile << L9 << " " << L10 << " " << L11 << " " << std::endl;
			outfile << "畸变参数" << std::endl;
			outfile << "k1 " << k1 << " k2 " << k2 << " p1 " << p1 << " p2 " << p2 << " " << std::endl;
			outfile << "像主点坐标(pix)" << std::endl;
			outfile << "x0 " << x0 << " y0 " << y0 << " " << std::endl;
			outfile << "dfx(pix) " << tmpNew - tmpOld << std::endl;


			if (fabs(tmpNew - tmpOld) < 1e-4)
			{
				break;
			}
			tmpOld = tmpNew;
		}
		if (i < ITER_TIMES)
		{
			outfile << "======== 迭代成功！" << "共迭代" << i + 1 << "次 ======== " << std::endl;
			outfile << "参与平差的控制点个数:" << count << std::endl;
			MatrixXd matV = matA * matX - matL;
			outfile << "像点残差(pix):" << std::endl;
			int j = 0;
			for (auto it = this->mRightCtrImgCoor.begin(); it != mRightCtrImgCoor.end(); it++, j++)
			{
				outfile << it->first << " vx: " << matV(j * 2, 0) << " vy: " << matV(j * 2 + 1, 0) << std::endl;
			}
			outfile << "======== 像点最或是值(pix) ========" << std::endl;
			j = 0;
			for (auto it = this->mRightCtrImgCoor.begin(); it != mRightCtrImgCoor.end(); it++, j++)
			{
				outfile << it->first << " " << matV(j * 2, 0) + it->second.mx << " " << matV(j * 2 + 1, 0) + it->second.my << std::endl;
			}

			outfile << "======== 精度统计 ========" << std::endl;
			// 单位权中误差
			int r = num * 2 - 15;
			double sigma = sqrt((matV.transpose() * matV)(0, 0) / r);
			outfile << "单位权中误差：" << sigma << std::endl;
			//协因数阵
			MatrixXd matQxhxh = (matA.transpose() * matA).inverse();
			//未知数中误差
			MatrixXd matMxhxh = matQxhxh.diagonal().cwiseSqrt() * sigma;
			//std::cout << "diag(Qxhxh)" << std::endl << matQxhxh.diagonal() << std::endl;
			outfile << "未知数中误差：" << std::endl;
			outfile << "L系数" << std::endl;
			outfile << matMxhxh(0, 0) << " " << matMxhxh(1, 0) << " " << matMxhxh(2, 0) << " " << matMxhxh(3, 0) << " " << std::endl;
			outfile << matMxhxh(4, 0) << " " << matMxhxh(5, 0) << " " << matMxhxh(6, 0) << " " << matMxhxh(7, 0) << " " << std::endl;
			outfile << matMxhxh(8, 0) << " " << matMxhxh(9, 0) << " " << matMxhxh(10, 0) << " " << std::endl;
			outfile << "畸变参数" << std::endl;
			outfile << "k1 " << matMxhxh(11, 0) << " k2 " << matMxhxh(12, 0) << " p1 " << matMxhxh(13, 0) << " p2 " << matMxhxh(14, 0) << " " << std::endl;



			this->mRightMatL(0, 0) = L1, this->mRightMatL(0, 1) = L2, this->mRightMatL(0, 2) = L3, this->mRightMatL(0, 3) = L4;
			this->mRightMatL(1, 0) = L5, this->mRightMatL(1, 1) = L6, this->mRightMatL(1, 2) = L7, this->mRightMatL(1, 3) = L8;
			this->mRightMatL(2, 0) = L9, this->mRightMatL(2, 1) = L10, this->mRightMatL(2, 2) = L11;
			this->mRight_x0 = x0, this->mRight_y0 = y0;

			this->mRightParamK[0] = k1, this->mRightParamK[1] = k2;
			this->mRightParamP[0] = p1, this->mRightParamP[1] = p2;

			outfile << "======== 计算结果 ========" << std::endl;
			outfile << "L系数" << std::endl;
			outfile << L1 << " " << L2 << " " << L3 << " " << L4 << " " << std::endl;
			outfile << L5 << " " << L6 << " " << L7 << " " << L8 << " " << std::endl;
			outfile << L9 << " " << L10 << " " << L11 << " " << std::endl;
			outfile << "畸变参数" << std::endl;
			outfile << "k1 " << k1 << " k2 " << k2 << " p1 " << p1 << " p2 " << p2 << " " << std::endl;
			outfile << "像主点坐标(pix)" << std::endl;
			outfile << "x0 " << x0 << " y0 " << y0 << " " << std::endl;

			{
				outfile << "外方位线元素(mm)" << std::endl;
				MatrixXd tmpXYZ(3, 1);
				Vector3d tmpL(3, 1);
				tmpL(0, 0) = -L4; tmpL(1, 0) = -L8; tmpL(2, 0) = -1;
				Matrix3d tmpA;
				tmpA(0, 0) = L1; tmpA(0, 1) = L2; tmpA(0, 2) = L3;
				tmpA(1, 0) = L5; tmpA(1, 1) = L6; tmpA(1, 2) = L7;
				tmpA(2, 0) = L9; tmpA(2, 1) = L10; tmpA(2, 2) = L11;
				tmpXYZ = tmpA.inverse() * tmpL;
				outfile << "Xs " << tmpXYZ(0, 0) << " Ys " << tmpXYZ(1, 0) << " Zs " << tmpXYZ(2, 0) << std::endl;


				double gamma3sq = 1 / (L9 * L9 + L10 * L10 + L11 * L11);
				double A = gamma3sq * (L1 * L1 + L2 * L2 + L3 * L3) - x0 * x0;
				double B = gamma3sq * (L5 * L5 + L6 * L6 + L7 * L7) - y0 * y0;
				double C = gamma3sq * (L1 * L5 + L2 * L6 + L3 * L7) - x0 * y0;
				double fx = sqrt((A * B - C * C) / B);
				outfile << "主距fx(pix) " << fx << std::endl;

				double dbeta = C < 0 ? asin(sqrt(C * C / (A * B))) : asin(-sqrt(C * C / (A * B)));
				double ds = sqrt(A / B) - 1;
				outfile << "比例尺不一系数ds " << ds << " 不正交性角dbeta " << dbeta << std::endl;


				outfile << "外方位角元素(rad)" << std::endl;
				double a3 = sqrt(gamma3sq) * L9;
				double b3 = sqrt(gamma3sq) * L10;
				double c3 = sqrt(gamma3sq) * L11;
				double b2 = sqrt(gamma3sq) * (L6 + L10 * y0) * (1 + ds) * cos(dbeta) / fx;
				double b1 = (L2 * sqrt(gamma3sq) + b3 * x0 + b2 * fx * tan(dbeta)) / fx;

				outfile << "phi " << atan(-a3 / c3) << " omega " << atan(-b3) << " kappa " << atan(b1 / b2) << std::endl;

			}


		}
		outfile.close();
	}
}
void DLT::correctImgCoor()
{
	//合并待定点与检查点的map
	this->mLeftUknImgCoor.insert(this->mLeftChkImgCoor.begin(), this->mLeftChkImgCoor.end());
	this->mRightUknImgCoor.insert(this->mRightChkImgCoor.begin(), this->mRightChkImgCoor.end());

	std::map<int, Point2f> tmpLeft;
	std::map<int, Point2f> tmpRight;

	//先改正
	double x, y;
	double x0, y0;
	double r, dx, dy;
	double k1, k2;
	double p1, p2;
	for (auto it = this->mLeftUknImgCoor.begin(); it != this->mLeftUknImgCoor.end(); it++)
	{
		x = it->second.mx;
		y = it->second.my;
		x0 = this->mLeft_x0;
		y0 = this->mLeft_y0;
		
		r = sqrt(pow(x - x0, 2) + pow(y - y0, 2));
		k1 = this->mLeftParamK[0], k2 = this->mLeftParamK[1];
		p1 = this->mLeftParamP[0], p2 = this->mLeftParamP[1];
		
		dx = k1 * (x - x0) * pow(r, 2) + k2 * (x - x0) * pow(r, 4) + p1 * (pow(r, 2) + 2 * pow(x - x0, 2)) + 2 * p2 * (x - x0) * (y - y0);
		dy = k1 * (y - y0) * pow(r, 2) + k2 * (y - y0) * pow(r, 4) + p2 * (pow(r, 2) + 2 * pow(y - y0, 2)) + 2 * p1 * (x - x0) * (y - y0);

		tmpLeft[it->first] = Point2f(x + dx, y + dy);
	}
	for (auto it = this->mRightUknImgCoor.begin(); it != this->mRightUknImgCoor.end(); it++)
	{
		x = it->second.mx;
		y = it->second.my;
		x0 = this->mRight_x0;
		y0 = this->mRight_y0;

		r = sqrt(pow(x - x0, 2) + pow(y - y0, 2));
		k1 = this->mRightParamK[0], k2 = this->mRightParamK[1];
		p1 = this->mRightParamP[0], p2 = this->mRightParamP[1];

		dx = k1 * (x - x0) * pow(r, 2) + k2 * (x - x0) * pow(r, 4) + p1 * (pow(r, 2) + 2 * pow(x - x0, 2)) + 2 * p2 * (x - x0) * (y - y0);
		dy = k1 * (y - y0) * pow(r, 2) + k2 * (y - y0) * pow(r, 4) + p2 * (pow(r, 2) + 2 * pow(y - y0, 2)) + 2 * p1 * (x - x0) * (y - y0);

		tmpRight[it->first] = Point2f(x + dx, y + dy);
	}

	//再配对
	this->mLeftCorImgCoor.clear();
	this->mRightCorImgCoor.clear();


	auto it1 = tmpLeft.begin();
	auto it2 = tmpRight.begin(); 
	while (it1 != tmpLeft.end() && it2 != tmpRight.end())
	{
		if (it1->first < it2->first)
			it1++;
		else if (it1->first > it2->first)
			it2++;
		else
		{
			this->mLeftCorImgCoor[it1->first] = it1->second;
			this->mRightCorImgCoor[it2->first] = it2->second;
			it1++;
			it2++;
		}
	}
}
void DLT::initUknObjCoorValue()
{
	this->mUknObjCoor.clear();

	MatrixXd matA();
	int ID;
	double Left_x, Left_y;
	double Right_x, Right_y;

	double LeftL1, LeftL2, LeftL3, LeftL4;
	double LeftL5, LeftL6, LeftL7, LeftL8;
	double LeftL9, LeftL10, LeftL11;

	double RightL1, RightL2, RightL3, RightL4;
	double RightL5, RightL6, RightL7, RightL8;
	double RightL9, RightL10, RightL11;

	LeftL1 = this->mLeftMatL(0, 0), LeftL2 = this->mLeftMatL(0, 1), LeftL3 = this->mLeftMatL(0, 2), LeftL4 = this->mLeftMatL(0, 3);
	LeftL5 = this->mLeftMatL(1, 0), LeftL6 = this->mLeftMatL(1, 1), LeftL7 = this->mLeftMatL(1, 2), LeftL8 = this->mLeftMatL(1, 3);
	LeftL9 = this->mLeftMatL(2, 0), LeftL10 = this->mLeftMatL(2, 1), LeftL11 = this->mLeftMatL(2, 2);

	RightL1 = this->mRightMatL(0, 0), RightL2 = this->mRightMatL(0, 1), RightL3 = this->mRightMatL(0, 2), RightL4 = this->mRightMatL(0, 3);
	RightL5 = this->mRightMatL(1, 0), RightL6 = this->mRightMatL(1, 1), RightL7 = this->mRightMatL(1, 2), RightL8 = this->mRightMatL(1, 3);
	RightL9 = this->mRightMatL(2, 0), RightL10 = this->mRightMatL(2, 1), RightL11 = this->mRightMatL(2, 2);


	//待定点物方坐标近似值计算
	{
		MatrixXd matA = MatrixXd::Zero(4, 3);
		MatrixXd matL = MatrixXd::Zero(4, 1);

		for (auto it = this->mLeftCorImgCoor.begin(); it != this->mLeftCorImgCoor.end(); it++)
		{
			ID = it->first;
			
			//左片
			Left_x = it->second.mx;
			Left_y = it->second.my;

			matA(0, 0) = LeftL1 + Left_x * LeftL9;
			matA(0, 1) = LeftL2 + Left_x * LeftL10;
			matA(0, 2) = LeftL3 + Left_x * LeftL11;
			matL(0, 0) = -(LeftL4 + Left_x);

			matA(1, 0) = LeftL5 + Left_y * LeftL9;
			matA(1, 1) = LeftL6 + Left_y * LeftL10;
			matA(1, 2) = LeftL7 + Left_y * LeftL11;
			matL(1, 0) = -(LeftL8 + Left_y);

			//右片
			Right_x = this->mRightCorImgCoor.at(ID).mx;
			Right_y = this->mRightCorImgCoor.at(ID).my;

			matA(2, 0) = RightL1 + Right_x * RightL9;
			matA(2, 1) = RightL2 + Right_x * RightL10;
			matA(2, 2) = RightL3 + Right_x * RightL11;
			matL(2, 0) = -(RightL4 + Right_x);

			matA(3, 0) = RightL5 + Right_y * RightL9;
			matA(3, 1) = RightL6 + Right_y * RightL10;
			matA(3, 2) = RightL7 + Right_y * RightL11;
			matL(3, 0) = -(RightL8 + Right_y);

			MatrixXd temp = (matA.transpose() * matA).inverse() * matA.transpose() * matL;
			this->mUknObjCoor[ID] = Point3f(temp(0, 0), temp(1, 0), temp(2, 0));
		}
	}
}
void DLT::calculateUknObjCoor()
{
	//改正待定点的像点坐标
	correctImgCoor();
	//计算待定点的物方坐标初始值
	initUknObjCoorValue();
	
	std::ofstream outfile;   //输出流
	std::map<int, double> tmpM; //记录单位权中误差
	std::map<int, MatrixXd> tmpMx; //记录未知数中误差

	outfile.open("output/" + this->mstr + "_uknCoor.rep", std::ios::trunc);
	if (!outfile.is_open())
	{
		std::cout << "output failed" << std::endl;
		return;
	}

	//计算待定点的物方坐标精确值
	{
		int ID;
		int i = 0; //计数变量
		double X = 0, Y = 0, Z = 0;
		double Left_x, Left_y;
		double Right_x, Right_y;
		double LeftA, RightA;

		double LeftL1, LeftL2, LeftL3, LeftL4;
		double LeftL5, LeftL6, LeftL7, LeftL8;
		double LeftL9, LeftL10, LeftL11;

		double RightL1, RightL2, RightL3, RightL4;
		double RightL5, RightL6, RightL7, RightL8;
		double RightL9, RightL10, RightL11;

		LeftL1 = this->mLeftMatL(0, 0), LeftL2 = this->mLeftMatL(0, 1), LeftL3 = this->mLeftMatL(0, 2), LeftL4 = this->mLeftMatL(0, 3);
		LeftL5 = this->mLeftMatL(1, 0), LeftL6 = this->mLeftMatL(1, 1), LeftL7 = this->mLeftMatL(1, 2), LeftL8 = this->mLeftMatL(1, 3);
		LeftL9 = this->mLeftMatL(2, 0), LeftL10 = this->mLeftMatL(2, 1), LeftL11 = this->mLeftMatL(2, 2);

		RightL1 = this->mRightMatL(0, 0), RightL2 = this->mRightMatL(0, 1), RightL3 = this->mRightMatL(0, 2), RightL4 = this->mRightMatL(0, 3);
		RightL5 = this->mRightMatL(1, 0), RightL6 = this->mRightMatL(1, 1), RightL7 = this->mRightMatL(1, 2), RightL8 = this->mRightMatL(1, 3);
		RightL9 = this->mRightMatL(2, 0), RightL10 = this->mRightMatL(2, 1), RightL11 = this->mRightMatL(2, 2);


		MatrixXd matA = MatrixXd::Zero(4, 3);
		MatrixXd matL = MatrixXd::Zero(4, 1);
		MatrixXd matX;
		double tmpOld, tmpNew;

		outfile << "======== 逐点迭代计算结果 ========" << std::endl;
		outfile << "共解算" << this->mLeftCorImgCoor.size() << "个物方点，其中待定点" << this->mLeftCorImgCoor.size() - this->mLeftChkImgCoor.size() << "个"
			<< "，检查点" << this->mLeftChkImgCoor.size() << "个" << std::endl;
		for (auto it = this->mLeftCorImgCoor.begin(); it != this->mLeftCorImgCoor.end(); it++)
		{
			ID = it->first;

			Left_x = it->second.mx;
			Left_y = it->second.my;

			Right_x = this->mRightCorImgCoor.at(ID).mx;
			Right_y = this->mRightCorImgCoor.at(ID).my;

			tmpOld = this->mUknObjCoor.at(ID).mx;


			outfile << "==== " << ID << " ====" << std::endl;
			outfile << "初值X,Y,Z(mm)" << std::endl;
			outfile << this->mUknObjCoor.at(ID).mx << " " << this->mUknObjCoor.at(ID).my << " " << this->mUknObjCoor.at(ID).mz << " " << std::endl;

			for (i = 0; i < ITER_TIMES; i++)
			{
				X = this->mUknObjCoor.at(ID).mx;
				Y = this->mUknObjCoor.at(ID).my;
				Z = this->mUknObjCoor.at(ID).mz;

				LeftA = LeftL9 * X + LeftL10 * Y + LeftL11 * Z + 1;
				RightA = RightL9 * X + RightL10 * Y + RightL11 * Z + 1;

				//左片
				matA(0, 0) = -(LeftL1 + Left_x * LeftL9) / LeftA;
				matA(0, 1) = -(LeftL2 + Left_x * LeftL10) / LeftA;
				matA(0, 2) = -(LeftL3 + Left_x * LeftL11) / LeftA;
				matL(0, 0) = (LeftL4 + Left_x) / LeftA;

				matA(1, 0) = -(LeftL5 + Left_y * LeftL9) / LeftA;
				matA(1, 1) = -(LeftL6 + Left_y * LeftL10) / LeftA;
				matA(1, 2) = -(LeftL7 + Left_y * LeftL11) / LeftA;
				matL(1, 0) = (LeftL8 + Left_y) / LeftA;

				//右片
				matA(2, 0) = -(RightL1 + Right_x * RightL9) / RightA;
				matA(2, 1) = -(RightL2 + Right_x * RightL10) / RightA;
				matA(2, 2) = -(RightL3 + Right_x * RightL11) / RightA;
				matL(2, 0) = (RightL4 + Right_x) / RightA;

				matA(3, 0) = -(RightL5 + Right_y * RightL9) / RightA;
				matA(3, 1) = -(RightL6 + Right_y * RightL10) / RightA;
				matA(3, 2) = -(RightL7 + Right_y * RightL11) / RightA;
				matL(3, 0) = (RightL8 + Right_y) / RightA;

				matX = (matA.transpose() * matA).inverse() * matA.transpose() * matL;
				tmpNew = matX(0, 0);

				outfile << "---- 第" << i + 1 << "次迭代 ----" << std::endl;
				outfile << "物方坐标X,Y,Z(mm)" << std::endl;
				outfile << matX(0, 0) << " " << matX(1, 0) << " " << matX(2, 0) << std::endl;

				if (fabs(tmpNew - tmpOld) < 1e-10)
				{
					break;
				}
				tmpOld = tmpNew;

				X = matX(0, 0);
				Y = matX(1, 0);
				Z = matX(2, 0);

			}
			if (i < ITER_TIMES)
			{
				this->mUknObjCoor[ID] = Point3f(matX(0, 0), matX(1, 0), matX(2, 0));
				outfile << "---- "<<"迭代成功！" << "共迭代" << i + 1 << "次 ----" << std::endl;

				MatrixXd matV = matA * matX - matL;

				outfile << "物方坐标X,Y,Z(mm)" << std::endl;
				outfile << matX(0, 0) << " " << matX(1, 0) << " " << matX(2, 0) << std::endl;

				outfile << "左片像点残差vx,vy(pix)" << std::endl;
				outfile << matV(0, 0) << " " << matV(1, 0) << std::endl;
				outfile << "右片像点残差vx,vy(pix)" << std::endl;
				outfile << matV(2, 0) << " " << matV(3, 0) << std::endl;

				// 单位权中误差
				int r = 2 * 2 - 3;
				double sigma = sqrt((matV.transpose() * matV)(0, 0) / r);
				tmpM[ID] = sigma;
				outfile << "单位权中误差：" << sigma << std::endl;

				//协因数阵
				MatrixXd matQxhxh = (matA.transpose() * matA).inverse();
				//未知数中误差
				MatrixXd matMxhxh = matQxhxh.diagonal().cwiseSqrt() * sigma;
				tmpMx[ID] = matMxhxh;

			}
		}
	}
	
	//输出物方坐标
	{
		int ID = 0;
		outfile << "======== 计算结果(mm) ========" << std::endl;
		outfile << "X Y Z m0 mX mY mZ [flag(is checking points)]" << std::endl;
		for (auto it = this->mUknObjCoor.begin(); it != this->mUknObjCoor.end(); it++)
		{
			ID = it->first;
			auto it1 = this->mChkObjCoor.find(ID);
			outfile << ID << " " << it->second.mx << " " << it->second.my << " " << it->second.mz << " " << tmpM.at(ID) << " "
				<< tmpMx.at(ID)(0, 0) << " " << tmpMx.at(ID)(1, 0) << " " << tmpMx.at(ID)(2, 0) << " " << (it1 != this->mChkObjCoor.end()) << std::endl;
		}
	}

	//计算外精度
	{
		int num = 0;
		int ID;
		double sumdX2 = 0, sumdY2 = 0, sumdZ2 = 0;
		outfile << "======== 检查点计算外精度(mm) ========" << std::endl;
		outfile << "ID vX vY vZ" << std::endl;
		for (auto it = this->mChkObjCoor.begin(); it != this->mChkObjCoor.end(); it++)
		{
			ID = it->first;
			auto it1 = this->mUknObjCoor.find(ID);

			if (it1 != this->mUknObjCoor.end())
			{
				outfile<<ID<<" " << it->second.mx - it1->second.mx << " " << it->second.my - it1->second.my << " " << it->second.mz - it1->second.mz<<std::endl;
				sumdX2 += pow(it->second.mx - it1->second.mx, 2);
				sumdY2 += pow(it->second.my - it1->second.my, 2);
				sumdZ2 += pow(it->second.mz - it1->second.mz, 2);
				num++;
				this->mUknObjCoor.erase(it1);

				auto it2 = this->mLeftUknImgCoor.find(ID);
				this->mLeftUknImgCoor.erase(it2);
				auto it3 = this->mRightUknImgCoor.find(ID);
				this->mRightUknImgCoor.erase(it3);
				auto it4 = this->mLeftCorImgCoor.find(ID);
				this->mLeftCorImgCoor.erase(it4);
				auto it5 = this->mRightCorImgCoor.find(ID);
				this->mRightCorImgCoor.erase(it5);
			}
		}
		outfile << "X坐标中误差: " << sqrt(sumdX2 / num) << std::endl;
		outfile << "Y坐标中误差: " << sqrt(sumdY2 / num) << std::endl;
		outfile << "Z坐标中误差: " << sqrt(sumdZ2 / num) << std::endl;
		outfile << "点位中误差： " << sqrt((sumdX2 + sumdY2 + sumdZ2) / num) << std::endl;


	}
	outfile.close();
}


