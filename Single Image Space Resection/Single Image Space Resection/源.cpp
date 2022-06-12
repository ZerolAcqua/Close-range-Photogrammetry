#include <stdio.h>
#include <istream>
#include <fstream> 
#include <iostream>


#include "Eigen/Eigen"  

#include "point.h"
#include "Resection.h"
#define PI 3.14159265358979

using namespace std;
using namespace Eigen;


void leftHand2RightHand(std::map<int, Point3f>& ptMap);
void rightHand2LeftHand(std::map<int, Point3f>& ptMap);

int main()
{
    /*  共线方程:
        x - x0 + dx = -f * Xbar / Zbar
        y - y0 + dy = -f * Ybar / Zbar

        畸变公式:
        r = sqrt((x - x0)^2 + (y - y0)^2)
        dx = (x - x0)(k1 * r^2 + k2 * r^4) 
           + p1(r^2 + 2(x - x0)^2 + 2 * p2(x - x0)(y - y0)) 
        dy = (y - y0)(k1 * r^2 + k2 * r^4) +
           + p2(r^2 + 2(y - y0)^2 + 2 * p1(x - x0)(y - y0)) 
    */


    string leftImgDataPath = "data/NEW_IMG_8620.JPG.dat";
    string rightImgDataPath = "data/NEW_IMG_8621.JPG.dat";
    string controlDataPath = "data/2022实习-近景控制场-20220520.txt";

    double dpix = 22.2 / 4272;
    int w = 4272;
    int h = 2848;


    std::map<int, Point2f> leftImgCoor;//单位：pix
    std::map<int, Point2f> rightImgCoor;//单位：pix
    std::map<int, Point3f> objCoor;//单位：mm

    ifstream infile;   //输入流

    //读取左影像像点坐标
    {
        int ID = 0;
        double x = 0, y = 0, height = 0, width = 0, angle = 0;
        infile.open(leftImgDataPath, ios::in);
        if (!infile.is_open())
        {
            cout << "left image coordinate file not found." << endl;
        }
        else {
            char temp[128];
            infile.getline(temp, 128);		//跳过注释信息
            while (!infile.eof())           // 若未到文件结束一直循环
            {
                infile >> ID >> x >> y >> width >> height >> angle;
                // 像点坐标量测时是以左上角为原点的，要转换为以左下角为原点的坐标
                if (ID >= 100)
                {
                    leftImgCoor[ID] = Point2f(x, h - 1 - y);
                }
            }
            infile.close();   //关闭文件
        }
    }

    //读取右影像像点坐标
    {
        int ID = 0;
        double x = 0, y = 0, height = 0, width = 0, angle = 0;
        infile.open(rightImgDataPath, ios::in);
        if (!infile.is_open())
        {
            cout << "right image coordinate file not found." << endl;
        }
        else {
            char temp[128];
            infile.getline(temp, 128);		//跳过注释信息
            while (!infile.eof())           // 若未到文件结束一直循环
            {
                infile >> ID >> x >> y >> width >> height >> angle;
                // 像点坐标量测时是以左上角为原点的，要转换为以左下角为原点的坐标
                if (ID >= 100)
                {
                    rightImgCoor[ID] = Point2f(x, h - 1 - y);
                }
            }
            infile.close();   //关闭文件
        }
    }

    //读取控制点坐标
    {
        int ID = 0;
        int num = 0;
        double x = 0, y = 0, z = 0;
        bool flag = false;

        infile.open(controlDataPath, ios::in);
        if (!infile.is_open())
        {
            cout << "control point coordinate file not found." << endl;
        }
        else {
            infile >> num;
            for (int i = 0; i < num; i++)
            {
                infile >> ID >> x >> y >> z >> flag;
                objCoor[ID] = Point3f(x, y, z);               
            }
            infile.close();   //关闭文件
        }
    }

    // 将控制点坐标转换为右手坐标系坐标
    leftHand2RightHand(objCoor);


    double lEOP[6] = { 800,0,-1000, 15.0 / 180 * PI,0,0 };        //使用右手坐标系！
    double rEOP[6] = { 5500 - 800,0,-1300,-15.0 / 180 * PI,0,0 }; //使用右手坐标系！
    double IOP[3] = { 20 / dpix, w / 2,h / 2 };

    //输出结果
    std::ofstream outfile;   //输出流

    //======== 左像片 ========
    string str = "left2";
    outfile.open("output/" + str + ".res", std::ios::trunc);
    if (!outfile.is_open())
    {
        std::cout << "output failed" << std::endl;
        return -1;
    }

    Resection leftImg;
    leftImg.setEOP(lEOP, 6);
    leftImg.setIOP(IOP, 6);
    leftImg.setCoor(leftImgCoor, objCoor);

    outfile << "======== 左像片 ========" << endl;
    leftImg.calculate(str);
    double lXsYxZx[3] = { 0 };

    lXsYxZx[0] = -leftImg.mEOP[2];
    lXsYxZx[1] = leftImg.mEOP[0];
    lXsYxZx[2] = leftImg.mEOP[1];

    outfile << "---- 外方位线元素(mm) ----" << endl;
    outfile << "Xs "<<lXsYxZx[0] << endl;
    outfile << "Ys " << lXsYxZx[1] << endl;
    outfile << "Zs " << lXsYxZx[2] << endl;

    outfile << "---- 外方位角元素(rad) ----" << endl;
    outfile << "phi " << leftImg.mEOP[3] << endl;
    outfile << "omega " << leftImg.mEOP[4] << endl;
    outfile << "kappa " << leftImg.mEOP[5] << endl;

    outfile << "---- 内方位元素(pix) ----" << endl;
    outfile << "f " << leftImg.mIOP[0] << endl;
    outfile << "x0 " << leftImg.mIOP[1] << endl;
    outfile << "y0 " << h - 1 - leftImg.mIOP[2] << endl;

    outfile << "---- 畸变系数 ----" << endl;
    outfile << "k1(pix^-2) " << leftImg.mParamK[0] << endl;
    outfile << "k2(pix^-4) " << leftImg.mParamK[1] << endl;
    outfile << "p1(pix^-1) " << leftImg.mParamP[0] << endl;
    outfile << "p2(pix^-1) " << leftImg.mParamP[1] << endl;
    outfile << endl;

    outfile.close();

    //======== 右像片 ========
    str = "right2";
    outfile.open("output/" + str + ".res", std::ios::trunc);
    if (!outfile.is_open())
    {
        std::cout << "output failed" << std::endl;
        return -1;
    }

    Resection rightImg;
    rightImg.setEOP(rEOP, 6);
    rightImg.setIOP(IOP, 6);
    rightImg.setCoor(rightImgCoor, objCoor);

    outfile << "======== 右像片 ========" << endl;
    rightImg.calculate(str);
    double rXsYxZx[3] = { 0 };

    rXsYxZx[0] = -rightImg.mEOP[2];
    rXsYxZx[1] = rightImg.mEOP[0];
    rXsYxZx[2] = rightImg.mEOP[1];

    outfile << "---- 外方位线元素(mm) ----" << endl;
    outfile << "Xs " << rXsYxZx[0] << endl;
    outfile << "Ys " << rXsYxZx[1] << endl;
    outfile << "Zs " << rXsYxZx[2] << endl;

    outfile << "---- 外方位角元素(rad) ----" << endl;
    outfile << "phi " << rightImg.mEOP[3] << endl;
    outfile << "omega " << rightImg.mEOP[4] << endl;
    outfile << "kappa " << rightImg.mEOP[5] << endl;

    outfile << "---- 内方位元素(pix) ----" << endl;
    outfile << "f " << rightImg.mIOP[0] << endl;
    outfile << "x0 " << rightImg.mIOP[1] << endl;
    outfile << "y0 " << h - 1 - rightImg.mIOP[2] << endl;

    outfile << "---- 畸变系数 ----" << endl;
    outfile << "k1(pix^-2) " << rightImg.mParamK[0] << endl;
    outfile << "k2(pix^-4) " << rightImg.mParamK[1] << endl;
    outfile << "p1(pix^-1) " << rightImg.mParamP[0] << endl;
    outfile << "p2(pix^-1) " << rightImg.mParamP[1] << endl;
    outfile << endl;

    outfile.close();

    return 0;
}

/*
 *   Z                                          Y
 *   |    X                                     |
 *   |   /                                      |
 *   |  /                     =====>            |
 *   | /                                        |----------- X
 *   |/                                        /
 *   |----------------- Y                     /
 *                                          Z/
 */
void leftHand2RightHand(std::map<int, Point3f>& objCoor)
{
    double x, y, z;
    for (auto it = objCoor.begin(); it != objCoor.end(); it++)
    {
        x = it->second.mx;
        y = it->second.my;
        z = it->second.mz;

        it->second.mx = y;
        it->second.my = z;
        it->second.mz = -x;
    }
}

/*
 *        Y                                     Z
 *        |                                     |    X 
 *        |                                     |   / 
 *        |                 =====>              |  /   
 *        |----------- X                        | /                           
 *       /                                      |/            
 *      /                                       |----------------- Y
 *    Z/
 */
void rightHand2LeftHand(std::map<int, Point3f>& objCoor)
{
    double x, y, z;
    for (auto it = objCoor.begin(); it != objCoor.end(); it++)
    {
        x = it->second.mx;
        y = it->second.my;
        z = it->second.mz;

        it->second.mx = -z;
        it->second.my = x;
        it->second.mz = y;
    }
}