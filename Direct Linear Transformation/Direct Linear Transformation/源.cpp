#include <stdio.h>
#include <istream>
#include <fstream> 
#include <iostream>


#include "Eigen/Eigen"  

#include "point.h"
#include "DLT.h"
#define PI 3.14159265358979

using namespace std;
using namespace Eigen;


void leftHand2RightHand(std::map<int, Point3f>& ptMap);
void rightHand2LeftHand(std::map<int, Point3f>& ptMap);

int main()
{
    /*  直接线性变换:
    x + vx + dx + ( L1 * X + L2 * Y + L3 * Z + L4)/( L9 * X + L10 * Y + L11 * Z + 1 ) = 0
    y + vy + dy + ( L5 * X + L6 * Y + L7 * Z + L8)/( L9 * X + L10 * Y + L11 * Z + 1 ) = 0 

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
    string checkConfigPath = "data/checkpoint.cig";
    int w = 4272;
    int h = 2848;




    //左像片控制点、待定点、检查点的像点坐标(pix)
    std::map<int, Point2f> leftCtrImgCoor;
    std::map<int, Point2f> leftUknImgCoor;
    std::map<int, Point2f> leftChkImgCoor;
    //右像片控制点和待定点、检查点的像点坐标(pix)
    std::map<int, Point2f> rightCtrImgCoor;
    std::map<int, Point2f> rightUknImgCoor;
    std::map<int, Point2f> rightChkImgCoor;
    //控制点的物方坐标(mm)
    std::map<int, Point3f> ctrObjCoor;
    //检查点的物方坐标(mm)
    std::map<int, Point3f> chkObjCoor;
    //待定点的物方坐标(mm)
    std::map<int, Point3f> uknObjCoor;
    ifstream infile;   //输入流

    leftCtrImgCoor.clear();
    leftUknImgCoor.clear();
    leftChkImgCoor.clear();
    rightCtrImgCoor.clear();
    rightUknImgCoor.clear();
    rightChkImgCoor.clear();
    ctrObjCoor.clear();
    chkObjCoor.clear();
    uknObjCoor.clear();

    //读取左影像控制点、待定点像点坐标
    {
        int ID = 0;
        double x = 0, y = 0, height = 0, width = 0, angle = 0;

        //先读控制点和待定点
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
                    leftCtrImgCoor[ID] = Point2f(x, h - 1 - y);
                }
                else
                {
                    leftUknImgCoor[ID] = Point2f(x, h - 1 - y);
                }
            }
            infile.close();   //关闭文件
        }
    }

    //读取右影像控制点、待定点像点坐标
    {
        int ID = 0;
        double x = 0, y = 0, height = 0, width = 0, angle = 0;

        //先读控制点和待定点
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
                    rightCtrImgCoor[ID] = Point2f(x, h - 1 - y);
                }
                else
                {
                    rightUknImgCoor[ID] = Point2f(x, h - 1 - y);
                }
            }
            infile.close();   //关闭文件
        }
    }

    //再选检查点
    {
        int ID = 0;
        infile.open(checkConfigPath, ios::in);
        if (!infile.is_open())
        {
            cout << "check point config file not found." << endl;
        }
        else {
            char temp[128];
            int num = 0;
            infile >> num;
            for (int i = 0; i < num; i++)
            {
                infile >> ID;
                if (ID >= 100)
                {
                    auto it1 = leftCtrImgCoor.find(ID);
                    auto it2 = rightCtrImgCoor.find(ID);
                    if (it1 != leftCtrImgCoor.end() && it2 != rightCtrImgCoor.end())
                    {
                        leftChkImgCoor[ID] = it1->second;
                        rightChkImgCoor[ID] = it2->second;
                        //移出控制点
                        leftCtrImgCoor.erase(it1);
                        rightCtrImgCoor.erase(it2);
                    }
                }
            }
            infile.close();   //关闭文件
        }
    }

    //读取控制点、检查点物方坐标
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
                ctrObjCoor[ID] = Point3f(x, y, z);
            }
            infile.close();   //关闭文件
        }
        //选检查点
        for (auto it = leftChkImgCoor.begin(); it != leftChkImgCoor.end(); it++)
        {
            auto it1 = ctrObjCoor.find(it->first);
            chkObjCoor[it->first] = it1->second;
            ctrObjCoor.erase(it1);
        }

    }

    // 将控制点坐标转换为右手坐标系坐标
    leftHand2RightHand(ctrObjCoor);
    leftHand2RightHand(chkObjCoor);

    DLT dlt;
    string str = "DLT2";

    dlt.setCtrImgCoor(leftCtrImgCoor, rightCtrImgCoor);
    dlt.setUknImgCoor(leftUknImgCoor, rightUknImgCoor);
    dlt.setChkImgCoor(leftChkImgCoor, rightChkImgCoor);
    dlt.setCtrObjCoor(ctrObjCoor);
    dlt.setChkObjCoor(chkObjCoor);
    dlt.setStr(str);

    dlt.calculateLValue();
    dlt.calculateUknObjCoor();

    //======== 计算结果 ========
    std::ofstream outfile;   //输出流
    outfile.open("output/" + str + ".res", std::ios::trunc);
    if (!outfile.is_open())
    {
        std::cout << "output failed" << std::endl;
        return -1;
    }

    std::map<int, Point3f> temp = dlt.mUknObjCoor;
    rightHand2LeftHand(temp);
    auto it52 = temp.find(52);

    outfile << "======== 计算待定点坐标X,Y,Z(mm)和与52号点间的距离(mm) ========" << std::endl;
    for (auto it = temp.begin(); it != temp.end(); it++)
    {
        outfile << it->first << " " << it->second.mx << " " << it->second.my << " " << it->second.mz << " " << sqrt(pow(it->second.mx - it52->second.mx, 2)
            + pow(it->second.my - it52->second.my, 2) + pow(it->second.mz - it52->second.mz, 2)) << std::endl;
    }
    outfile.close();


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