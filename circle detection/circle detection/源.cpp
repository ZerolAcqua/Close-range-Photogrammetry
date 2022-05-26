#include <stdio.h>
#include <istream>
#include <fstream> 
#include <iostream>

#include <map> 

#include <opencv2/core/core.hpp>   
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgproc/types_c.h>


using namespace std;
using namespace cv;



Mat rawImg;			//原彩色图,不叠加其他东西的，用于各种图像的处理
Mat src4showImg;	//原彩色图,但会叠加其他东西的图像
Mat scalesImg;		//缩放后的彩色图
Mat grayImg;		//灰度图
Mat preImg;			//预览用图
Mat focusImg;		//框选区域用图,不叠加其他东西
Mat focus4showImg;	//框选区域用图,会叠加其他东西


double scale = 1;	//缩放比例

Rect select;				//框选的矩形
bool select_flag = false;	//框选状态
Point origin;				//框选的起点

#define SHOW_RESULT 0



void onMouse(int event, int x, int y, int, void*)
{
	if (event == EVENT_LBUTTONDOWN)
	{
		select_flag = true;			//左键按下赋真值
		origin = Point(x, y);		//保存左键单击的坐标
		select = Rect(x, y, 0, 0);

	}
	else if (event == EVENT_LBUTTONUP)
	{
		select_flag = false;
	}
	if (select_flag)
	{
		select.x = MIN(origin.x, x);					//鼠标按下开始到弹起这段时间实时计算所选矩形框左上角点坐标
		select.y = MIN(origin.y, y);
		select.width = abs(x - origin.x);				//算矩形宽度和高度
		select.height = abs(y - origin.y);
		select &= Rect(0, 0, rawImg.cols * scale, rawImg.rows * scale);	//保证所选矩形框在图片区域之内
	}
}

int main()
{
	scale = 1.0 / 4;

	int ID = 0;			//编号ID
	int maxArea = 5000;
	int minArea = 200;
	double roundness = 0.85;
	double distance = 5;
	int limitSize = 600;

	//点位点号数据
	map<int, RotatedRect> posMap;

	RotatedRect boundBox;
	string str1 = "data/IMG_8621.JPG";

	cout << "回车确定选框，Z取消选框，ESC退出，D删点" << endl;

	//cout << "文件路径:";
	//cin >> str1;
	//cout << "缩放比例(0.25):";
	//cin >> scale;
	//cout << "最小圆面积(150):";
	//cin >> minArea;
	//cout << "最大圆面积(1700):";
	//cin >> maxArea;
	//cout << "最小圆度(0.8):";
	//cin >> roundness;
	//cout << "最小圆心间距(10):";
	//cin >> distance;

	//读取点位数据
	double x, y, height, width, angle;
	ifstream infile;   //输入流


	infile.open(str1 + ".dat", ios::in);
	if (!infile.is_open())
	{
		cout << "position file not found." << endl;
	}
	else {
		char temp[128];
		infile.getline(temp, 128);		//跳过注释信息
		while (!infile.eof())           // 若未到文件结束一直循环
		{
			infile >> ID >> x >> y >> width >> height >> angle;
			posMap[ID] = RotatedRect(Point2f(x, y), Size2f(width, height), angle);
		}
		infile.close();   //关闭文件
	}



	rawImg = imread(str1);
	rawImg.copyTo(src4showImg);
	focusImg = rawImg;
	focus4showImg = src4showImg;

	resize(rawImg, preImg, Size(), scale, scale);
	imshow("原图", preImg);


	//捕捉鼠标事件
	setMouseCallback("原图", onMouse, 0);
	char c;		//用于退出

	while (true)
	{

		rawImg.copyTo(src4showImg);

		//设置绘制文本和圆
		for (auto it = posMap.begin(); it != posMap.end(); it++)
		{
			//绘制圆
			//椭圆拟合
			//画出拟合的椭圆
			ellipse(src4showImg, it->second, Scalar(0, 0, 255), 5, FILLED);


			string text = to_string(it->first);
			int font_face = cv::FONT_HERSHEY_COMPLEX;
			double font_scale = 0.5 / scale;
			int thickness = 0.5 / scale;
			int baseline;
			//获取文本框的长宽
			cv::Size text_size = cv::getTextSize(text, font_face, font_scale, thickness, &baseline);
			//将文本框居中绘制
			cv::putText(src4showImg, text, it->second.center, font_face, font_scale, cv::Scalar(255, 0, 255), thickness, 8, 0);
		}



		resize(src4showImg, preImg, Size(), scale, scale);


		//画出矩形框
		rectangle(preImg, select, Scalar(255, 0, 0), 1, 8, 0);//显示在画矩形窗口时的痕迹

		//显示图片到窗口
		imshow("原图", preImg);

		if (select.width > 0 && select.height > 0)
		{
			focusImg = rawImg(Rect(select.x / scale, select.y / scale, select.width / scale, select.height / scale));
			focus4showImg = src4showImg(Rect(select.x / scale, select.y / scale, select.width / scale, select.height / scale));
		}

		//键盘响应
		c = (char)waitKey(20);

		if (c == 13 && (select.width > 0 && select.height > 0))//回车键选定
		{
			//暂时解除鼠标事件
			setMouseCallback("原图", NULL, NULL);



			/////////////////////以下是处理提取的核心部分/////////////////////
			cvtColor(focusImg, grayImg, CV_BGR2GRAY);

			//高斯模糊
			GaussianBlur(grayImg, grayImg, Size(3, 3), 0);

			////直方图均衡化
			//equalizeHist(grayImg, grayImg);

			//边缘提取
			Mat edgeImg;
			Canny(grayImg, edgeImg, 30, 70);

			//追踪轮廓
			vector<vector<Point>> contours;
			findContours(edgeImg, contours, CV_RETR_LIST, CV_CHAIN_APPROX_NONE);

			double area = 0;
			double perimeter = 0;
			Mat cirImg = Mat::zeros(edgeImg.size(), CV_8UC3);
			Point center = Point(-999, -999);

			//找出并筛选"圆形"
			for (auto it1 = contours.begin(); it1 != contours.end(); )
			{
				//拟合点个数筛选
				if (it1->size() < 6)
				{
					it1 = contours.erase(it1);
					continue;
				}
				//面积筛选与圆度筛选
				//这里面积的最小阈值非常关键
				area = contourArea(*it1);
				perimeter = arcLength(*it1, true);
				if (area <minArea || area>maxArea)
				{
					it1 = contours.erase(it1);
					continue;
				}
				if (4 * 3.1415926 * area / pow(perimeter, 2) < roundness)
				{
					//圆形度不符的被剔除
					it1 = contours.erase(it1);
					continue;
				}

				//简单排除过于接近的圆，主要是同一个标识的内外轮廓有重复检测
				Point temp = fitEllipse(*it1).center;
				if (sqrtf(powf((center.x - temp.x), 2) + powf((center.y - temp.y), 2)) < distance)
				{
					it1 = contours.erase(it1);
					continue;
				}
				else
				{
					center = temp;
				}

				//迭代器后移
				it1++;
			}

			/////////////////////以上是处理提取的核心部分/////////////////////
			//预览点位，并调整窗口
			for (auto it = contours.begin(); it != contours.end(); it++)
			{
				//椭圆拟合
				boundBox = fitEllipse(*it);
				//画出拟合的椭圆
				ellipse(focus4showImg, boundBox, Scalar(0, 0, 255), 5, FILLED);
			}
			if (MAX(focus4showImg.cols, focus4showImg.rows) > limitSize)
			{
				double tempScale = 1.0 * limitSize / MAX(focus4showImg.cols, focus4showImg.rows);
				resize(focus4showImg, scalesImg, Size(), tempScale, tempScale);
				imshow("拟合结果", scalesImg);
			}
			else
			{
				imshow("拟合结果", focus4showImg);
			}
			while ((c = waitKey(20)) != 13 && c != 'z')
			{

			}
			if (c == 'z')
			{
				//恢复鼠标事件
				setMouseCallback("原图", onMouse, 0);
				continue;
			}


			//标记圆形
			auto it1 = contours.begin();
			for (auto it1 = contours.begin(); it1 != contours.end(); it1++)
			{
				for (auto it2 = contours.begin(); it2 != contours.end(); it2++)
				{
					if (it2 == it1)
						continue;
					//椭圆拟合
					boundBox = fitEllipse(*it2);
					//画出拟合的椭圆
					ellipse(focus4showImg, boundBox, Scalar(0, 0, 255), 5, FILLED);
				}
				//标记当前要标号的圆形
				//椭圆拟合
				boundBox = fitEllipse(*it1);
				//画出拟合的椭圆
				ellipse(focus4showImg, boundBox, Scalar(255, 0, 0), 5, FILLED);

				if (MAX(focus4showImg.cols, focus4showImg.rows) > limitSize)
				{
					double tempScale = 1.0 * limitSize / MAX(focus4showImg.cols, focus4showImg.rows);
					resize(focus4showImg, scalesImg, Size(), tempScale, tempScale);
					imshow("拟合结果", scalesImg);
					waitKey(30);
				}
				else
				{
					imshow("拟合结果", focus4showImg);
					waitKey(30);
				}


				cout << "输入蓝色圆的号码：";
				cin >> ID;

				RotatedRect temp = boundBox;
				temp.center.x += select.x / scale;
				temp.center.y += select.y / scale;
				posMap[ID] = temp;
			}

			for (auto it1 = contours.begin(); it1 != contours.end(); it1++)
			{
				//椭圆拟合
				boundBox = fitEllipse(*it1);
				//画出拟合的椭圆
				ellipse(focus4showImg, boundBox, Scalar(0, 0, 255), 5, FILLED);
			}


			if (MAX(focus4showImg.cols, focus4showImg.rows) > limitSize)
			{
				double tempScale = 1.0 * limitSize / MAX(focus4showImg.cols, focus4showImg.rows);
				resize(focus4showImg, scalesImg, Size(), tempScale, tempScale);
				imshow("拟合结果", scalesImg);
				waitKey(30);
			}
			else
			{
				imshow("拟合结果", focus4showImg);
				waitKey(30);
			}




			//恢复鼠标事件
			setMouseCallback("原图", onMouse, 0);
		}
		else if (c == 'd')
		{
			//d键删除点
			cout << "输入要删除的号码：";
			cin >> ID;
			auto it = posMap.find(ID);
			if (it != posMap.end())
			{
				posMap.erase(it);
			}
		}
		else if (c == 's')
		{
			imwrite("拟合结果.jpg", src4showImg);
			ofstream outfile;   //输出流
			outfile.open(str1 + ".dat", ios::trunc);
			if (!outfile.is_open())
				cout << "output failed" << endl;
			outfile << "# ID""\t" << "X""\t" << "Y""\t" << "width""\t" << "height""\t" << "angle""\t";
			for (auto it = posMap.begin(); it != posMap.end(); it++)
			{
				//在result.txt中写入结果
				outfile << endl
					<< it->first << "\t"
					<< it->second.center.x << "\t"
					<< it->second.center.y << "\t"
					<< it->second.size.width << "\t"
					<< it->second.size.height << "\t"
					<< it->second.angle;
			}
			outfile.close();
			system("cls");
			cout << "回车确定选框，Z取消选框，ESC退出，D删点" << endl;
		}
		else if (c == 27)
		{
			//ESC键退出
			break;
		}
	}
	imwrite("拟合结果.jpg", src4showImg);
	ofstream outfile;   //输出流
	outfile.open(str1 + ".dat", ios::trunc);
	if (!outfile.is_open())
		cout << "output failed" << endl;
	outfile << "# ID""\t" << "X""\t" << "Y""\t" << "width""\t" << "height""\t" << "angle""\t";
	for (auto it = posMap.begin(); it != posMap.end(); it++)
	{
		//在result.txt中写入结果
		outfile << endl
			<< it->first << "\t"
			<< it->second.center.x << "\t"
			<< it->second.center.y << "\t"
			<< it->second.size.width << "\t"
			<< it->second.size.height << "\t"
			<< it->second.angle;
	}
	outfile.close();

}




