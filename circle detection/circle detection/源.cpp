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



Mat rawImg;			//ԭ��ɫͼ,���������������ģ����ڸ���ͼ��Ĵ���
Mat src4showImg;	//ԭ��ɫͼ,�����������������ͼ��
Mat scalesImg;		//���ź�Ĳ�ɫͼ
Mat grayImg;		//�Ҷ�ͼ
Mat preImg;			//Ԥ����ͼ
Mat focusImg;		//��ѡ������ͼ,��������������
Mat focus4showImg;	//��ѡ������ͼ,�������������


double scale = 1;	//���ű���

Rect select;				//��ѡ�ľ���
bool select_flag = false;	//��ѡ״̬
Point origin;				//��ѡ�����

#define SHOW_RESULT 0



void onMouse(int event, int x, int y, int, void*)
{
	if (event == EVENT_LBUTTONDOWN)
	{
		select_flag = true;			//������¸���ֵ
		origin = Point(x, y);		//�����������������
		select = Rect(x, y, 0, 0);

	}
	else if (event == EVENT_LBUTTONUP)
	{
		select_flag = false;
	}
	if (select_flag)
	{
		select.x = MIN(origin.x, x);					//��갴�¿�ʼ���������ʱ��ʵʱ������ѡ���ο����Ͻǵ�����
		select.y = MIN(origin.y, y);
		select.width = abs(x - origin.x);				//����ο�Ⱥ͸߶�
		select.height = abs(y - origin.y);
		select &= Rect(0, 0, rawImg.cols * scale, rawImg.rows * scale);	//��֤��ѡ���ο���ͼƬ����֮��
	}
}

int main()
{
	scale = 1.0 / 4;

	int ID = 0;			//���ID
	int maxArea = 5000;
	int minArea = 200;
	double roundness = 0.85;
	double distance = 5;
	int limitSize = 600;

	//��λ�������
	map<int, RotatedRect> posMap;

	RotatedRect boundBox;
	string str1 = "data/IMG_8621.JPG";

	cout << "�س�ȷ��ѡ��Zȡ��ѡ��ESC�˳���Dɾ��" << endl;

	//cout << "�ļ�·��:";
	//cin >> str1;
	//cout << "���ű���(0.25):";
	//cin >> scale;
	//cout << "��СԲ���(150):";
	//cin >> minArea;
	//cout << "���Բ���(1700):";
	//cin >> maxArea;
	//cout << "��СԲ��(0.8):";
	//cin >> roundness;
	//cout << "��СԲ�ļ��(10):";
	//cin >> distance;

	//��ȡ��λ����
	double x, y, height, width, angle;
	ifstream infile;   //������


	infile.open(str1 + ".dat", ios::in);
	if (!infile.is_open())
	{
		cout << "position file not found." << endl;
	}
	else {
		char temp[128];
		infile.getline(temp, 128);		//����ע����Ϣ
		while (!infile.eof())           // ��δ���ļ�����һֱѭ��
		{
			infile >> ID >> x >> y >> width >> height >> angle;
			posMap[ID] = RotatedRect(Point2f(x, y), Size2f(width, height), angle);
		}
		infile.close();   //�ر��ļ�
	}



	rawImg = imread(str1);
	rawImg.copyTo(src4showImg);
	focusImg = rawImg;
	focus4showImg = src4showImg;

	resize(rawImg, preImg, Size(), scale, scale);
	imshow("ԭͼ", preImg);


	//��׽����¼�
	setMouseCallback("ԭͼ", onMouse, 0);
	char c;		//�����˳�

	while (true)
	{

		rawImg.copyTo(src4showImg);

		//���û����ı���Բ
		for (auto it = posMap.begin(); it != posMap.end(); it++)
		{
			//����Բ
			//��Բ���
			//������ϵ���Բ
			ellipse(src4showImg, it->second, Scalar(0, 0, 255), 5, FILLED);


			string text = to_string(it->first);
			int font_face = cv::FONT_HERSHEY_COMPLEX;
			double font_scale = 0.5 / scale;
			int thickness = 0.5 / scale;
			int baseline;
			//��ȡ�ı���ĳ���
			cv::Size text_size = cv::getTextSize(text, font_face, font_scale, thickness, &baseline);
			//���ı�����л���
			cv::putText(src4showImg, text, it->second.center, font_face, font_scale, cv::Scalar(255, 0, 255), thickness, 8, 0);
		}



		resize(src4showImg, preImg, Size(), scale, scale);


		//�������ο�
		rectangle(preImg, select, Scalar(255, 0, 0), 1, 8, 0);//��ʾ�ڻ����δ���ʱ�ĺۼ�

		//��ʾͼƬ������
		imshow("ԭͼ", preImg);

		if (select.width > 0 && select.height > 0)
		{
			focusImg = rawImg(Rect(select.x / scale, select.y / scale, select.width / scale, select.height / scale));
			focus4showImg = src4showImg(Rect(select.x / scale, select.y / scale, select.width / scale, select.height / scale));
		}

		//������Ӧ
		c = (char)waitKey(20);

		if (c == 13 && (select.width > 0 && select.height > 0))//�س���ѡ��
		{
			//��ʱ�������¼�
			setMouseCallback("ԭͼ", NULL, NULL);



			/////////////////////�����Ǵ�����ȡ�ĺ��Ĳ���/////////////////////
			cvtColor(focusImg, grayImg, CV_BGR2GRAY);

			//��˹ģ��
			GaussianBlur(grayImg, grayImg, Size(3, 3), 0);

			////ֱ��ͼ���⻯
			//equalizeHist(grayImg, grayImg);

			//��Ե��ȡ
			Mat edgeImg;
			Canny(grayImg, edgeImg, 30, 70);

			//׷������
			vector<vector<Point>> contours;
			findContours(edgeImg, contours, CV_RETR_LIST, CV_CHAIN_APPROX_NONE);

			double area = 0;
			double perimeter = 0;
			Mat cirImg = Mat::zeros(edgeImg.size(), CV_8UC3);
			Point center = Point(-999, -999);

			//�ҳ���ɸѡ"Բ��"
			for (auto it1 = contours.begin(); it1 != contours.end(); )
			{
				//��ϵ����ɸѡ
				if (it1->size() < 6)
				{
					it1 = contours.erase(it1);
					continue;
				}
				//���ɸѡ��Բ��ɸѡ
				//�����������С��ֵ�ǳ��ؼ�
				area = contourArea(*it1);
				perimeter = arcLength(*it1, true);
				if (area <minArea || area>maxArea)
				{
					it1 = contours.erase(it1);
					continue;
				}
				if (4 * 3.1415926 * area / pow(perimeter, 2) < roundness)
				{
					//Բ�ζȲ����ı��޳�
					it1 = contours.erase(it1);
					continue;
				}

				//���ų����ڽӽ���Բ����Ҫ��ͬһ����ʶ�������������ظ����
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

				//����������
				it1++;
			}

			/////////////////////�����Ǵ�����ȡ�ĺ��Ĳ���/////////////////////
			//Ԥ����λ������������
			for (auto it = contours.begin(); it != contours.end(); it++)
			{
				//��Բ���
				boundBox = fitEllipse(*it);
				//������ϵ���Բ
				ellipse(focus4showImg, boundBox, Scalar(0, 0, 255), 5, FILLED);
			}
			if (MAX(focus4showImg.cols, focus4showImg.rows) > limitSize)
			{
				double tempScale = 1.0 * limitSize / MAX(focus4showImg.cols, focus4showImg.rows);
				resize(focus4showImg, scalesImg, Size(), tempScale, tempScale);
				imshow("��Ͻ��", scalesImg);
			}
			else
			{
				imshow("��Ͻ��", focus4showImg);
			}
			while ((c = waitKey(20)) != 13 && c != 'z')
			{

			}
			if (c == 'z')
			{
				//�ָ�����¼�
				setMouseCallback("ԭͼ", onMouse, 0);
				continue;
			}


			//���Բ��
			auto it1 = contours.begin();
			for (auto it1 = contours.begin(); it1 != contours.end(); it1++)
			{
				for (auto it2 = contours.begin(); it2 != contours.end(); it2++)
				{
					if (it2 == it1)
						continue;
					//��Բ���
					boundBox = fitEllipse(*it2);
					//������ϵ���Բ
					ellipse(focus4showImg, boundBox, Scalar(0, 0, 255), 5, FILLED);
				}
				//��ǵ�ǰҪ��ŵ�Բ��
				//��Բ���
				boundBox = fitEllipse(*it1);
				//������ϵ���Բ
				ellipse(focus4showImg, boundBox, Scalar(255, 0, 0), 5, FILLED);

				if (MAX(focus4showImg.cols, focus4showImg.rows) > limitSize)
				{
					double tempScale = 1.0 * limitSize / MAX(focus4showImg.cols, focus4showImg.rows);
					resize(focus4showImg, scalesImg, Size(), tempScale, tempScale);
					imshow("��Ͻ��", scalesImg);
					waitKey(30);
				}
				else
				{
					imshow("��Ͻ��", focus4showImg);
					waitKey(30);
				}


				cout << "������ɫԲ�ĺ��룺";
				cin >> ID;

				RotatedRect temp = boundBox;
				temp.center.x += select.x / scale;
				temp.center.y += select.y / scale;
				posMap[ID] = temp;
			}

			for (auto it1 = contours.begin(); it1 != contours.end(); it1++)
			{
				//��Բ���
				boundBox = fitEllipse(*it1);
				//������ϵ���Բ
				ellipse(focus4showImg, boundBox, Scalar(0, 0, 255), 5, FILLED);
			}


			if (MAX(focus4showImg.cols, focus4showImg.rows) > limitSize)
			{
				double tempScale = 1.0 * limitSize / MAX(focus4showImg.cols, focus4showImg.rows);
				resize(focus4showImg, scalesImg, Size(), tempScale, tempScale);
				imshow("��Ͻ��", scalesImg);
				waitKey(30);
			}
			else
			{
				imshow("��Ͻ��", focus4showImg);
				waitKey(30);
			}




			//�ָ�����¼�
			setMouseCallback("ԭͼ", onMouse, 0);
		}
		else if (c == 'd')
		{
			//d��ɾ����
			cout << "����Ҫɾ���ĺ��룺";
			cin >> ID;
			auto it = posMap.find(ID);
			if (it != posMap.end())
			{
				posMap.erase(it);
			}
		}
		else if (c == 's')
		{
			imwrite("��Ͻ��.jpg", src4showImg);
			ofstream outfile;   //�����
			outfile.open(str1 + ".dat", ios::trunc);
			if (!outfile.is_open())
				cout << "output failed" << endl;
			outfile << "# ID""\t" << "X""\t" << "Y""\t" << "width""\t" << "height""\t" << "angle""\t";
			for (auto it = posMap.begin(); it != posMap.end(); it++)
			{
				//��result.txt��д����
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
			cout << "�س�ȷ��ѡ��Zȡ��ѡ��ESC�˳���Dɾ��" << endl;
		}
		else if (c == 27)
		{
			//ESC���˳�
			break;
		}
	}
	imwrite("��Ͻ��.jpg", src4showImg);
	ofstream outfile;   //�����
	outfile.open(str1 + ".dat", ios::trunc);
	if (!outfile.is_open())
		cout << "output failed" << endl;
	outfile << "# ID""\t" << "X""\t" << "Y""\t" << "width""\t" << "height""\t" << "angle""\t";
	for (auto it = posMap.begin(); it != posMap.end(); it++)
	{
		//��result.txt��д����
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




