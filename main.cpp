/*
#include<stdio.h>
#include<stdlib.h>

#include <windows.h>
#include <GL/glut.h>
#include <cmath>
#include <iostream>
#include <cstdlib>   // rand and srand
#include <ctime>
#include <vector>
*/

#include <windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

// #include <stdlib.h>
#include "1605066_classes.h"
#include "bitmap_image.hpp"
// #include <vector>
// #include <bits/stdc++.h>
// #include <sstream>

using namespace std;

/// ==================variables=============================
// vector<Object *> objects;
// vector<Light> lights;
// int lightSourceNo, levelOfRecursion, pixels, objectNo;
// int imageWidth, imageHeight;
const double fovY = 90.0;

const int axisLength = 100;
const int windowHeight = 500;
const int windowWidth = 500;
const float scaling = 1.3;
const float angleScaling = 0.08;

float cameraAngle = 1.0;
int drawaxes;
float value = (1 / sqrt(2));

const float fullBodyCutOffAngle = 45.0;
float fullBodyAngle = 0.0;
float gunRotateAngle = 0;
float gunUpDownAngle = 0;
float HalfBodyUpDownAngle = 0;

// Point3D pos = Point3D(150.0, 150.0, 50.0); //eye position
Point3D u = Point3D(0, 0, 1);
Point3D r = Point3D(-value, value, 0);
Point3D l = Point3D(-value, -value, 0);
Point3D center(0, 0, 0);

//float moveX = 100, moveY = 100.0, moveZ = 0.0;

class AngleClass
{
public:
	float fullBodyAngle;
	float HalfBodyUpDownAngle;
	float gunUpDownAngle;
	float gunRotateAngle;

	AngleClass()
	{
		fullBodyAngle = 0.0;
		HalfBodyUpDownAngle = 0.0;
		gunUpDownAngle = 0.0;
		gunRotateAngle = 0.0;
	}
};

vector<AngleClass> saveCurrentAngles;

void drawAxes()
{
	if (drawaxes == 1)
	{
		glBegin(GL_LINES);
		{
			glColor3f(1.0, 0, 0); //x axis red
			glVertex3f(axisLength, 0, 0);
			glVertex3f(-axisLength, 0, 0);

			glColor3f(0, 0, 1); //y axis blue
			glVertex3f(0, axisLength, 0);
			glVertex3f(0, -axisLength, 0);

			glColor3f(1, 1, 1); //z axis white
			glVertex3f(0, 0, axisLength);
			glVertex3f(0, 0, -axisLength);
		}
		glEnd();
	}
}

void getNexPoistionVector(Point3D *v, int option)
{

	if (option == 1)
	{
		pos.x = pos.x + scaling * v->x;
		pos.y = pos.y + scaling * v->y;
		pos.z = pos.z + scaling * v->z;
	}
	else if (option == 2)
	{
		pos.x = pos.x - scaling * v->x;
		pos.y = pos.y - scaling * v->y;
		pos.z = pos.z - scaling * v->z;
	}
}

void rotateVector(Point3D *L, Point3D *U, float angle)
{
	Point3D L_prime = Point3D();

	L_prime.x = L->x * cos(pi * angle / 180) + U->x * sin(pi * angle / 180);
	L_prime.y = L->y * cos(pi * angle / 180) + U->y * sin(pi * angle / 180);
	L_prime.z = L->z * cos(pi * angle / 180) + U->z * sin(pi * angle / 180);

	U->x = U->x * cos(pi * angle / 180) - L->x * sin(pi * angle / 180);
	U->y = U->y * cos(pi * angle / 180) - L->y * sin(pi * angle / 180);
	U->z = U->z * cos(pi * angle / 180) - L->z * sin(pi * angle / 180);

	L->x = L_prime.x;
	L->y = L_prime.y;
	L->z = L_prime.z;
	//cout<<"----------------x: "<<L->x<<" y: "<<L->y<< " z: "<<L->z<<endl;
	//cout<<"------DOT PRODUCT--: "<<l.getDotProduct(u)<<endl;
}

void capture()
{
	cout << "Capturing Method called" << endl;
	bitmap_image image(imageWidth, imageHeight);
	double planeDistance = (windowHeight / 2.0) / tan((cameraAngle * .5));
	Point3D topLeft;

	Point3D L = scaleVector(l, planeDistance);
	Point3D R = scaleVector(r, windowWidth / 2);
	Point3D U = scaleVector(u, windowHeight / 2);

	topLeft = substractTwoPointVector(addTwoPointVector(addTwoPointVector(pos, L), U), R);
	double du = (1.0 * windowWidth) / imageWidth;	//dx
	double dv = (1.0 * windowHeight) / imageHeight; //dy

	topLeft = addTwoPointVector(topLeft, substractTwoPointVector(scaleVector(r, 0.5 * du), scaleVector(u, 0.5 * dv)));
	Point3D currentMiddleCornerPoint;
	Color dummyColor = Color(0.1, 0.1, 0.1);

	int nearestObjectIndex = -1;
	bool isIntersect = false;

	for (int i = 0; i < pixels; i++)
	{
		for (int j = 0; j < pixels; j++)
		{
			Color intersectedObjectColor;
			isIntersect = false;
			nearestObjectIndex = -1;
			double t, tMin = DBL_MAX;
			currentMiddleCornerPoint = addTwoPointVector(topLeft, substractTwoPointVector(scaleVector(r, i * du), scaleVector(u, j * dv)));
			Point3D Rdir = substractTwoPointVector(currentMiddleCornerPoint, pos);

			Rdir.normalizePointVector();

			/// create ray from the eye to view point
			//-------add RayObject---
			Ray ray = Ray(pos, Rdir);

			for (int k = 0; k < objects.size(); k++)
			{
				//cout << "i " << i << "j: " << j << "k: " << k << endl;

				t = objects[k]->intersect(ray, &dummyColor, 0, k);

				if (t > 0 && t < tMin)
				{
					tMin = t;
					nearestObjectIndex = k;
					isIntersect = true;
				}
			}
			if (isIntersect)
			{
				intersectedObjectColor = objects[nearestObjectIndex]->color;
				Color objectColor = Color(objects[nearestObjectIndex]->color.r, objects[nearestObjectIndex]->color.g, objects[nearestObjectIndex]->color.b);

				objects[nearestObjectIndex]->intersect(ray, &objectColor, 1, nearestObjectIndex);
				image.set_pixel(i, j, objectColor.r * 255, objectColor.g * 255, objectColor.b * 255);

				//  Color color;
			}
		}
	}

	//image.save_image("D:\\BUET\\L-4 T-1\\CSE 410 (Graphics)\\Offline-3\\1605066_1\\1605066_out.bmp");
	image.save_image("D:\\BUET\\L-4 T-1\\CSE 410 (Graphics)\\Offline-3\\1605066_1\\1605066_out.bmp");

	cout << "Image captured" << endl;
}

void keyboardListener(unsigned char key, int x, int y)
{
	switch (key)
	{
	case '0':
		//capture() will be call here
		capture();
		break;
	case '1':
		///rotate/LookLeft--->U fixed
		//cameraAngle +=angleScaling;
		rotateVector(&r, &l, cameraAngle);
		cout << "---rotate/LookLeft\n---DOT PRODUCT--: " << l.getDotProduct(r) << endl;
		break;
	case '2':
		///rotate/LookRight --->U fixed
		//cameraAngle -=angleScaling;
		rotateVector(&r, &l, -cameraAngle);
		cout << "---rotate/LookRight\n---DOT PRODUCT--: " << l.getDotProduct(r) << endl;
		break;
	case '3':
		///--->R fixed
		//cameraAngle +=angleScaling;
		rotateVector(&l, &u, cameraAngle);
		break;
	case '4':
		///--->R fixed
		rotateVector(&l, &u, -cameraAngle);
		break;
	case '5':
		///tilt Clockwise--->L fixed
		
		rotateVector(&u, &r, cameraAngle);
		cout<<"---tilt Clockwise\n---DOT PRODUCT--: "<<u.getDotProduct(r)<<endl;
		break;
	case '6':
		///tilt AntiClockwise--->L fixed
		rotateVector(&u, &r, -cameraAngle);
		//cout<<"---tilt AntiClockwise\n---DOT PRODUCT--: "<<u.getDotProduct(r)<<endl;
		break;

	case 'q':

		break;
	case 'w':

		break;
	case 'e':

		break;
	case 'r':

		break;
	case 'a':

		break;
	case 's':

		break;
	case 'd':

		break;
	case 'f':

		break;

	default:
		break;
	}
}

void specialKeyListener(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_DOWN: //down arrow key
		//cameraHeight -= 3.0;
		//pos.y += 2;
		getNexPoistionVector(&l, 2);
		break;
	case GLUT_KEY_UP: // up arrow key
		//cameraHeight += 3.0;
		//pos.y -= 2;
		getNexPoistionVector(&l, 1);
		break;

	case GLUT_KEY_RIGHT:
		//cameraAngle += 0.03;
		//pos.x -= 2;
		getNexPoistionVector(&r, 1);
		break;
	case GLUT_KEY_LEFT:
		//cameraAngle -= 0.03;
		//pos.x += 2;
		getNexPoistionVector(&r, 2);
		break;

	case GLUT_KEY_PAGE_UP:
		//pos.z += 2;
		getNexPoistionVector(&u, 1);
		break;
	case GLUT_KEY_PAGE_DOWN:
		//pos.z -= 2;
		getNexPoistionVector(&u, 2);
		break;

	case GLUT_KEY_INSERT:
		break;

	case GLUT_KEY_HOME:
		break;
	case GLUT_KEY_END:
		break;

	default:
		break;
	}
}

void mouseListener(int button, int state, int x, int y)
{ //x, y is the x-y of the screen (2D)
	switch (button)
	{
	case GLUT_LEFT_BUTTON:
		if (state == GLUT_DOWN)
		{
			/*
                cout<<"Current Angle----\nfullBodyAngle:"<<fullBodyAngle<<"\tHalfBodyUpDownAngle:"
                <<HalfBodyUpDownAngle<<"\tgunUpDownAngle:"<<gunUpDownAngle<<"\tgunRotateAngle:"<<gunRotateAngle<<endl;
                */
			AngleClass currentAngle;
			currentAngle.fullBodyAngle = fullBodyAngle;
			currentAngle.HalfBodyUpDownAngle = HalfBodyUpDownAngle;
			currentAngle.gunUpDownAngle = gunUpDownAngle;
			currentAngle.gunRotateAngle = gunRotateAngle;

			saveCurrentAngles.push_back(currentAngle);
			//getAllCurrentAngles();
		}
		break;

	case GLUT_RIGHT_BUTTON:
		if (state == GLUT_DOWN)
		{ // 2 times?? in ONE click? -- solution is checking DOWN or UP
			drawaxes = 1 - drawaxes;
		}
		break;

	case GLUT_MIDDLE_BUTTON:
		//........
		break;

	default:
		break;
	}
}

void display()
{

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0, 0, 0, 0); //color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);

	drawAxes();
	// floor.draw();
	for (int i = 0; i < objects.size(); i++)
	{
		/* code */
		objects[i]->draw();
	}

	glutSwapBuffers();
}

void animate()
{
	///codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init()
{

	drawaxes = 1;
	glClearColor(0, 0, 0, 0);

	/************************
	/ set-up projection here
	************************/
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	///give PERSPECTIVE parameters
	gluPerspective(80, 1, 1, 1000.0);
	//glFrustum(-1,1,-1,1,2,10);
	glMatrixMode(GL_MODELVIEW);
}

void loadData()
{
	// string myInputText;
	ifstream MyReadFile;
	MyReadFile.open("D:\\BUET\\L-4 T-1\\CSE 410 (Graphics)\\Offline-3\\1605066_1\\scene.txt");
	MyReadFile >> levelOfRecursion;
	MyReadFile >> pixels; //imageWidth = imageHeight
	imageWidth = pixels;
	imageHeight = pixels;
	MyReadFile >> objectNo;

	double a, d, s, r; //co-efficients
	int shine;
	Object *temp;

	for (int i = 0; i < objectNo; i++)
	{
		string objectType;
		MyReadFile >> objectType;
		cout << i + 1 << " : " << objectType << endl;
		if (objectType == "sphere")
		{
			Point3D center;
			Color color;
			double radius;
			MyReadFile >> center.x >> center.y >> center.z;
			MyReadFile >> radius;
			MyReadFile >> color.r >> color.g >> color.b;
			temp = new Sphere(center, radius, color);
			MyReadFile >> a >> d >> s >> r;
			temp->setCoEfficients(a, d, s, r);
			MyReadFile >> shine;
			temp->setShine(shine);

			objects.push_back(temp);
		}
		else if (objectType == "triangle")
		{
			// cout << "in triangle" << endl;
			Point3D A, B, C;
			Color color;
			MyReadFile >> A.x >> A.y >> A.z;
			MyReadFile >> B.x >> B.y >> B.z;
			MyReadFile >> C.x >> C.y >> C.z;

			MyReadFile >> color.r >> color.g >> color.b;

			temp = new Triangle(A, B, C, color);
			MyReadFile >> a >> d >> s >> r;
			temp->setCoEfficients(a, d, s, r);
			MyReadFile >> shine;
			temp->setShine(shine);

			objects.push_back(temp);
		}
		else if (objectType == "general")
		{
			double A, B, C, D, E, F, G, H, I, J;
			Point3D refPoint;

			// length, width, height (0 indicates no clipping along this dimension)
			double l, w, h; //length,width, height-->cutOff
			Color color;
			double tempValue;
			vector<double> quadCoefficients;

			for (int i = 0; i < 10; i++)
			{
				MyReadFile >> tempValue;
				quadCoefficients.push_back(tempValue);
			}

			// MyReadFile >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
			MyReadFile >> refPoint.x >> refPoint.y >> refPoint.z >> l >> w >> h;
			MyReadFile >> color.r >> color.g >> color.b;
			MyReadFile >> a >> d >> s >> r;
			MyReadFile >> shine;
			temp = new Quadric(quadCoefficients, refPoint, color, l, w, h);
			temp->setCoEfficients(a, d, s, r);
			temp->setShine(shine);

			objects.push_back(temp);
		}
	}

	MyReadFile >> lightSourceNo;
	for (int i = 0; i < lightSourceNo; i++)
	{
		Light light;
		MyReadFile >> light.lightPosition.x >> light.lightPosition.y >> light.lightPosition.z;
		MyReadFile >> light.color.r >> light.color.g >> light.color.b;
		lights.push_back(light);
	}

	temp = new Floor(1000, 20);
	temp->setCoEfficients(0.5, 0.4, 0.3, 0.1);
	temp->setShine(5);
	objects.push_back(temp);

	// cout << "TEXT:" << myInputText << endl;

	///----------------Test Code -----------
	/*
	Object *tmp;
	//tmp = new Floor(1000,20);
	//objects.push_back(tmp);

	tmp = new Sphere(Point3D(10.0, 50.0, 5.0),5.0,Color(0.30, 1.0, 0.0));
	// tmp = new Sphere(10.0, 50.0, 5.0,5.0,Color(0.30, 1.0, 0.0));
	// tmp = new Triangle(Point3D(-20.0, -20.0, 0.0), Point3D(20.0, -20.0, 0.0), Point3D(0.0, 0.0, 20.0), Color(1.0, 0.0, 0.0));
	tmp->setCoEfficients(0.4,0.2,0.1,0.3);
	tmp->setShine(5);
	objects.push_back(tmp);

	*/

	MyReadFile.close();
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitWindowSize(windowHeight, windowWidth);
	glutInitWindowPosition(600, 150);

	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Camera-1");

	init();
	loadData();

	glEnable(GL_DEPTH_TEST);  //enable Depth Testing
	glutDisplayFunc(display); //display callback function
	glutIdleFunc(animate);	  //what you want to do in the idle time (when no drawing is occuring)
	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();
	return 0;
}
