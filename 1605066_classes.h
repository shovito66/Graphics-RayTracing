#include <stdio.h>
#include <iostream>
#include <windows.h>
#include <GL/glut.h>
#include <bits/stdc++.h>

#define pi (2 * acos(0.0))
using namespace std;

///-----
#define AMBIENT 0
#define DIFFUSE 1
#define SPECULAR 2
#define REFLECTION 3
///------------------------

int stacks = 40;
int slices = 60;

class Point3D
{
public:
    double x, y, z;
    Point3D()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    Point3D(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    double getDotProduct(Point3D B)
    {
        return (x * B.x + y * B.y + z * B.z);
    }

    void normalizePointVector()
    {
        double magnitude = sqrt(pow(this->x, 2) + pow(this->y, 2) + pow(this->z, 2));
        this->x = this->x / magnitude;
        this->y = this->y / magnitude;
        this->z = this->z / magnitude;
    }

    void printPoint()
    {
        cout << "printing the POINT:" << this->x << " " << this->y << " " << this->z << endl;
    }
};

Point3D pos = Point3D(150.0, 150.0, 50.0); //eye position

Point3D scaleVector(Point3D A, double scaleFactor)
{
    return Point3D(A.x * scaleFactor, A.y * scaleFactor, A.z * scaleFactor);
}

Point3D addTwoPointVector(Point3D A, Point3D B)
{
    return Point3D(A.x + B.x, A.y + B.y, A.z + B.z);
}

Point3D substractTwoPointVector(Point3D A, Point3D B)
{
    return Point3D(A.x - B.x, A.y - B.y, A.z - B.z);
}

Point3D crossProductVector(Point3D a, Point3D b)
{
    ///cross: (axb)
    Point3D result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = -(a.x * b.z - a.z * b.x);
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

vector<double> convertStringToVector(string myInputText)
{
    stringstream ss(myInputText);
    double token;
    vector<double> tokenizeVector;
    while (ss >> token)
    {
        tokenizeVector.push_back(token);
    }
    return tokenizeVector;
}

void drawSquare(double a)
{
    glBegin(GL_QUADS);
    {
        glVertex3f(a, a, 0);
        glVertex3f(a, -a, 0);
        glVertex3f(-a, -a, 0);
        glVertex3f(-a, a, 0);
    }
    glEnd();
}

struct point
{
    double x, y, z;
};

///-----------------------------CLASS-----------------------------
class Color
{
public:
    double r, g, b;
    Color()
    {
    }
    Color(double r, double g, double b)
    {
        this->r = r;
        this->g = g;
        this->b = b;
    }

    void printColor()
    {
        cout << " red:" << r << " green:" << g << " blue:" << b << endl;
    }
};

class Ray
{
public:
    Point3D rStart;
    Point3D rDir;

    Ray() {}

    Ray(Point3D rStart, Point3D rDir)
    {
        this->rStart = rStart;
        this->rDir = rDir;
    }
};

class Object
{
public:
    Point3D refPoint;
    double height, width, length;
    Color color;
    double coEfficients[4];
    int shine;

    bool isSphere = false, isTriangle = false, isQuad = false, isFloor = false;

    Object() {}
    virtual void draw() {}

    virtual double intersect(Ray ray, Color *currentColor, int level, int objectOwnIndex)
    {
        return -1;
    }

    void setReferencePoint(double x, double y, double z)
    {
        this->refPoint.x = x;
        this->refPoint.y = y;
        this->refPoint.z = z;
    }

    Point3D getReferencePoint()
    {
        return refPoint;
    }

    void setColor(double r, double g, double b)
    {
        this->color.r = r;
        this->color.g = g;
        this->color.b = b;
    }

    Color getColor()
    {
        return color;
    }

    void setShine(int shine)
    {
        this->shine = shine;
    }

    void setCoEfficients(double a, double d, double s, double r)
    {
        this->coEfficients[AMBIENT] = a;
        this->coEfficients[DIFFUSE] = d;
        this->coEfficients[SPECULAR] = s;
        this->coEfficients[REFLECTION] = r;
    }
};

class Light
{
public:
    Point3D lightPosition;
    Color color;
    Light(){};

    Light(Point3D pos, double r, double g, double b)
    {
        this->lightPosition = pos;
        this->color.r = r;
        this->color.g = g;
        this->color.b = b;
    }
};

/// ==================variables=============================
vector<Object *> objects;
vector<Light> lights;
int lightSourceNo, levelOfRecursion, pixels, objectNo;
int imageWidth, imageHeight;

class Sphere : public Object
{
public:
    Sphere()
    {
    }

    Sphere(Point3D center, double radius, Color color)
    {
        this->refPoint = center;
        this->height = radius;
        this->width = radius;
        this->length = radius;
        this->color = color;
        this->isSphere = true;
    }

    Sphere(double centerX, double centerY, double centerZ, double radius, Color color)
    {
        this->refPoint.x = centerX;
        this->refPoint.y = centerY;
        this->refPoint.z = centerZ;
        this->height = radius;
        this->width = radius;
        this->length = radius;
        this->color = color;
        this->isSphere = true;
    }

    Point3D getNormal(Point3D intersectedPoint)
    {
        Point3D normalVector = substractTwoPointVector(intersectedPoint, this->refPoint);
        normalVector.normalizePointVector();
        return normalVector;
    }

    void draw()
    {
        //cout << "Inside Draw:Sphere" << endl;
        // struct point points[100][100];
        // int i, j;
        // double h, r;
        glPushMatrix();
        glTranslated(refPoint.x, refPoint.y, refPoint.z);
        glColor3f(color.r, color.g, color.b);
        glutSolidSphere(this->height, 100, 100);
        glPopMatrix();
        //generate points
        /*
        for (i = 0; i <= stacks; i++)
        {
            h = this->length * sin(((double)i / (double)stacks) * (pi / 2));
            r = this->length * cos(((double)i / (double)stacks) * (pi / 2));
            for (j = 0; j <= slices; j++)
            {
                points[i][j].x = r * cos(((double)j / (double)slices) * 2 * pi);
                points[i][j].y = r * sin(((double)j / (double)slices) * 2 * pi);
                points[i][j].z = h;
            }
        }
        //draw quads using generated points
        for (i = 0; i < stacks; i++)
        {
            glColor3f(this->color.r, this->color.g, this->color.b);

            for (j = 0; j < slices; j++)
            {
                glBegin(GL_QUADS);
                {
                    //upper hemisphere

                    glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
                    //lower hemisphere

                    glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
                }
                glEnd();
            }
        }

        */
    }

    double intersect(Ray ray, Color *currentColor, int level, int objectOwnIndex)
    {
        // cout << "Inside Spehere Intersect Method" << endl;
        ///logic from: https://viclw17.github.io/2018/07/16/raytracing-ray-sphere-intersection/
        //cout << "ray.rDir:" << ray.rDir.x << " " << ray.rDir.y << " " << ray.rDir.z << endl;
        //cout << "ray.rStart:" << ray.rStart.x << " " << ray.rStart.y << " " << ray.rStart.z << endl;

        double a = ray.rDir.getDotProduct(ray.rDir);
        Point3D rStart_To_reference = substractTwoPointVector(this->refPoint, ray.rStart);
        double b = 2 * ray.rDir.getDotProduct(rStart_To_reference);
        double c = rStart_To_reference.getDotProduct(rStart_To_reference) - (this->height * this->height);
        double discriminant = (b * b) - (4 * a * c);
        // cout << "discriminant: " << discriminant << endl;
        double t_root, root1, root2;

        if (discriminant >= 0)
        {
            /*         cout << "discriminant: " << discriminant << endl;
            double t1 = (-b + sqrt(discriminant)) / (2 * a);
            double t2 = (-b - sqrt(discriminant)) / (2 * a); //always minimun but sometimes negative

            double t = min(t1, t2);
            if (t > 0)
            {
                t_root = t;
            }
            else
            {
                t_root = max(t1, t2);
            }
            root is set */

            root1 = (-b - sqrt(discriminant)) / (2 * a);
            root2 = (-b + sqrt(discriminant)) / (2 * a);

            t_root == -1;
            if (root1 * root2 < 0)
                t_root = root2;
            if (root1 > 0)
                t_root = root1;

            if (root1 == 0)
                t_root = root2;

            // Point3D intersectPoint = addTwoPointVector(ray.rStart, scaleVector(ray.rDir, t_root));
            // intersectPoint.printPoint();
            if (level == 0 || t_root == -1)
                return t_root;

            ///need to complete for other levels
        }
        else
        {
            //the line of the ray does not intersect the sphere (missed);
            return -1;
        }

        Point3D intersectPoint = addTwoPointVector(ray.rStart, scaleVector(ray.rDir, t_root));
        ///--------------- Illumination Part -----------------------------
        /*Ambient Light*/
        Color outputColor = Color(currentColor->r * this->coEfficients[AMBIENT], currentColor->g * this->coEfficients[AMBIENT], currentColor->b * this->coEfficients[AMBIENT]);
        // // currentColor.r = ;
        // currentColor.g = currentColor.g * this->coEfficients[AMBIENT];
        // currentColor.b = currentColor.b * this->coEfficients[AMBIENT];
        Point3D normal;
        normal = getNormal(intersectPoint);

        Point3D cameraPosition = Point3D(pos.x, pos.y, pos.z);
        /*View Vector*/
        Point3D viewVector = substractTwoPointVector(cameraPosition, intersectPoint);
        viewVector.normalizePointVector();

        /*Diffuse & Specular Reflection */
        for (int i = 0; i < lights.size(); i++)
        {
            Light currentLight = lights[i];

            Ray lightRay = Ray(currentLight.lightPosition, substractTwoPointVector(currentLight.lightPosition, intersectPoint));
            Point3D reflectionVector = substractTwoPointVector(scaleVector(normal, normal.getDotProduct(lightRay.rDir) * 2), lightRay.rDir);
            reflectionVector.normalizePointVector();
            /// Check if the light actually reaches the object or if the light should illuminate the pixel
            bool lightIlluminated = true;
            Ray shadowRay = Ray(intersectPoint, lightRay.rDir);

            for (int j = 0; j < objects.size(); j++)
            {

                if (j != objectOwnIndex && objects[j]->intersect(shadowRay, &color, 0, -999) > 0)
                {
                    lightIlluminated = false; //ray is obstructed by anothre object
                }
                else if (j == objectOwnIndex && objects[j]->isSphere == true && lightRay.rDir.getDotProduct(normal) < 0)
                {
                    lightIlluminated = false; //light ray is opposite to intersectPoint
                }
            }
            if (lightIlluminated)
            {
                outputColor.r = outputColor.r + currentColor->r * this->coEfficients[DIFFUSE] * lightRay.rDir.getDotProduct(normal) * lights[i].color.r +
                                currentColor->r * this->coEfficients[SPECULAR] * pow(reflectionVector.getDotProduct(viewVector), this->shine) * lights[i].color.r;
                outputColor.g = outputColor.g + currentColor->g * this->coEfficients[DIFFUSE] * lightRay.rDir.getDotProduct(normal) * lights[i].color.g +
                                currentColor->g * this->coEfficients[SPECULAR] * pow(reflectionVector.getDotProduct(viewVector), this->shine) * lights[i].color.g;
                outputColor.b = outputColor.b + currentColor->b * this->coEfficients[DIFFUSE] * lightRay.rDir.getDotProduct(normal) * lights[i].color.b +
                                currentColor->b * this->coEfficients[SPECULAR] * pow(reflectionVector.getDotProduct(viewVector), this->shine) * lights[i].color.b;
            }

            currentColor->r = outputColor.r;
            currentColor->g = outputColor.g;
            currentColor->b = outputColor.b;

            if (level > levelOfRecursion)
            {
                return t_root;
            }
        }

        currentColor->r = outputColor.r;
        currentColor->g = outputColor.g;
        currentColor->b = outputColor.b;
    }
};

class Triangle : public Object
{
public:
    Point3D a, b, c;

    Triangle() {}

    Triangle(Point3D a, Point3D b, Point3D c, Color color)
    {
        this->a = a;
        this->b = b;
        this->c = c;
        this->color = color;
        this->isTriangle = true;
    }

    void draw()
    {
        glColor3f(this->color.r, this->color.g, this->color.b);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x, a.y, a.z);
            glVertex3f(b.x, b.y, b.z);
            glVertex3f(c.x, c.y, c.z);
        }
        glEnd();
    }

    double intersect(Ray ray, Color *currentColor, int level, int objectOwnIndex)
    {
        //Möller–Trumbore intersection algorithm:
        //https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm

        Point3D edge1 = substractTwoPointVector(this->b, this->a);
        Point3D edge2 = substractTwoPointVector(this->c, this->a);

        Point3D h = crossProductVector(ray.rDir, edge2);
        double A = edge1.getDotProduct(h);

        if (A > -DBL_EPSILON && A < DBL_EPSILON)
        {
            // output.push_back(-1);
            // return output;
            return -1; // This ray is parallel to this triangle.
        }

        double f = 1.0 / A;
        Point3D s = substractTwoPointVector(ray.rStart, this->a);

        double u = f * s.getDotProduct(h);

        if (u < 0.0 || u > 1.0)
        {
            // output.push_back(-1);
            // return output;
            return -1;
        }
        //    cout<<"AFter u"<<endl;

        Point3D q = crossProductVector(s, edge1);
        double v = f * ray.rDir.getDotProduct(q);

        if (v < 0.0 || u + v > 1.0)
        {
            // output.push_back(-1);
            // return output;
            return -1;
        }

        // At this stage we can compute t to find out where the intersection point is on the line.
        double root = f * edge2.getDotProduct(q);

        if (root < 0) // This means that there is a line intersection but not a ray intersection
        {
            return -1;
        }

        Point3D intersectPoint = addTwoPointVector(ray.rStart, scaleVector(ray.rDir, root));

        if (level == 0)
        {
            //cout << "root:" << root << endl;
            return root;
        }

        ///--------------- Illumination Part -----------------------------
        /*Ambient Light*/
        Color outputColor = Color(currentColor->r * this->coEfficients[AMBIENT], currentColor->g * this->coEfficients[AMBIENT], currentColor->b * this->coEfficients[AMBIENT]);
        Point3D normal;
        normal = crossProductVector(edge1, edge2);
        normal.normalizePointVector();
        Point3D cameraPosition = Point3D(pos.x, pos.y, pos.z);
        /*View Vector*/
        Point3D viewVector = substractTwoPointVector(cameraPosition, intersectPoint);
        viewVector.normalizePointVector();

        /*Diffuse & Specular Reflection*/
        // for (int i = 0; i < lights.size(); i++){

        // }
        currentColor->r = outputColor.r;
        currentColor->g = outputColor.g;
        currentColor->b = outputColor.b;
        if (level > levelOfRecursion)
        {
            return root;
        }
    }
};

class Floor : public Object
{
public:
    double floorWidth;
    double length; //tule length or size
    int tilesNo;

    Floor()
    {
    }

    Floor(double floorWidth, double length)
    {
        this->floorWidth = floorWidth;
        this->length = length;
        this->height = 0;
        this->tilesNo = floorWidth / length;
        this->refPoint = Point3D(-floorWidth / 2.0, -floorWidth / 2.0, 0);
        this->isFloor = true;
    }

    void draw()
    {

        for (int i = 0; i < tilesNo; i++)
        {
            for (int j = 0; j < tilesNo; j++)
            {
                double colorCode = (i + j) % 2;

                glColor3f(colorCode, colorCode, colorCode);

                glBegin(GL_QUADS);
                // square ABCD
                glVertex3f(refPoint.x + length * i, refPoint.y + length * j, refPoint.z);             //a
                glVertex3f(refPoint.x + length * (i + 1), refPoint.y + length * j, refPoint.z);       //b
                glVertex3f(refPoint.x + length * (i + 1), refPoint.y + length * (j + 1), refPoint.z); //c
                glVertex3f(refPoint.x + length * i, refPoint.y + length * (j + 1), refPoint.z);       //d
                glEnd();
            }
        }
    }

    Point3D getNormal(Point3D p)
    {
        p.x = 0;
        p.y = 0;
        p.z = 1;
        return p;
    }

    double intersect(Ray ray, Color *currentColor, int level, int objectOwnIndex)
    {
        Point3D normal = getNormal(refPoint);
        double up = normal.getDotProduct(ray.rStart);
        double down = normal.getDotProduct(ray.rDir);
        double t = (-1.0 * up) / down;
        if (t <= 0)
            return -1;
        if (level == 0)
            return t;
        Point3D intersectPoint = addTwoPointVector(ray.rStart, scaleVector(ray.rDir, t));

        //find out the tile : for the intersection point
        int row = abs(intersectPoint.x - refPoint.x) / this->length;
        int column = abs(intersectPoint.y - refPoint.y) / this->length;
        double colorCode = (row + column) % 2;
        Color tileColor = Color(colorCode, colorCode, colorCode);
        currentColor->r = tileColor.r * this->coEfficients[AMBIENT];
        currentColor->g = tileColor.g * this->coEfficients[AMBIENT];
        currentColor->b = tileColor.b * this->coEfficients[AMBIENT];
    }
};

class Quadric : public Object
{
public:
    vector<double> quadCoefficients;
    Quadric()
    {
    }

    Quadric(vector<double> quadCoefficients, Point3D refPoint, Color color, double l, double w, double h)
    {
        this->quadCoefficients = quadCoefficients;
        this->refPoint = refPoint;
        this->color = color;
        this->length = l;
        this->width = w;
        this->height = h;
        isQuad = true;
    }

    void draw()
    {
        double value = 20;
        glColor3f(color.r, color.g, color.b);
        glBegin(GL_QUADS);
        {
            glVertex3f(refPoint.x, refPoint.y, refPoint.z);
            glVertex3f(refPoint.x + value, refPoint.y, refPoint.z);
            glVertex3f(refPoint.x + value, refPoint.y + value, refPoint.z);
            glVertex3f(refPoint.x, refPoint.y + value, refPoint.z);
        }
        glEnd();
    }

    Point3D getIntersectPoint(Ray ray, double root)
    {
        return addTwoPointVector(ray.rStart, scaleVector(ray.rDir, root));
    }
    /*
    bool checkLengthCondition(Point3D intersectPoint)
    {
        if ((this->length > 0 && (intersectPoint.x < this->refPoint.x || intersectPoint.x > refPoint.x + this->length)) || (this->width > 0 && (intersectPoint.y < this->refPoint.y || intersectPoint.y > refPoint.y + this->width)) || (this->height > 0 && (intersectPoint.z < this->refPoint.z || intersectPoint.z > refPoint.z + this->height)))
            return true;
        else
            return false;
    } */
    bool checkLengthCondition(Point3D intersectPoint)
    {
        if ((this->length > 0 && (intersectPoint.x < this->refPoint.x || intersectPoint.x > refPoint.x + this->length)))
            return true;
        else
            return false;
    }

    bool checkWidthCondition(Point3D intersectPoint)
    {
        if ((this->width > 0 && (intersectPoint.y < this->refPoint.y || intersectPoint.y > refPoint.x + this->width)))
            return true;
        else
            return false;
    }

    bool checkHeightCondition(Point3D intersectPoint)
    {
        if ((this->height > 0 && (intersectPoint.z < this->refPoint.z || intersectPoint.z > refPoint.z + this->height)))
            return true;
        else
            return false;
    }

    Point3D getNormal(Point3D intersectPoint)
    {
        Point3D normal;
        normal.x = 2 * quadCoefficients[0] * intersectPoint.x + quadCoefficients[3] * intersectPoint.y + quadCoefficients[4] * intersectPoint.z + quadCoefficients[6];
        normal.y = 2 * quadCoefficients[1] * intersectPoint.y + quadCoefficients[3] * intersectPoint.x + quadCoefficients[5] * intersectPoint.z + quadCoefficients[7];
        normal.z = 2 * quadCoefficients[2] * intersectPoint.z + quadCoefficients[4] * intersectPoint.x + quadCoefficients[5] * intersectPoint.y + quadCoefficients[8];
        return normal;
    }

    double intersect(Ray ray, Color *currentColor, int level, int objectOwnIndex)
    {
        /* sourceLink: http://skuld.bmsc.washington.edu/people/merritt/graphics/quadrics.html */

        double Aq = quadCoefficients[0] * pow(ray.rDir.x, 2) +
                    quadCoefficients[1] * pow(ray.rDir.y, 2) +
                    quadCoefficients[2] * pow(ray.rDir.z, 2) +
                    quadCoefficients[3] * ray.rDir.x * ray.rDir.y +
                    quadCoefficients[4] * ray.rDir.x * ray.rDir.z +
                    quadCoefficients[5] * ray.rDir.y * ray.rDir.z;

        double Bq = 2 * quadCoefficients[0] * ray.rStart.x * ray.rDir.x +
                    2 * quadCoefficients[1] * ray.rStart.y * ray.rDir.y +
                    2 * quadCoefficients[2] * ray.rStart.z * ray.rDir.z +
                    quadCoefficients[3] * (ray.rStart.x * ray.rDir.y + ray.rStart.y * ray.rDir.x) +
                    quadCoefficients[4] * (ray.rStart.x * ray.rDir.z + ray.rStart.z * ray.rDir.x) +
                    quadCoefficients[5] * (ray.rStart.y * ray.rDir.z + ray.rStart.z * ray.rDir.y) +
                    quadCoefficients[6] * ray.rDir.x +
                    quadCoefficients[7] * ray.rDir.y +
                    quadCoefficients[8] * ray.rDir.z;

        double Cq = quadCoefficients[0] * ray.rStart.x * ray.rStart.x +
                    quadCoefficients[1] * ray.rStart.y * ray.rStart.y +
                    quadCoefficients[2] * ray.rStart.z * ray.rStart.z +
                    quadCoefficients[3] * ray.rStart.x * ray.rStart.y +
                    quadCoefficients[4] * ray.rStart.x * ray.rStart.z +
                    quadCoefficients[5] * ray.rStart.y * ray.rStart.z +
                    quadCoefficients[6] * ray.rStart.x +
                    quadCoefficients[7] * ray.rStart.y +
                    quadCoefficients[8] * ray.rStart.z +
                    quadCoefficients[9];

        double t1, t2, t_root;
        bool backIntersect = false;
        Point3D intersectPoint;

        if (Aq == 0)
        {
            t_root = -(Cq / Bq);
        }
        else if ((Bq * Bq - 4 * Aq * Cq) < 0)
        {
            return -1; //no intersection
        }
        else
        {
            t1 = (-Bq - (pow((Bq * Bq - 4 * Aq * Cq), 0.5))) / (2 * Aq);
            intersectPoint = getIntersectPoint(ray, t1);
            t_root = t1;

            if (t1 <= 0 || checkLengthCondition(intersectPoint) || checkWidthCondition(intersectPoint) || checkHeightCondition(intersectPoint))
            {
                t2 = (-Bq + (pow((Bq * Bq - 4 * Aq * Cq), 0.5))) / (2 * Aq);
                intersectPoint = getIntersectPoint(ray, t1);
                t_root = t2;
                if (t2 <= 0 || checkLengthCondition(intersectPoint) || checkWidthCondition(intersectPoint) || checkHeightCondition(intersectPoint))
                {
                    return -1;
                }
                else
                {
                    backIntersect = true;
                }
            }
        }

        if (level == 0)
            return t_root;

        /*Ambient Light*/
        Color outputColor = Color(currentColor->r * this->coEfficients[AMBIENT], currentColor->g * this->coEfficients[AMBIENT], currentColor->b * this->coEfficients[AMBIENT]);
        Point3D normal;
        normal = getNormal(intersectPoint);
        normal.normalizePointVector();
        Point3D cameraPosition = Point3D(pos.x, pos.y, pos.z);
        /*View Vector*/
        Point3D viewVector = substractTwoPointVector(cameraPosition, intersectPoint);
        viewVector.normalizePointVector();
    }
};

Point3D
convertStringToPoint(string myInputText)
{
    stringstream ss(myInputText);
    Point3D myPoint3D;
    int index = 0;
    double token;
    while (ss >> token)
    {
        if (index == 0)
            myPoint3D.x = token;
        else if (index == 1)
            myPoint3D.y = token;
        else if (index == 2)
            myPoint3D.z = token;
        index++;
    }
    //cout<<myPoint3D.x<<" "<<myPoint3D.y<<" "<<myPoint3D.z<<endl;
    return myPoint3D;
}
