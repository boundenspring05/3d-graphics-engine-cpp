#define WIN32_LEAN_AND_MEAN
#define UNICODE
#include <windows.h>
#include <math.h>
#include <stdint.h>
#include <conio.h>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

typedef uint32_t u32;

extern "C"
{
    #define STB_IMAGE_IMPLEMENTATION
    #include "stb_image.h"
}

using namespace std;

inline u32 RGBtoHex(short R, short G, short B)
{

    return R * 65536 + G * 256 + B;
}

bool load_image(std::vector<std::vector<unsigned int>> &image, const std::string &filename, int &x, int &y)
{
    int n;
    unsigned char *data = stbi_load(filename.c_str(), &x, &y, &n, 4);
    if (data != nullptr)
    {

        std::vector<unsigned char> image2 = std::vector<unsigned char>(data, data + x * y * 4);
        int index = 0;
        image = std::vector<std::vector<unsigned int>>(y);
        for (int i = 0; i < y; i++)
        {
            image[i] = std::vector<unsigned int>(x);
            for (int j = 0; j < x; j++)
            {
                image[i][j] = (RGBtoHex(static_cast<int>(image2[index]), static_cast<int>(image2[index + 1]), static_cast<int>(image2[index + 2])));
                index += 4;
            }
        }
        image2.clear();
    }
    stbi_image_free(data);
    return (data != nullptr);
}

bool load_image_invert(std::vector<std::vector<unsigned int>> &image, const std::string &filename, int &x, int &y)
{
    int n;
    unsigned char *data = stbi_load(filename.c_str(), &x, &y, &n, 4);
    if (data != nullptr)
    {

        std::vector<unsigned char> image2 = std::vector<unsigned char>(data, data + x * y * 4);
        int index = 0;
        image = std::vector<std::vector<unsigned int>>(y);
        for (int i = 0; i < y; i++)
        {
            image[i] = std::vector<unsigned int>(x);
            for (int j = 0; j < x; j++)
            {
                image[i][j] = (RGBtoHex(255 - static_cast<int>(image2[index]), 255 - static_cast<int>(image2[index + 1]), 255 - static_cast<int>(image2[index + 2])));
                index += 4;
            }
        }
        image2.clear();
    }
    stbi_image_free(data);
    return (data != nullptr);
}

bool load_image_lin(std::vector<unsigned char> &image, const std::string &filename, int &x, int &y)
{
    int n;
    unsigned char *data = stbi_load(filename.c_str(), &x, &y, &n, 4);
    if (data != nullptr)
    {
        /*
        std::vector<unsigned char>image2 = std::vector<unsigned char>(data, data + x * y * 4);
        int index=0;
        for(int i=0; i<x*y; i++)
        {
            image.push_back(RGBtoHex(static_cast<int>(image2[index]),static_cast<int>(image2[index+1]),static_cast<int>(image2[index+2]) ));
            index+=4;
        }
        image2.clear();
        */
        image = std::vector<unsigned char>(data, data + x * y * 4);
    }
    stbi_image_free(data);
    return (data != nullptr);
}

void *memory;
int client_width;
int client_height;
int width, height;
float uw = 1, uh = 0.00f, ud = 0.00f;
float invuw, invuh, invud;
float nearuh, nearuw;
float invnearuh, invnearuw;
float uznear = 0.01f;
float camxdist = 0, camydist = 0, camzdist = 0;
float anglex = 0, angley = 0, anglez = 0;
float lx = 0, ly = 0, lz = 1;
float pxlx = 0, pxly = 0, pxlz = 0;
float pylx = 0, pyly = 0, pylz = 0;
float sinx = 0, siny = 0, cosx = 0, cosy = 0;
bool moving = false;
bool backmoving = false;
float movespeed = 5.0f;
bool rightmov = false;
bool leftmov = false;
float fov = 60, verfov = 0;
const float PI = 3.1415926f;
std::vector<std::vector<float>> depthbuffer;
bool init = true;
DWORD currentTime, lasttime = 0;
POINT P;
int xOrigin = -1, yOrigin = -1;
float sensitivity = 20;
bool readonce = false;
struct texture
{
    std::vector<std::vector<unsigned int>> image;
    int imagewidth, imageheight;
    std::string filename;
    bool success;
};
struct texture2
{
    std::vector<unsigned char> image;
    int imagewidth, imageheight;
    std::string filename;
    bool success;
};
texture2 image1;
texture2 image2;

inline long long PerformanceCounter() noexcept
{
    LARGE_INTEGER li;
    ::QueryPerformanceCounter(&li);
    return li.QuadPart;
}

inline long long PerformanceFrequency() noexcept
{
    LARGE_INTEGER li;
    ::QueryPerformanceFrequency(&li);
    return li.QuadPart;
}

struct vertex
{
    float x, y, z;
    vertex()
    {
        x = 0, y = 0, z = 0;
    }
    vertex(float a, float b, float c)
    {
        x = a;
        y = b;
        z = c;
    }
    vertex(vertex a, vertex b, float lamda)
    {
        x = a.x + lamda * (b.x - a.x);
        y = a.y + lamda * (b.y - a.y);
        z = a.z + lamda * (b.z - a.z);
    }
    vertex(vertex a, float alpha)
    {
        x = a.x * alpha;
        y = a.y * alpha;
        z = a.z * alpha;
    }
    vertex(vertex a, vertex b)
    {
        x = a.x - b.x;
        y = a.y - b.y;
        z = a.z - b.z;
    }
};
vertex campos(0, 2, -10);
/*vertex pointdir(0.57735,-0.57735,0.57735f);*/
float invmag = 1 / sqrt(1.0f * 1.0f + 0.0f * 0.0f + 1.0f * 1.0f);
vertex pointdir(1 * invmag, 0.0 * invmag, 1 * invmag);
vertex topnormal, bottomnormal, rightnormal, leftnormal;

struct line
{
    float a, b, c;
    line()
    {
        a = 0, b = 0, c = 0;
    }
    line(float x, float y, float z)
    {
        a = x;
        b = y;
        c = z;
    }
};

vertex transformmat(vertex a)
{
    float dx = cosy * a.x - siny * a.z;
    float dy = sinx * (cosy * a.z + siny * a.x) + cosx * a.y;
    float dz = cosx * (cosy * a.z + siny * a.x) - sinx * a.y;
    return vertex(dx, dy, dz);
}

vertex crossprod(vertex a, vertex b)
{
    vertex result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

int inmod(int a)
{
    if (a < 0)
    {
        return -a;
    }
    else
    {
        return a;
    }
}
int signum(int a)
{
    if (a < 0)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}
float mod(float a)
{
    if (a < 0)
    {
        return -a;
    }
    else
    {
        return a;
    }
}
inline double dmod(double a)
{
    if (a < 0)
    {
        return -a;
    }
    else
    {
        return a;
    }
}
struct coordinate
{
    double x;
    double y;
    double z;
    coordinate()
    {
        x = 0;
        y = 0;
    }
    coordinate(double a, double b)
    {
        x = a;
        y = b;
    }
    coordinate(coordinate a, coordinate b, double lamda)
    {
        x = a.x + lamda * (b.x - a.x);
        y = a.y + lamda * (b.y - a.y);
    }
};
struct normal
{
    float x;
    float y;
    float z;
};

struct texturecoord
{
    float s;
    float t;
};

struct faces
{
    int v1 = 0, v2 = 0, v3 = 0;
    int vn1 = 0, vn2 = 0, vn3 = 0;
    int vt1 = 0, vt2 = 0, vt3 = 0;
    normal globalnormal;
    float area;
    float sineplanenormal;
    bool changetex = false;
};

int searchcharacter(string line, char a)
{
    int countl = 0;
    for (int i = 0; i < line.size(); i++)
    {
        if (line[i] == a)
        {
            countl++;
        }
    }
    return countl;
}

void replacecharacter(string &line, char a, char b)
{
    for (int i = 0; i < line.size(); i++)
    {
        if (line[i] == a)
        {
            line[i] = b;
        }
    }
}

bool searchstring(string a, string b, int pos)
{
    if (a.size() < b.size() + pos)
        return false;

    for (int i = 0; i < b.size(); i++)
    {
        if (a[i + pos] != b[i])
        {
            return false;
        }
    }
    return true;
}

float crossprodtrianglearea(vertex a, vertex b, vertex c)
{
    vertex unitvector;
    vertex vector1(a.x - b.x, a.y - b.y, a.z - b.z);

    vector1.x = vector1.x;
    vector1.y = vector1.y;
    vector1.z = vector1.z;

    vertex vector2(c.x - b.x, c.y - b.y, c.z - b.z);
    vector2.x = vector2.x;
    vector2.y = vector2.y;
    vector2.z = vector2.z;

    vertex crossprod(vector1.y * vector2.z - vector1.z * vector2.y, vector1.z * vector2.x - vector1.x * vector2.z, vector1.x * vector2.y - vector1.y * vector2.x);
    return sqrt(crossprod.x * crossprod.x + crossprod.y * crossprod.y + crossprod.z * crossprod.z) * 0.5f;
}

vertex matrirxmulti(vertex a, float matrix[3][3])
{
    float result[3];
    for (int i = 0; i < 3; i++)
    {
        result[i] = matrix[i][0] * a.x + matrix[i][1] * a.y + matrix[i][2] * a.z;
    }
    return vertex(result[0], result[1], result[2]);
}

inline float dotprod(vertex a, vertex b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vertex getBarycentricCoordinates(vertex a, vertex b, vertex c, vertex P, vertex normal, float invareaABC)
{
    vertex bary;

    // The area of a triangle is
    /*
    float areaABC = dotprod( normal, crossprod( vertex(b , a), vertex(c , a) )  ) ;
    */
    float areaPBC = dotprod(normal, crossprod(vertex(b, P), vertex(c, P)));
    float areaPCA = dotprod(normal, crossprod(vertex(c, P), vertex(a, P)));

    bary.x = dmod(areaPBC) * invareaABC; // alpha
    bary.y = dmod(areaPCA) * invareaABC; // beta
    bary.z = 1.0f - bary.x - bary.y;     // gamma

    return bary;
}

inline void draw_pixel(int x, int y, u32 color)
{
    u32 *pixel = (u32 *)memory;
    pixel += y * client_width + x;
    *pixel = color;
}

void drawline2(coordinate a, coordinate b)
{
    float dirx = b.x - a.x;
    float diry = b.y - a.y;
    float distance = sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
    for (int i = 0; i <= distance; i++)
    {
        draw_pixel(a.x + i * dirx / distance, a.y + i * diry / distance, 0xffffff);
    }
}
void drawline(coordinate a, coordinate b, u32 color)
{
    float dirx = b.x - a.x;
    float diry = b.y - a.y;
    if (mod(diry) > mod(dirx))
    {
        if (a.y < b.y)
        {
            for (int i = a.y; i <= b.y; i++)
            {
                draw_pixel(a.x + ((i - a.y) / (diry)) * dirx, i, color);
            }
        }
        else
        {
            for (int i = a.y; i >= b.y; i--)
            {
                draw_pixel(a.x + ((i - a.y) / (diry)) * dirx, i, color);
            }
        }
    }
    else
    {
        if (a.x < b.x)
        {
            for (int i = a.x; i <= b.x; i++)
            {
                draw_pixel(i, a.y + ((i - a.x) / (dirx)) * diry, color);
            }
        }
        else
        {
            for (int i = a.x; i >= b.x; i--)
            {
                draw_pixel(i, a.y + ((i - a.x) / (dirx)) * diry, color);
            }
        }
    }
}

void BresenhamInt(int x1, int y1, int x2, int y2)
{
    bool changed = false;
    int x = x1;
    int y = y1;
    int dx = inmod(x2 - x1);
    int dy = inmod(y2 - y1);
    int signx = signum(x2 - x1);
    int signy = signum(y2 - y1);
    if (dy > dx)
    {
        std::swap(dx, dy);
        changed = true;
    }
    float e = (dy + dy) - dx;
    for (int i = 1; i <= dx; i++)
    {
        draw_pixel(x, y, 0xffffff);
        while (e >= 0)
        {
            if (changed)
            {
                x = x + signx;
            }
            else
            {
                y = y + signy;
            }
            e = e - (dx + dx);
        }
        if (changed)
        {
            y += signy;
        }
        else
        {
            x += signx;
        }
        e = e + (dy + dy);
    }
}

void drawsquare(int x, int y, int width)
{
    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= width; j++)
        {
            draw_pixel(x + j, y + i, 0xff0000);
        }
    }
}

inline bool check2linelineintersection(vertex line1, vertex line2, vertex point, vertex point2)
{
    vertex dir(line2.x - line1.x, line2.y - line1.y, line2.z - line1.z);
    vertex line(point2.x - point.x, point2.y - point.y, point2.z - point.z);
    vertex intdir(line1.x - point.x, line1.y - point.y, line1.z - point.z);
    float d1321 = intdir.x * dir.x + intdir.y * dir.y + intdir.z * dir.z;
    float d2121 = dir.x * dir.x + dir.y * dir.y + dir.z * dir.z;
    float d4321 = line.x * dir.x + line.y * dir.y + line.z * dir.z;
    float d1343 = intdir.x * line.x + intdir.y * line.y + intdir.z * line.z;
    float d4343 = line.x * line.x + line.y * line.y + line.z * line.z;
    if ((d2121 * d4343 - d4321 * d4321) == 0)
    {
        return false;
    }
    float k = (d1343 * d4321 - d1321 * d4343) / (d2121 * d4343 - d4321 * d4321);
    if (k < 1 && k > 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline bool pointtriangle(vertex p1, vertex p2, vertex p3, vertex cp, float area)
{
    if (area < 0.0001)
    {
        return false;
    }
    vertex centroid((p1.x + p2.x + p3.x) / 3, (p1.y + p2.y + p3.y) / 3, (p1.z + p2.z + p3.z) / 3);
    if (check2linelineintersection(centroid, cp, p1, p2) == true || check2linelineintersection(centroid, cp, p2, p3) == true || check2linelineintersection(centroid, cp, p3, p1) == true)
    {
        return false;
    }
    else
    {
        return true;
    }
}

double Cross(coordinate a, coordinate b)
{
    return a.x * b.y - a.y * b.x;
}

inline bool GetWeights(coordinate A, coordinate B, coordinate C, coordinate point)
{
    double w_A = 0, w_B = 0, w_C = 0;
    double crossab = Cross(A, B);
    double crossbc = Cross(B, C);
    double crossca = Cross(C, A);
    double xd = crossab + crossbc + crossca;
    double invxd = 1 / xd;
    if (dmod(xd) > 1e-13)
    {
        double xa = crossbc + Cross(point, coordinate(B.x - C.x, B.y - C.y));
        double xb = crossca + Cross(point, coordinate(C.x - A.x, C.y - A.y));

        w_A = xa * invxd;
        w_B = xb * invxd;
        w_C = 1 - (w_A + w_B);
        if (w_A >= 0 && w_A <= 1 && w_B >= 0 && w_B <= 1 && w_C >= 0 && w_C <= 1)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    w_A = 0;
    w_B = 0;
    w_C = 0;
    return false;
}

vertex rayplane(float nx, float ny, float nz, vertex point, vertex v1, vertex v2)
{
    float dx = v2.x - v1.x;
    float dy = v2.y - v1.y;
    float dz = v2.z - v1.z;
    float lamda = dotprod(vertex(point.x - v1.x, point.y - v1.y, point.z - v1.z), vertex(nx, ny, nz)) / dotprod(vertex(dx, dy, dz), vertex(nx, ny, nz));
    return vertex(v1.x + lamda * dx, v1.y + lamda * dy, v1.z + lamda * dz);
}
float rayplanelamda(float nx, float ny, float nz, vertex point, vertex v1, vertex v2)
{
    float dx = v2.x - v1.x;
    float dy = v2.y - v1.y;
    float dz = v2.z - v1.z;
    float lamda = dotprod(vertex(point.x - v1.x, point.y - v1.y, point.z - v1.z), vertex(nx, ny, nz)) / dotprod(vertex(dx, dy, dz), vertex(nx, ny, nz));
    return lamda;
}

vertex rayplaneIntersect(vertex planeN, vertex planeP, vertex rayP, vertex rayD)
{
    float d = dotprod(planeP, vertex(-planeN.x, -planeN.y, -planeN.z));
    float t = -(d + rayP.z * planeN.z + rayP.y * planeN.y + rayP.x * planeN.x) / (rayD.z * planeN.z + rayD.y * planeN.y + rayD.x * planeN.x);
    return vertex(rayP.x + t * rayD.x, rayP.z + t * rayD.z, rayP.z + t * rayD.z);
}

void setfrustrum(vertex camerapos, float camanglex2, float camangley2, float camanglez2, float fov, float znear, float aspect)
{
    cosx = cos(camanglex2 * PI / 180);
    sinx = sin(camanglex2 * PI / 180);
    cosy = cos(camangley2 * PI / 180);
    siny = sin(camangley2 * PI / 180);

    ud = uw / tan((fov / 2) * PI / 180);
    uh = uw / aspect;
    invuh = 1 / uh;
    invuw = 1 / uw;
    invud = 1 / ud;
    lx = 0;
    ly = 0;
    lz = 1;
    nearuh = uh * znear * invud;
    nearuw = uw * znear * invud;
    invnearuh = 1 / nearuh;
    invnearuw = 1 / nearuw;
    uznear = znear;
}

void draw2dtriangle2(coordinate a, coordinate b, coordinate c)
{
    coordinate tri[3] = {a, b, c};
    int maxy = tri[0].y;
    int miny = tri[0].y;
    int maxyi = 0;
    int minyi = 0;
    int maxx = tri[0].x;
    int minx = tri[0].x;
    int maxxi = 0;
    int minxi = 0;
    for (int i = 0; i < 3; i++)
    {
        if (tri[i].y > maxy)
        {
            maxyi = i;
            maxy = tri[i].y;
        }
        if (tri[i].y < miny)
        {
            minyi = i;
            miny = tri[i].y;
        }
        if (tri[i].x > maxx)
        {
            maxxi = i;
            maxx = tri[i].x;
        }
        if (tri[i].x < minx)
        {
            minxi = i;
            minx = tri[i].x;
        }
    }

    for (int i = miny; i <= maxy; i++)
    {
        for (int j = minx; j <= maxx; j++)
        {
            if (pointtriangle(vertex(a.x, a.y, 0), vertex(b.x, b.y, 0), vertex(c.x, c.y, 0), vertex(j, i, 0), 1) == true)
            {
                draw_pixel(j, i, 0xffffff);
            }
        }
    }
}

void draw2dtriangle3(coordinate a, coordinate b, coordinate c, int winwidth, int winheight)
{
    coordinate tri[3] = {a, b, c};
    int maxy = tri[0].y;
    int miny = tri[0].y;
    int maxyi = 0;
    int minyi = 0;
    int maxx = tri[0].x;
    int minx = tri[0].x;
    int maxxi = 0;
    int minxi = 0;
    for (int i = 1; i < 3; i++)
    {
        if (tri[i].y > maxy)
        {
            maxyi = i;
            maxy = tri[i].y;
        }
        if (tri[i].y < miny)
        {
            minyi = i;
            miny = tri[i].y;
        }
        if (tri[i].x > maxx)
        {
            maxxi = i;
            maxx = tri[i].x;
        }
        if (tri[i].x < minx)
        {
            minxi = i;
            minx = tri[i].x;
        }
    }
    if (minx < 0)
    {
        minx = 0;
    }
    if (maxx > width - 1)
    {
        maxx = width - 1;
    }
    if (miny < 0)
    {
        miny = 0;
    }
    if (maxy > height - 1)
    {
        maxy = height - 1;
    }
    double w_A = 0, w_B = 0, w_C = 0;
    double crossab = Cross(a, b);
    double crossbc = Cross(b, c);
    double crossca = Cross(c, a);
    double xd = crossab + crossbc + crossca;
    double invxd = 1 / xd;
    for (int i = miny; i <= maxy; i++)
    {
        for (int j = minx; j <= maxx; j++)
        {
            bool pass = false;
            coordinate point = coordinate(j, i);
            if (dmod(xd) > 1e-13)
            {
                double xa = crossbc + Cross(point, coordinate(b.x - c.x, b.y - c.y));
                double xb = crossca + Cross(point, coordinate(c.x - a.x, c.y - a.y));

                w_A = xa * invxd;
                w_B = xb * invxd;
                w_C = 1 - (w_A + w_B);
                if (w_A >= 0 && w_A <= 1 && w_B >= 0 && w_B <= 1 && w_C >= 0 && w_C <= 1)
                {
                    pass = true;
                }
                else
                {
                    pass = false;
                }
            }
            else
            {
                w_A = 0;
                w_B = 0;
                w_C = 0;
            }
            if (pass == true)
            {

                draw_pixel(j, i, 0xffffff);
            }
        }
    }
}

void draw2dtriangleup(coordinate a, coordinate b, coordinate c, int winwidth, int winheight, u32 color)
{

    coordinate tri[3] = {a, b, c};
    int maxy = tri[0].y;
    int miny = tri[0].y;
    int maxyi = 0;
    int minyi = 0;
    int maxx = tri[0].x;
    int minx = tri[0].x;
    int maxxi = 0;
    int minxi = 0;
    for (int i = 1; i < 3; i++)
    {
        if (tri[i].y > maxy)
        {
            maxyi = i;
            maxy = tri[i].y;
        }
        if (tri[i].y < miny)
        {
            minyi = i;
            miny = tri[i].y;
        }
        if (tri[i].x > maxx)
        {
            maxxi = i;
            maxx = tri[i].x;
        }
        if (tri[i].x < minx)
        {
            minxi = i;
            minx = tri[i].x;
        }
    }

    if (tri[maxyi].y - tri[minyi].y < 2 || tri[maxxi].x - tri[minxi].x < 2)
    {
    }
    else
    {

        int mediumyi = 3 - (minyi + maxyi);
        int mediumy = tri[mediumyi].y;

        bool achanged = false;
        int ax1 = 0;
        int ay1 = 0;
        float az1 = 0;
        if (tri[minyi].x < tri[mediumyi].x)
        {
            ax1 = tri[minyi].x;
            ay1 = tri[minyi].y;
            az1 = tri[minyi].z;
        }
        else
        {
            ax1 = tri[mediumyi].x;
            ay1 = tri[mediumyi].y;
            az1 = tri[mediumyi].z;
        }
        int ax2 = tri[maxyi].x;
        int ay2 = tri[maxyi].y;
        float az2 = tri[maxyi].z;
        int ax = ax1;
        int ay = ay1;
        float az = az1;
        float invaz = 1 / az;
        int recax = ax;
        int recay = ay;
        float recaz = az;
        float invaz1 = 1 / az1;
        float invaz2 = 1 / az2;

        int adx = inmod(ax2 - ax1);
        int ady = inmod(ay2 - ay1);
        int asignx = signum(ax2 - ax1);
        int asigny = signum(ay2 - ay1);
        int recai = 1;
        if (ady > adx)
        {
            std::swap(adx, ady);
            achanged = true;
        }
        float invadx;
        if (achanged == false)
        {
            invadx = 1 / float(adx);
        }
        else
        {
            invadx = 1 / float(adx);
        }
        float aconst = invadx * (invaz2 - invaz1);
        float ae = (ady + ady) - adx;

        bool bchanged = false;
        int bx1 = 0;
        int by1 = 0;
        float bz1 = 0;
        if (tri[minyi].x < tri[mediumyi].x)
        {
            bx1 = tri[mediumyi].x;
            by1 = tri[mediumyi].y;
            bz1 = tri[mediumyi].z;
        }
        else
        {
            bx1 = tri[minyi].x;
            by1 = tri[minyi].y;
            bz1 = tri[minyi].z;
        }
        int bx2 = tri[maxyi].x;
        int by2 = tri[maxyi].y;
        float bz2 = tri[maxyi].z;
        int bx = bx1;
        int by = by1;
        float bz = bz1;
        float invbz1 = 1 / bz1;
        float invbz2 = 1 / bz2;
        float invbz = 1 / bz;
        int bdx = inmod(bx2 - bx1);
        int bdy = inmod(by2 - by1);
        int bsignx = signum(bx2 - bx1);
        int bsigny = signum(by2 - by1);
        int recbi = 1;

        if (bdy > bdx)
        {
            std::swap(bdx, bdy);
            bchanged = true;
        }
        float invbdx;
        if (bchanged == false)
        {
            invbdx = 1 / float(bdx);
        }
        else
        {
            invbdx = 1 / float(bdx);
        }
        float bconst = invbdx * (invbz2 - invbz1);
        float be = (bdy + bdy) - bdx;
        if (miny < 0)
        {
            miny = 0;
        }
        if (miny > height - 1)
        {
            miny = height - 1;
        }
        if (maxy < 0)
        {
            maxy = 0;
        }
        if (maxy > height - 1)
        {
            maxy = height - 1;
        }
        for (int i = miny; i <= maxy; i++)
        {

            for (int ai = recai; ai <= adx; ai++)
            {
                if (ay == i)
                {
                    recai = ai;
                    break;
                }
                invaz = invaz + aconst;
                while (ae >= 0)
                {
                    if (achanged)
                    {
                        ax = ax + asignx;
                    }
                    else
                    {
                        ay = ay + asigny;
                    }
                    ae = ae - (adx + adx);
                }
                if (achanged)
                {
                    ay += asigny;
                }
                else
                {
                    ax += asignx;
                }
                ae = ae + (ady + ady);
            }

            for (int bi = recbi; bi <= bdx; bi++)
            {
                if (by == i)
                {
                    recbi = bi;
                    break;
                }
                invbz = invbz + bconst;
                while (be >= 0)
                {
                    if (bchanged)
                    {
                        bx = bx + bsignx;
                    }
                    else
                    {
                        by = by + bsigny;
                    }
                    be = be - (bdx + bdx);
                }
                if (bchanged)
                {
                    by += bsigny;
                }
                else
                {
                    bx += bsignx;
                }
                be = be + (bdy + bdy);
            }
            int tempax = ax;
            int tempbx = bx;
            if (tempax > winwidth - 1)
            {
                tempax = winwidth - 1;
            }
            if (tempax < 0)
            {
                tempax = 0;
            }
            if (tempbx > winwidth - 1)
            {
                tempbx = winwidth - 1;
            }
            if (tempbx < 0)
            {
                tempbx = 0;
            }

            float tempinvdx = 1;
            if (tempax != tempbx)
            {
                tempinvdx = 1 / float(tempbx - tempax);
            }
            float tempinvz = invaz;
            float tempconst = tempinvdx * (invbz - invaz);
            float tempz;
            for (int j = tempax; j < tempbx; j++)
            {
                tempz = 1 / tempinvz;

                if (tempinvz > depthbuffer[i][j])
                {
                    depthbuffer[i][j] = tempinvz;

                    draw_pixel(j, i, color);
                }
                tempinvz = tempinvz + tempconst;
            }
            recax = ax;
            recay = ay;
            recaz = az;
        }
    }
}

void draw2dtriangledown(coordinate a, coordinate b, coordinate c, int winwidth, int winheight, u32 color)
{

    coordinate tri[3] = {a, b, c};
    int maxy = tri[0].y;
    int miny = tri[0].y;
    int maxyi = 0;
    int minyi = 0;
    int maxx = tri[0].x;
    int minx = tri[0].x;
    int maxxi = 0;
    int minxi = 0;
    for (int i = 1; i < 3; i++)
    {
        if (tri[i].y > maxy)
        {
            maxyi = i;
            maxy = tri[i].y;
        }
        if (tri[i].y < miny)
        {
            minyi = i;
            miny = tri[i].y;
        }
        if (tri[i].x > maxx)
        {
            maxxi = i;
            maxx = tri[i].x;
        }
        if (tri[i].x < minx)
        {
            minxi = i;
            minx = tri[i].x;
        }
    }
    if (tri[maxyi].y - tri[minyi].y >= 2 && tri[maxxi].x - tri[minxi].x >= 2)
    {
        int mediumyi = 3 - (minyi + maxyi);
        int mediumy = tri[mediumyi].y;
        bool achanged = false;
        int ax1 = 0;
        int ay1 = 0;
        float az1 = 0;
        if (tri[maxyi].x < tri[mediumyi].x)
        {
            ax1 = tri[maxyi].x;
            ay1 = tri[maxyi].y;
            az1 = tri[maxyi].z;
        }
        else
        {
            ax1 = tri[mediumyi].x;
            ay1 = tri[mediumyi].y;
            az1 = tri[mediumyi].z;
        }

        int ax2 = tri[minyi].x;
        int ay2 = tri[minyi].y;
        float az2 = tri[minyi].z;

        int ax = ax1;
        int ay = ay1;
        float az = az1;
        float invaz = 1 / az;
        float invaz1 = 1 / az1;
        float invaz2 = 1 / az2;
        int adx = inmod(ax2 - ax1);
        int ady = inmod(ay2 - ay1);
        int asignx = signum(ax2 - ax1);
        int asigny = signum(ay2 - ay1);
        int recai = 1;
        if (ady > adx)
        {
            std::swap(adx, ady);
            achanged = true;
        }

        float invadx;
        if (achanged == false)
        {
            invadx = 1 / float(adx);
        }
        else
        {
            invadx = 1 / float(adx);
        }
        float aconst = invadx * (invaz2 - invaz1);

        float ae = (ady + ady) - adx;

        bool bchanged = false;
        int bx1 = 0;
        int by1 = 0;
        float bz1 = 0;
        if (tri[maxyi].x < tri[mediumyi].x)
        {
            bx1 = tri[mediumyi].x;
            by1 = tri[mediumyi].y;
            bz1 = tri[mediumyi].z;
        }
        else
        {
            bx1 = tri[maxyi].x;
            by1 = tri[maxyi].y;
            bz1 = tri[maxyi].z;
        }
        int bx2 = tri[minyi].x;
        int by2 = tri[minyi].y;
        float bz2 = tri[minyi].z;
        int bx = bx1;
        int by = by1;
        float bz = bz1;
        float invbz = 1 / bz;
        float invbz1 = 1 / bz1;
        float invbz2 = 1 / bz2;
        int bdx = inmod(bx2 - bx1);
        int bdy = inmod(by2 - by1);
        int bsignx = signum(bx2 - bx1);
        int bsigny = signum(by2 - by1);
        int recbi = 1;

        if (bdy > bdx)
        {
            std::swap(bdx, bdy);
            bchanged = true;
        }
        float invbdx = 0;
        if (bchanged == false)
        {
            invbdx = 1 / float(bdx);
        }
        else
        {
            invbdx = 1 / float(bdx);
        }
        float bconst = invbdx * (invbz2 - invbz1);

        float be = (bdy + bdy) - bdx;
        if (miny < 0)
        {
            miny = 0;
        }
        if (miny > height - 1)
        {
            miny = height - 1;
        }
        if (maxy < 0)
        {
            maxy = 0;
        }
        if (maxy > height - 1)
        {
            maxy = height - 1;
        }
        for (int i = maxy; i >= miny; i--)
        {

            for (int ai = recai; ai <= adx; ai++)
            {
                if (ay == i)
                {
                    recai = ai;
                    break;
                }
                invaz = invaz + aconst;
                while (ae >= 0)
                {
                    if (achanged)
                    {
                        ax = ax + asignx;
                    }
                    else
                    {
                        ay = ay + asigny;
                    }
                    ae = ae - (adx + adx);
                }
                if (achanged)
                {
                    ay += asigny;
                }
                else
                {
                    ax += asignx;
                }
                ae = ae + (ady + ady);
            }

            for (int bi = recbi; bi <= bdx; bi++)
            {
                if (by == i)
                {
                    recbi = bi;
                    break;
                }
                invbz = invbz + bconst;
                while (be >= 0)
                {
                    if (bchanged)
                    {
                        bx = bx + bsignx;
                    }
                    else
                    {
                        by = by + bsigny;
                    }
                    be = be - (bdx + bdx);
                }
                if (bchanged)
                {
                    by += bsigny;
                }
                else
                {
                    bx += bsignx;
                }
                be = be + (bdy + bdy);
            }
            int tempax = ax;
            int tempbx = bx;
            if (tempax > winwidth - 1)
            {
                tempax = winwidth - 1;
            }
            if (tempax < 0)
            {
                tempax = 0;
            }
            if (tempbx > winwidth - 1)
            {
                tempbx = winwidth - 1;
            }
            if (tempbx < 0)
            {
                tempbx = 0;
            }
            float tempinvdx = 1;
            if (tempax != tempbx)
            {
                tempinvdx = 1 / float(tempbx - tempax);
            }
            float tempinvz = invaz;
            float tempconst = tempinvdx * (invbz - invaz);
            float tempz;
            for (int j = tempax; j < tempbx; j++)
            {
                tempz = 1 / tempinvz;

                if (tempz < depthbuffer[i][j])
                {

                    depthbuffer[i][j] = tempz;
                    draw_pixel(j, i, color);
                }
                tempinvz = tempinvz + tempconst;
            }
        }
    }
}

void draw2dtriangle(coordinate a, coordinate b, coordinate c, int winwidth, int winheight, u32 color)
{
    if (a.x < 0 && b.x < 0 && c.x < 0 || a.x > winwidth && b.x > winwidth && c.x > winwidth || a.y < 0 && b.y < 0 && c.y < 0 || a.y > winheight && b.y > winheight && c.y > winheight)
    {
    }
    else
    {
        coordinate tri[3] = {a, b, c};
        int maxy = tri[0].y;
        int miny = tri[0].y;
        int maxyi = 0;
        int minyi = 0;
        int maxx = tri[0].x;
        int minx = tri[0].x;
        int maxxi = 0;
        int minxi = 0;
        for (int i = 1; i < 3; i++)
        {
            if (tri[i].y > maxy)
            {
                maxyi = i;
                maxy = tri[i].y;
            }
            if (tri[i].y < miny)
            {
                minyi = i;
                miny = tri[i].y;
            }
            if (tri[i].x > maxx)
            {
                maxxi = i;
                maxx = tri[i].x;
            }
            if (tri[i].x < minx)
            {
                minxi = i;
                minx = tri[i].x;
            }
        }

        int mediumyi = 3 - (minyi + maxyi);
        int mediumy = tri[mediumyi].y;
        coordinate d;
        d.y = tri[mediumyi].y;
        d.x = tri[maxyi].x + (tri[minyi].x - tri[maxyi].x) * ((tri[mediumyi].y - tri[maxyi].y) / (tri[minyi].y - tri[maxyi].y));
        d.x = int(d.x);
        if (inmod(tri[maxyi].x - tri[minyi].x) > inmod(tri[maxyi].y - tri[minyi].y))
        {
            float alpha = (d.x - tri[minyi].x) / (tri[maxyi].x - tri[minyi].x);
            d.z = alpha * (1 / tri[maxyi].z) + (1 - alpha) * (1 / tri[minyi].z);
            d.z = 1 / d.z;
        }
        else
        {
            float alpha = (d.y - tri[minyi].y) / (tri[maxyi].y - tri[minyi].y);
            d.z = alpha * (1 / tri[maxyi].z) + (1 - alpha) * (1 / tri[minyi].z);
            d.z = 1 / d.z;
        }
        draw2dtriangleup(tri[maxyi], tri[mediumyi], d, winwidth, winheight, color);

        draw2dtriangledown(tri[minyi], tri[mediumyi], d, winwidth, winheight, color);
    }
}

void draw3dtriangle(vertex a2, vertex b2, vertex c2, int winwidth, int winheight, u32 color)
{

    vertex tempa = transformmat(vertex(a2.x - campos.x, a2.y - campos.y, a2.z - campos.z));
    vertex tempb = transformmat(vertex(b2.x - campos.x, b2.y - campos.y, b2.z - campos.z));
    vertex tempc = transformmat(vertex(c2.x - campos.x, c2.y - campos.y, c2.z - campos.z));
    /*
    vertex tempa=matrirxmulti(vertex(a2.x-campos.x,a2.y-campos.y,a2.z-campos.z),matrix);
    vertex tempb=matrirxmulti(vertex(b2.x-campos.x,b2.y-campos.y,b2.z-campos.z),matrix);
    vertex tempc=matrirxmulti(vertex(c2.x-campos.x,c2.y-campos.y,c2.z-campos.z),matrix);
    */

    coordinate a2d, b2d, c2d;
    vertex posa, posb, posc;
    float adistx, adisty;
    float bdistx, bdisty;
    float cdistx, cdisty;
    float adist = tempa.z;
    float bdist = tempb.z;
    float cdist = tempc.z;

    if (adist <= uznear || bdist <= uznear || cdist <= uznear)
    {

        if (adist <= uznear && bdist <= uznear && cdist <= uznear || adist <= 0 && bdist <= 0 && cdist <= 0)
        {
        }
        else
        {
            vertex tri[3] = {tempa, tempb, tempc};
            int minzi = 0;
            float minz = tempa.z;
            int maxzi = 0;
            float maxz = tempa.z;
            for (int i = 0; i < 3; i++)
            {
                if (tri[i].z < minz)
                {
                    minz = tri[i].z;
                    minzi = i;
                }
                if (tri[i].z > maxz)
                {
                    maxz = tri[i].z;
                    maxzi = i;
                }
            }
            int mediumzi = 3 - minzi - maxzi;
            float mediumz = tri[mediumzi].z;
            if (maxz > uznear)
            {

                if (mediumz > uznear)
                {
                    coordinate min2d, medium2d, max2d;
                    float invmediumz = 1 / mediumz;
                    float invmaxz = 1 / maxz;
                    vertex tempposmed = vertex(tri[mediumzi].x * ud * invmediumz, tri[mediumzi].y * ud * invmediumz, ud);
                    vertex tempposmax = vertex(tri[maxzi].x * ud * invmaxz, tri[maxzi].y * ud * invmaxz, ud);
                    float meddistx = tempposmed.x * invuw;
                    float meddisty = tempposmed.y * invuh;
                    float maxdistx = tempposmax.x * invuw;
                    float maxdisty = tempposmax.y * invuh;
                    medium2d.x = winwidth * (meddistx + 1) * 0.5;
                    medium2d.y = winheight * (meddisty + 1) * 0.5;
                    max2d.x = winwidth * (maxdistx + 1) * 0.5;
                    max2d.y = winheight * (maxdisty + 1) * 0.5;
                    medium2d.z = mediumz;
                    max2d.z = maxz;

                    float lamdamax = rayplanelamda(0, 0, 1, vertex(0, 0, uznear), tri[maxzi], tri[minzi]);
                    float lamdamed = rayplanelamda(0, 0, 1, vertex(0, 0, uznear), tri[mediumzi], tri[minzi]);

                    vertex maxnear = vertex(tri[maxzi], tri[minzi], lamdamax);
                    coordinate nearmax2d;
                    float nearmax2ddistx = maxnear.x * invnearuw;
                    float nearmax2ddisty = maxnear.y * invnearuh;
                    nearmax2d.x = winwidth * (nearmax2ddistx + 1) * 0.5;
                    nearmax2d.y = winheight * (nearmax2ddisty + 1) * 0.5;
                    nearmax2d.z = maxnear.z;

                    vertex mednear = vertex(tri[mediumzi], tri[minzi], lamdamed);
                    coordinate nearmed2d;
                    float nearmed2ddistx = mednear.x * invnearuw;
                    float nearmed2ddisty = mednear.y * invnearuh;
                    nearmed2d.x = winwidth * (nearmed2ddistx + 1) * 0.5;
                    nearmed2d.y = winheight * (nearmed2ddisty + 1) * 0.5;
                    nearmed2d.z = mednear.z;
                    draw2dtriangle(nearmed2d, medium2d, max2d, winwidth, winheight, color);
                    draw2dtriangle(nearmed2d, nearmax2d, max2d, winwidth, winheight, color);
                }
                else
                {
                    coordinate max2d;
                    float invmaxz = 1 / maxz;
                    vertex tempposmax = vertex(tri[maxzi].x * ud * invmaxz, tri[maxzi].y * ud * invmaxz, ud);
                    float maxdistx = tempposmax.x * invuw;
                    float maxdisty = tempposmax.y * invuh;
                    max2d.x = winwidth * (maxdistx + 1) * 0.5;
                    max2d.y = winheight * (maxdisty + 1) * 0.5;
                    max2d.z = maxz;

                    float lamdamed = rayplanelamda(0, 0, 1, vertex(0, 0, uznear), tri[maxzi], tri[mediumzi]);
                    float lamdamin = rayplanelamda(0, 0, 1, vertex(0, 0, uznear), tri[maxzi], tri[minzi]);

                    vertex mednear = vertex(tri[maxzi], tri[mediumzi], lamdamed);
                    coordinate nearmed2d;
                    float nearmed2ddistx = mednear.x * invnearuw;
                    float nearmed2ddisty = mednear.y * invnearuh;
                    nearmed2d.x = winwidth * (nearmed2ddistx + 1) * 0.5;
                    nearmed2d.y = winheight * (nearmed2ddisty + 1) * 0.5;
                    nearmed2d.z = mednear.z;

                    vertex minnear = vertex(tri[maxzi], tri[minzi], lamdamin);
                    coordinate nearmin2d;
                    float nearmin2ddistx = minnear.x * invnearuw;
                    float nearmin2ddisty = minnear.y * invnearuh;
                    nearmin2d.x = winwidth * (nearmin2ddistx + 1) * 0.5;
                    nearmin2d.y = winheight * (nearmin2ddisty + 1) * 0.5;
                    nearmin2d.z = minnear.z;

                    draw2dtriangle(max2d, nearmin2d, nearmed2d, winwidth, winheight, color);
                }
            }
            else
            {
            }
        }
    }
    else
    {

        /*
        posa=rayplane(lx,ly, lz, camzpoint,campos, vertex(a.x,a.y,a.z));
        posb=rayplane(lx,ly, lz, camzpoint,campos, vertex(b.x,b.y,b.z));
        posc=rayplane(lx,ly, lz, camzpoint,campos, vertex(c.x,c.y,c.z));
        */
        float invadist = 1 / adist;
        float invbdist = 1 / bdist;
        float invcdist = 1 / cdist;
        vertex tempposa = vertex(tempa.x * ud * invadist, tempa.y * ud * invadist, ud);
        vertex tempposb = vertex(tempb.x * ud * invbdist, tempb.y * ud * invbdist, ud);
        vertex tempposc = vertex(tempc.x * ud * invcdist, tempc.y * ud * invcdist, ud);

        adistx = tempposa.x * invuw;
        adisty = tempposa.y * invuh;
        bdistx = tempposb.x * invuw;
        bdisty = tempposb.y * invuh;
        cdistx = tempposc.x * invuw;
        cdisty = tempposc.y * invuh;
        a2d.x = winwidth * (adistx + 1) * 0.5;
        a2d.y = winheight * (adisty + 1) * 0.5;
        b2d.x = winwidth * (bdistx + 1) * 0.5;
        b2d.y = winheight * (bdisty + 1) * 0.5;
        c2d.x = winwidth * (cdistx + 1) * 0.5;
        c2d.y = winheight * (cdisty + 1) * 0.5;
        a2d.z = tempa.z;
        b2d.z = tempb.z;
        c2d.z = tempc.z;
        draw2dtriangle(a2d, b2d, c2d, winwidth, winheight, color);
    }
}

void draw2dtriangleuptex(coordinate a, coordinate b, coordinate c, int winwidth, int winheight, coordinate uvtex1, coordinate uvtex2, coordinate uvtex3, vertex normal, vertex lightdir, float diffuse, float ambient, std::vector<unsigned char> &image, float imgwidth, float imgheight)
{
    float cosine = -dotprod(normal, lightdir);
    vertex reflection;
    reflection.x = cosine * 2.0f * normal.x + lightdir.x;
    reflection.y = cosine * 2.0f * normal.y + lightdir.y;
    reflection.z = cosine * 2.0f * normal.z + lightdir.z;
    float invreflect_mag = reflection.x * reflection.x + reflection.y * reflection.y + reflection.z * reflection.z;
    invreflect_mag = 1 / invreflect_mag;
    float invheight = 1 / float(height);
    float invwidth = 1 / float(width);
    if (cosine < 0)
    {
        cosine = 0;
    }

    /*
    float coeff=ambient+cosine*(diffuse-ambient);
    */
    float coeff = ambient + cosine * (diffuse - ambient);
    if (coeff > 1)
    {
        coeff = 1;
    }
    coordinate tex1(uvtex1.x * imgwidth, uvtex1.y * imgheight);
    coordinate tex2(uvtex2.x * imgwidth, uvtex2.y * imgheight);
    coordinate tex3(uvtex3.x * imgwidth, uvtex3.y * imgheight);
    coordinate tri[3] = {a, b, c};
    coordinate tritex[3] = {tex1, tex2, tex3};
    int maxy = tri[0].y;
    int miny = tri[0].y;
    BYTE maxyi = 0;
    BYTE minyi = 0;
    int maxx = tri[0].x;
    int minx = tri[0].x;
    BYTE maxxi = 0;
    BYTE minxi = 0;
    for (BYTE i = 1; i < 3; i++)
    {
        if (tri[i].y > maxy)
        {
            maxyi = i;
            maxy = tri[i].y;
        }
        if (tri[i].y < miny)
        {
            minyi = i;
            miny = tri[i].y;
        }
        if (tri[i].x > maxx)
        {
            maxxi = i;
            maxx = tri[i].x;
        }
        if (tri[i].x < minx)
        {
            minxi = i;
            minx = tri[i].x;
        }
    }

    if (a.x < 0 && b.x < 0 && c.x < 0 || a.x > winwidth && b.x > winwidth && c.x > winwidth || a.y < 0 && b.y < 0 && c.y < 0 || a.y > winheight && b.y > winheight && c.y > winheight || tri[maxyi].y - tri[minyi].y < 1 || tri[maxxi].x - tri[minxi].x < 1)
    {
    }
    else
    {

        BYTE mediumyi = 3 - (minyi + maxyi);
        int mediumy = tri[mediumyi].y;
        const size_t RGBA = 4;
        bool achanged = false;
        int ax1 = 0;
        int ay1 = 0;
        float az1 = 0;
        coordinate texa1, texa2;
        if (tri[minyi].x < tri[mediumyi].x)
        {
            ax1 = tri[minyi].x;
            ay1 = tri[minyi].y;
            az1 = tri[minyi].z;
            texa1 = tritex[minyi];
        }
        else
        {
            ax1 = tri[mediumyi].x;
            ay1 = tri[mediumyi].y;
            az1 = tri[mediumyi].z;
            texa1 = tritex[mediumyi];
        }
        int ax2 = tri[maxyi].x;
        int ay2 = tri[maxyi].y;
        float az2 = tri[maxyi].z;
        texa2 = tritex[maxyi];
        int ax = ax1;
        int ay = ay1;
        float az = az1;
        float invaz = 1 / az;
        int recax = ax;
        int recay = ay;
        float recaz = az;
        double invaz1 = 1 / double(az1);
        double invaz2 = 1 / double(az2);
        coordinate intertexa1(texa1.x * invaz1, texa1.y * invaz1);
        coordinate intertexa2(texa2.x * invaz2, texa2.y * invaz2);
        coordinate intertexa = intertexa1;

        int adx = inmod(ax2 - ax1);
        int ady = inmod(ay2 - ay1);
        int asignx = signum(ax2 - ax1);
        int asigny = signum(ay2 - ay1);
        int recai = 1;
        if (ady > adx)
        {
            std::swap(adx, ady);
            achanged = true;
        }
        double invadx;
        if (achanged == false)
        {
            invadx = 1 / double(adx);
        }
        else
        {
            invadx = 1 / double(adx);
        }
        double aconst = invadx * (invaz2 - invaz1);
        double texaconstx = invadx * (intertexa2.x - intertexa1.x);
        double texaconsty = invadx * (intertexa2.y - intertexa1.y);
        float ae = (ady + ady) - adx;

        bool bchanged = false;
        int bx1 = 0;
        int by1 = 0;
        coordinate texb1;
        float bz1 = 0;
        if (tri[minyi].x < tri[mediumyi].x)
        {
            bx1 = tri[mediumyi].x;
            by1 = tri[mediumyi].y;
            bz1 = tri[mediumyi].z;
            texb1 = tritex[mediumyi];
        }
        else
        {
            bx1 = tri[minyi].x;
            by1 = tri[minyi].y;
            bz1 = tri[minyi].z;
            texb1 = tritex[minyi];
        }
        int bx2 = tri[maxyi].x;
        int by2 = tri[maxyi].y;
        float bz2 = tri[maxyi].z;
        coordinate texb2 = tritex[maxyi];
        int bx = bx1;
        int by = by1;
        float bz = bz1;
        double invbz1 = 1 / bz1;
        double invbz2 = 1 / bz2;
        double invbz = 1 / bz;
        coordinate intertexb1(texb1.x * invbz1, texb1.y * invbz1);
        coordinate intertexb2(texb2.x * invbz2, texb2.y * invbz2);
        coordinate intertexb = intertexb1;
        int bdx = inmod(bx2 - bx1);
        int bdy = inmod(by2 - by1);
        int bsignx = signum(bx2 - bx1);
        int bsigny = signum(by2 - by1);
        int recbi = 1;

        if (bdy > bdx)
        {
            std::swap(bdx, bdy);
            bchanged = true;
        }
        double invbdx;
        if (bchanged == false)
        {
            invbdx = 1 / double(bdx);
        }
        else
        {
            invbdx = 1 / double(bdx);
        }
        double bconst = invbdx * (invbz2 - invbz1);
        double texbconstx = invbdx * (intertexb2.x - intertexb1.x);
        double texbconsty = invbdx * (intertexb2.y - intertexb1.y);
        float be = (bdy + bdy) - bdx;
        float maximumy = 0;
        float minimumy = 0;
        if (miny < 0)
        {
            minimumy = 0;
        }
        else
        {
            minimumy = miny;
        }
        if (maxy > winheight - 1)
        {
            maximumy = winheight - 1;
        }
        else
        {
            maximumy = maxy;
        }
        for (int i = miny; i < minimumy; i++)
        {

            for (int ai = recai; ai <= adx; ai++)
            {
                if (ay == i)
                {
                    recai = ai;
                    break;
                }
                invaz = invaz + aconst;
                intertexa.x = intertexa.x + texaconstx;
                intertexa.y = intertexa.y + texaconsty;
                while (ae >= 0)
                {
                    if (achanged)
                    {
                        ax = ax + asignx;
                    }
                    else
                    {
                        ay = ay + asigny;
                    }
                    ae = ae - (adx + adx);
                }
                if (achanged)
                {
                    ay += asigny;
                }
                else
                {
                    ax += asignx;
                }
                ae = ae + (ady + ady);
            }

            for (int bi = recbi; bi <= bdx; bi++)
            {
                if (by == i)
                {
                    recbi = bi;
                    break;
                }
                invbz = invbz + bconst;
                intertexb.x = intertexb.x + texbconstx;
                intertexb.y = intertexb.y + texbconsty;
                while (be >= 0)
                {
                    if (bchanged)
                    {
                        bx = bx + bsignx;
                    }
                    else
                    {
                        by = by + bsigny;
                    }
                    be = be - (bdx + bdx);
                }
                if (bchanged)
                {
                    by += bsigny;
                }
                else
                {
                    bx += bsignx;
                }
                be = be + (bdy + bdy);
            }
            int tempax = ax;
            int tempbx = bx;
            /*
            if(tempax>winwidth-1)
            {
                tempax=winwidth-1;
            }
            if(tempax<0)
            {
                tempax=0;
            }
            if(tempbx>winwidth-1)
            {
                tempbx=winwidth-1;
            }
            if(tempbx<0)
            {
                tempbx=0;
            }
            */

            /*
             double tempinvdx=1;
             if(tempax!=tempbx)
             {
                 tempinvdx=1/double(tempbx-tempax);
             }
             double tempinvz=invaz;
             coordinate tempintertex=intertexa;
             coordinate temptex;
             double temptexconstx=tempinvdx*(intertexb.x-intertexa.x);
             double temptexconsty=tempinvdx*(intertexb.y-intertexa.y);
             double tempconst=tempinvdx*(invbz-invaz);
             double tempz;
             */

            /*
            if(i>=0 && i<=winheight-1)
            {
                float minimumx=0;
                float maximumx=0;
                if(tempax<0)
                {
                    minimumx=0;
                }
                else
                {
                    minimumx=tempax;
                }
                if(tempbx>winwidth-1)
                {
                    maximumx=winwidth-1;
                }
                else
                {
                    maximumx=tempbx;
                }
            for(int j=tempax; j<minimumx; j++)
            {

               if(j>=0 && j<=winwidth-1)
               {

                if(tempinvz>depthbuffer[i][j])
                {

                depthbuffer[i][j]=tempinvz;
                tempz=1/tempinvz;
                temptex.x=tempintertex.x*tempz;
                temptex.y=tempintertex.y*tempz;

                size_t index=(int(temptex.y*imgheight)*int(imgwidth) + int((temptex.x)*imgwidth));

                draw_pixel(j,i, image[index]);
                }
               }

                tempinvz=tempinvz+tempconst;
                tempintertex.x=tempintertex.x+temptexconstx;
                tempintertex.y=tempintertex.y+temptexconsty;

            }
            for(int j=minimumx; j<maximumx; j++)
            {

                if(tempinvz>depthbuffer[i][j])
                {

                depthbuffer[i][j]=tempinvz;
                tempz=1/tempinvz;
                temptex.x=tempintertex.x*tempz;
                temptex.y=tempintertex.y*tempz;

                size_t index=(int(temptex.y*imgheight)*int(imgwidth) + int((temptex.x)*imgwidth));

                draw_pixel(j,i,RGBtoHex( static_cast<int>(image[index+0]),static_cast<int>(image[index+1]), static_cast<int>(image[index+2]) ) image[index]);
                }

                tempinvz=tempinvz+tempconst;
                tempintertex.x=tempintertex.x+temptexconstx;
                tempintertex.y=tempintertex.y+temptexconsty;

            }
            }
            */
            recax = ax;
            recay = ay;
            recaz = az;
        }
        for (int i = minimumy; i <= maximumy; i++)
        {

            for (int ai = recai; ai <= adx; ai++)
            {
                if (ay == i)
                {
                    recai = ai;
                    break;
                }
                invaz = invaz + aconst;
                intertexa.x = intertexa.x + texaconstx;
                intertexa.y = intertexa.y + texaconsty;
                while (ae >= 0)
                {
                    if (achanged)
                    {
                        ax = ax + asignx;
                    }
                    else
                    {
                        ay = ay + asigny;
                    }
                    ae = ae - (adx + adx);
                }
                if (achanged)
                {
                    ay += asigny;
                }
                else
                {
                    ax += asignx;
                }
                ae = ae + (ady + ady);
            }

            for (int bi = recbi; bi <= bdx; bi++)
            {
                if (by == i)
                {
                    recbi = bi;
                    break;
                }
                invbz = invbz + bconst;
                intertexb.x = intertexb.x + texbconstx;
                intertexb.y = intertexb.y + texbconsty;
                while (be >= 0)
                {
                    if (bchanged)
                    {
                        bx = bx + bsignx;
                    }
                    else
                    {
                        by = by + bsigny;
                    }
                    be = be - (bdx + bdx);
                }
                if (bchanged)
                {
                    by += bsigny;
                }
                else
                {
                    bx += bsignx;
                }
                be = be + (bdy + bdy);
            }
            int tempax = ax;
            int tempbx = bx;
            /*
            if(tempax>winwidth-1)
            {
                tempax=winwidth-1;
            }
            if(tempax<0)
            {
                tempax=0;
            }
            if(tempbx>winwidth-1)
            {
                tempbx=winwidth-1;
            }
            if(tempbx<0)
            {
                tempbx=0;
            }
            */

            double tempinvdx = 1;
            if (tempax != tempbx)
            {
                tempinvdx = 1 / double(tempbx - tempax);
            }
            double tempinvz = invaz;
            coordinate tempintertex = intertexa;
            coordinate temptex;
            double temptexconstx = tempinvdx * (intertexb.x - intertexa.x);
            double temptexconsty = tempinvdx * (intertexb.y - intertexa.y);
            double tempconst = tempinvdx * (invbz - invaz);
            double tempz;

            float minimumx = 0;
            float maximumx = 0;
            if (tempax < 0)
            {
                minimumx = 0;
            }
            else
            {
                minimumx = tempax;
            }
            if (tempbx > winwidth - 1)
            {
                maximumx = winwidth - 1;
            }
            else
            {
                maximumx = tempbx;
            }
            /*
            if(j>=0 && j<=winwidth-1)
            {

             if(tempinvz>depthbuffer[i][j])
             {

             depthbuffer[i][j]=tempinvz;
             tempz=1/tempinvz;
             temptex.x=tempintertex.x*tempz;
             temptex.y=tempintertex.y*tempz;

             size_t index=(int(temptex.y*imgheight)*int(imgwidth) + int((temptex.x)*imgwidth));

             draw_pixel(j,i, image[index]);
             }
            }
            */
            tempinvz = tempinvz + (minimumx - float(tempax)) * tempconst;
            tempintertex.x = tempintertex.x + (minimumx - float(tempax)) * temptexconstx;
            tempintertex.y = tempintertex.y + (minimumx - float(tempax)) * temptexconsty;

            coordinate tempinstcoord(float(minimumx) * invwidth, float(i) * invheight);

            for (int j = minimumx; j < maximumx; j++)
            {

                if (tempinvz > depthbuffer[size_t(i)][size_t(j)])
                {

                    depthbuffer[size_t(i)][size_t(j)] = tempinvz;
                    tempz = 1 / tempinvz;
                    temptex.x = tempintertex.x * tempz;
                    temptex.y = tempintertex.y * tempz;

                    size_t index = RGBA * (int(temptex.y) * int(imgwidth) + int(temptex.x));

                    float speccoeff;
                    if (cosine > 0)
                    {
                        coordinate tempcoord(tempinstcoord.x + tempinstcoord.x - 1, tempinstcoord.y + tempinstcoord.y - 1);
                        vertex tempver(-tempcoord.x * uw * invud * tempz, -tempcoord.y * uh * invud * tempz, -tempz);
                        speccoeff = dotprod(tempver, reflection);
                        if (speccoeff > 0.0f)
                        {
                            speccoeff = (speccoeff * speccoeff * invreflect_mag) / (tempver.x * tempver.x + tempver.y * tempver.y + tempver.z * tempver.z);
                            speccoeff = speccoeff * speccoeff;
                            speccoeff = speccoeff * speccoeff;
                            speccoeff = speccoeff * speccoeff;
                            speccoeff = speccoeff * speccoeff;
                            speccoeff = speccoeff * speccoeff;
                        }
                        else
                        {
                            speccoeff = 0;
                        }
                    }
                    else
                    {
                        speccoeff = 0;
                    }
                    draw_pixel(j, i, RGBtoHex((1 - speccoeff) * coeff * static_cast<int>(image[index + 0]) + speccoeff * 255, (1 - speccoeff) * coeff * static_cast<int>(image[index + 1]) + speccoeff * 255, (1 - speccoeff) * coeff * static_cast<int>(image[index + 2]) + speccoeff * 255) /*coeff*float(image[index])*/);
                }

                tempinstcoord.x = tempinstcoord.x + invwidth;
                tempinvz = tempinvz + tempconst;
                tempintertex.x = tempintertex.x + temptexconstx;
                tempintertex.y = tempintertex.y + temptexconsty;
            }

            recax = ax;
            recay = ay;
            recaz = az;
        }
    }
}

void draw2dtriangledowntex(coordinate a, coordinate b, coordinate c, int winwidth, int winheight, coordinate uvtex1, coordinate uvtex2, coordinate uvtex3, vertex normal, vertex lightdir, float diffuse, float ambient, std::vector<unsigned char> &image, float imgwidth, float imgheight)
{
    float cosine = -dotprod(normal, lightdir);
    vertex reflection;
    reflection.x = cosine * 2.0f * normal.x + lightdir.x;
    reflection.y = cosine * 2.0f * normal.y + lightdir.y;
    reflection.z = cosine * 2.0f * normal.z + lightdir.z;
    float invreflect_mag = reflection.x * reflection.x + reflection.y * reflection.y + reflection.z * reflection.z;
    invreflect_mag = 1 / invreflect_mag;
    float invheight = 1 / float(height);
    float invwidth = 1 / float(width);
    if (cosine < 0)
    {
        cosine = 0;
    }
    /*
    float coeff=ambient+cosine*(diffuse-ambient);
    */
    float coeff = ambient + cosine * (diffuse - ambient);
    if (coeff > 1)
    {
        coeff = 1;
    }
    coordinate tex1(uvtex1.x * imgwidth, uvtex1.y * imgheight);
    coordinate tex2(uvtex2.x * imgwidth, uvtex2.y * imgheight);
    coordinate tex3(uvtex3.x * imgwidth, uvtex3.y * imgheight);
    coordinate tri[3] = {a, b, c};
    coordinate tritex[3] = {tex1, tex2, tex3};
    int maxy = tri[0].y;
    int miny = tri[0].y;
    BYTE maxyi = 0;
    BYTE minyi = 0;
    int maxx = tri[0].x;
    int minx = tri[0].x;
    BYTE maxxi = 0;
    BYTE minxi = 0;
    for (BYTE i = 1; i < 3; i++)
    {
        if (tri[i].y > maxy)
        {
            maxyi = i;
            maxy = tri[i].y;
        }
        if (tri[i].y < miny)
        {
            minyi = i;
            miny = tri[i].y;
        }
        if (tri[i].x > maxx)
        {
            maxxi = i;
            maxx = tri[i].x;
        }
        if (tri[i].x < minx)
        {
            minxi = i;
            minx = tri[i].x;
        }
    }
    if (a.x < 0 && b.x < 0 && c.x < 0 || a.x > winwidth && b.x > winwidth && c.x > winwidth || a.y < 0 && b.y < 0 && c.y < 0 || a.y > winheight && b.y > winheight && c.y > winheight || tri[maxyi].y - tri[minyi].y < 1 || tri[maxxi].x - tri[minxi].x < 1)
    {
    }
    else
    {

        BYTE mediumyi = 3 - (minyi + maxyi);
        int mediumy = tri[mediumyi].y;
        const size_t RGBA = 4;
        bool achanged = false;
        int ax1 = 0;
        int ay1 = 0;
        float az1 = 0;
        coordinate texa1;
        if (tri[maxyi].x < tri[mediumyi].x)
        {
            ax1 = tri[maxyi].x;
            ay1 = tri[maxyi].y;
            az1 = tri[maxyi].z;
            texa1 = tritex[maxyi];
        }
        else
        {
            ax1 = tri[mediumyi].x;
            ay1 = tri[mediumyi].y;
            az1 = tri[mediumyi].z;
            texa1 = tritex[mediumyi];
        }

        int ax2 = tri[minyi].x;
        int ay2 = tri[minyi].y;
        float az2 = tri[minyi].z;
        coordinate texa2 = tritex[minyi];

        int ax = ax1;
        int ay = ay1;
        float az = az1;
        double invaz = 1 / az;
        double invaz1 = 1 / az1;
        double invaz2 = 1 / az2;
        coordinate intertexa1(texa1.x * invaz1, texa1.y * invaz1);
        coordinate intertexa2(texa2.x * invaz2, texa2.y * invaz2);
        coordinate intertexa = intertexa1;
        int adx = inmod(ax2 - ax1);
        int ady = inmod(ay2 - ay1);
        int asignx = signum(ax2 - ax1);
        int asigny = signum(ay2 - ay1);
        int recai = 1;
        if (ady > adx)
        {
            std::swap(adx, ady);
            achanged = true;
        }

        double invadx;
        if (achanged == false)
        {
            invadx = 1 / double(adx);
        }
        else
        {
            invadx = 1 / double(adx);
        }
        double aconst = invadx * (invaz2 - invaz1);
        double texaconstx = invadx * (intertexa2.x - intertexa1.x);
        double texaconsty = invadx * (intertexa2.y - intertexa1.y);
        float ae = (ady + ady) - adx;

        bool bchanged = false;
        int bx1 = 0;
        int by1 = 0;
        float bz1 = 0;
        coordinate texb1;
        if (tri[maxyi].x < tri[mediumyi].x)
        {
            bx1 = tri[mediumyi].x;
            by1 = tri[mediumyi].y;
            bz1 = tri[mediumyi].z;
            texb1 = tritex[mediumyi];
        }
        else
        {
            bx1 = tri[maxyi].x;
            by1 = tri[maxyi].y;
            bz1 = tri[maxyi].z;
            texb1 = tritex[maxyi];
        }
        int bx2 = tri[minyi].x;
        int by2 = tri[minyi].y;
        float bz2 = tri[minyi].z;
        coordinate texb2 = tritex[minyi];
        int bx = bx1;
        int by = by1;
        float bz = bz1;
        double invbz = 1 / bz;
        double invbz1 = 1 / bz1;
        double invbz2 = 1 / bz2;
        coordinate intertexb1(texb1.x * invbz1, texb1.y * invbz1);
        coordinate intertexb2(texb2.x * invbz2, texb2.y * invbz2);
        coordinate intertexb = intertexb1;
        int bdx = inmod(bx2 - bx1);
        int bdy = inmod(by2 - by1);
        int bsignx = signum(bx2 - bx1);
        int bsigny = signum(by2 - by1);
        int recbi = 1;

        if (bdy > bdx)
        {
            std::swap(bdx, bdy);
            bchanged = true;
        }
        double invbdx = 0;
        if (bchanged == false)
        {
            invbdx = 1 / double(bdx);
        }
        else
        {
            invbdx = 1 / double(bdx);
        }
        double bconst = invbdx * (invbz2 - invbz1);
        double texbconstx = invbdx * (intertexb2.x - intertexb1.x);
        double texbconsty = invbdx * (intertexb2.y - intertexb1.y);
        float be = (bdy + bdy) - bdx;
        float maximumy = 0;
        float minimumy = 0;
        if (miny < 0)
        {
            minimumy = 0;
        }
        else
        {
            minimumy = miny;
        }
        if (maxy > winheight - 1)
        {
            maximumy = winheight - 1;
        }
        else
        {
            maximumy = maxy;
        }

        for (int i = maxy; i > maximumy; i--)
        {

            for (int ai = recai; ai <= adx; ai++)
            {
                if (ay == i)
                {
                    recai = ai;
                    break;
                }
                invaz = invaz + aconst;
                intertexa.x = intertexa.x + texaconstx;
                intertexa.y = intertexa.y + texaconsty;
                while (ae >= 0)
                {
                    if (achanged)
                    {
                        ax = ax + asignx;
                    }
                    else
                    {
                        ay = ay + asigny;
                    }
                    ae = ae - (adx + adx);
                }
                if (achanged)
                {
                    ay += asigny;
                }
                else
                {
                    ax += asignx;
                }
                ae = ae + (ady + ady);
            }

            for (int bi = recbi; bi <= bdx; bi++)
            {
                if (by == i)
                {
                    recbi = bi;
                    break;
                }
                invbz = invbz + bconst;
                intertexb.x = intertexb.x + texbconstx;
                intertexb.y = intertexb.y + texbconsty;
                while (be >= 0)
                {
                    if (bchanged)
                    {
                        bx = bx + bsignx;
                    }
                    else
                    {
                        by = by + bsigny;
                    }
                    be = be - (bdx + bdx);
                }
                if (bchanged)
                {
                    by += bsigny;
                }
                else
                {
                    bx += bsignx;
                }
                be = be + (bdy + bdy);
            }
            int tempax = ax;
            int tempbx = bx;
            /*
            if(tempax>winwidth-1)
            {
                tempax=winwidth-1;
            }
            if(tempax<0)
            {
                tempax=0;
            }
            if(tempbx>winwidth-1)
            {
                tempbx=winwidth-1;
            }
            if(tempbx<0)
            {
                tempbx=0;
            }
            */

            /*
            double tempinvdx=1;
            if(tempax!=tempbx)
            {
                tempinvdx=1/double(tempbx-tempax);
            }
            double tempinvz=invaz;
            coordinate tempintertex=intertexa;
            coordinate temptex;
            double temptexconstx=tempinvdx*(intertexb.x-intertexa.x);
            double temptexconsty=tempinvdx*(intertexb.y-intertexa.y);
            double tempconst=tempinvdx*(invbz-invaz);
            double tempz;
            */

            /*
            if(i>=0 && i<=winheight-1)
            {
            float minimumx=0;
                float maximumx=0;
                if(tempax<0)
                {
                    minimumx=0;
                }
                else
                {
                    minimumx=tempax;
                }
                if(tempbx>winwidth-1)
                {
                    maximumx=winwidth-1;
                }
                else
                {
                    maximumx=tempbx;
                }
            for(int j=tempax; j<minimumx; j++)
            {

               if(j>=0 && j<=winwidth-1)
               {

                if(tempinvz>depthbuffer[i][j])
                {

                depthbuffer[i][j]=tempinvz;
                tempz=1/tempinvz;
                temptex.x=tempintertex.x*tempz;
                temptex.y=tempintertex.y*tempz;

                size_t index=(int(temptex.y*imgheight)*int(imgwidth) + int((temptex.x)*imgwidth));

                draw_pixel(j,i, image[index]);
                }
               }

                tempinvz=tempinvz+tempconst;
                tempintertex.x=tempintertex.x+temptexconstx;
                tempintertex.y=tempintertex.y+temptexconsty;

            }
            for(int j=minimumx; j<maximumx; j++)
            {

                if(tempinvz>depthbuffer[i][j])
                {

                depthbuffer[i][j]=tempinvz;
                tempz=1/tempinvz;
                temptex.x=tempintertex.x*tempz;
                temptex.y=tempintertex.y*tempz;

                size_t index=(int(temptex.y*imgheight)*int(imgwidth) + int((temptex.x)*imgwidth));

                draw_pixel(j,i,RGBtoHex( static_cast<int>(image[index+0]),static_cast<int>(image[index+1]), static_cast<int>(image[index+2]) ) image[index]);
                }

                tempinvz=tempinvz+tempconst;
                tempintertex.x=tempintertex.x+temptexconstx;
                tempintertex.y=tempintertex.y+temptexconsty;

            }
            }
            */
        }
        for (int i = maximumy; i >= minimumy; i--)
        {

            for (int ai = recai; ai <= adx; ai++)
            {
                if (ay == i)
                {
                    recai = ai;
                    break;
                }
                invaz = invaz + aconst;
                intertexa.x = intertexa.x + texaconstx;
                intertexa.y = intertexa.y + texaconsty;
                while (ae >= 0)
                {
                    if (achanged)
                    {
                        ax = ax + asignx;
                    }
                    else
                    {
                        ay = ay + asigny;
                    }
                    ae = ae - (adx + adx);
                }
                if (achanged)
                {
                    ay += asigny;
                }
                else
                {
                    ax += asignx;
                }
                ae = ae + (ady + ady);
            }

            for (int bi = recbi; bi <= bdx; bi++)
            {
                if (by == i)
                {
                    recbi = bi;
                    break;
                }
                invbz = invbz + bconst;
                intertexb.x = intertexb.x + texbconstx;
                intertexb.y = intertexb.y + texbconsty;
                while (be >= 0)
                {
                    if (bchanged)
                    {
                        bx = bx + bsignx;
                    }
                    else
                    {
                        by = by + bsigny;
                    }
                    be = be - (bdx + bdx);
                }
                if (bchanged)
                {
                    by += bsigny;
                }
                else
                {
                    bx += bsignx;
                }
                be = be + (bdy + bdy);
            }
            int tempax = ax;
            int tempbx = bx;
            /*
            if(tempax>winwidth-1)
            {
                tempax=winwidth-1;
            }
            if(tempax<0)
            {
                tempax=0;
            }
            if(tempbx>winwidth-1)
            {
                tempbx=winwidth-1;
            }
            if(tempbx<0)
            {
                tempbx=0;
            }
            */
            double tempinvdx = 1;
            if (tempax != tempbx)
            {
                tempinvdx = 1 / double(tempbx - tempax);
            }
            double tempinvz = invaz;
            coordinate tempintertex = intertexa;
            coordinate temptex;
            double temptexconstx = tempinvdx * (intertexb.x - intertexa.x);
            double temptexconsty = tempinvdx * (intertexb.y - intertexa.y);
            double tempconst = tempinvdx * (invbz - invaz);
            double tempz;

            float minimumx = 0;
            float maximumx = 0;
            if (tempax < 0)
            {
                minimumx = 0;
            }
            else
            {
                minimumx = tempax;
            }
            if (tempbx > winwidth - 1)
            {
                maximumx = winwidth - 1;
            }
            else
            {
                maximumx = tempbx;
            }
            /*
            if(j>=0 && j<=winwidth-1)
            {

             if(tempinvz>depthbuffer[i][j])
             {

             depthbuffer[i][j]=tempinvz;
             tempz=1/tempinvz;
             temptex.x=tempintertex.x*tempz;
             temptex.y=tempintertex.y*tempz;

             size_t index=(int(temptex.y*imgheight)*int(imgwidth) + int((temptex.x)*imgwidth));

             draw_pixel(j,i, image[index]);
             }
            }
            */
            tempinvz = tempinvz + (minimumx - float(tempax)) * tempconst;
            tempintertex.x = tempintertex.x + (minimumx - float(tempax)) * temptexconstx;
            tempintertex.y = tempintertex.y + (minimumx - float(tempax)) * temptexconsty;

            coordinate tempinstcoord(float(minimumx) * invwidth, float(i) * invheight);
            for (int j = minimumx; j < maximumx; j++)
            {

                if (tempinvz > depthbuffer[size_t(i)][size_t(j)])
                {

                    depthbuffer[size_t(i)][size_t(j)] = tempinvz;
                    tempz = 1 / tempinvz;
                    temptex.x = tempintertex.x * tempz;
                    temptex.y = tempintertex.y * tempz;

                    size_t index = RGBA * (int(temptex.y) * int(imgwidth) + int(temptex.x));

                    float speccoeff;
                    if (cosine > 0)
                    {
                        coordinate tempcoord(tempinstcoord.x + tempinstcoord.x - 1, tempinstcoord.y + tempinstcoord.y - 1);
                        vertex tempver(-tempcoord.x * uw * invud * tempz, -tempcoord.y * uh * invud * tempz, -tempz);
                        speccoeff = dotprod(tempver, reflection);
                        if (speccoeff > 0.0f)
                        {
                            speccoeff = (speccoeff * speccoeff * invreflect_mag) / (tempver.x * tempver.x + tempver.y * tempver.y + tempver.z * tempver.z);
                            speccoeff = speccoeff * speccoeff;
                            speccoeff = speccoeff * speccoeff;
                            speccoeff = speccoeff * speccoeff;
                            speccoeff = speccoeff * speccoeff;
                            speccoeff = speccoeff * speccoeff;
                        }
                        else
                        {
                            speccoeff = 0;
                        }
                    }
                    else
                    {
                        speccoeff = 0;
                    }
                    short difR = coeff * static_cast<int>(image[index + 0]);
                    short difG = coeff * static_cast<int>(image[index + 1]);
                    short difB = coeff * static_cast<int>(image[index + 2]);
                    draw_pixel(j, i, RGBtoHex(difR + speccoeff * (255 - difR), difG + speccoeff * (255 - difG), difB + speccoeff * (255 - difB)) /*coeff*float(image[index])*/);
                }

                tempinstcoord.x = tempinstcoord.x + invwidth;
                tempinvz = tempinvz + tempconst;
                tempintertex.x = tempintertex.x + temptexconstx;
                tempintertex.y = tempintertex.y + temptexconsty;
            }
        }
    }
}

void draw2dtriangletex(coordinate a, coordinate b, coordinate c, int winwidth, int winheight, coordinate tex1, coordinate tex2, coordinate tex3, vertex normal, vertex lightdir, float diffuse, float ambient, std::vector<unsigned char> &image, float imgwidth, float imgheight)
{
    if (a.x < 0 && b.x < 0 && c.x < 0 || a.x > winwidth && b.x > winwidth && c.x > winwidth || a.y < 0 && b.y < 0 && c.y < 0 || a.y > winheight && b.y > winheight && c.y > winheight)
    {
    }
    else
    {
        coordinate tri[3] = {a, b, c};
        coordinate tritex[3] = {tex1, tex2, tex3};
        float maxy = tri[0].y;
        float miny = tri[0].y;
        BYTE maxyi = 0;
        BYTE minyi = 0;
        float maxx = tri[0].x;
        float minx = tri[0].x;
        BYTE maxxi = 0;
        BYTE minxi = 0;
        for (BYTE i = 1; i < 3; i++)
        {
            if (tri[i].y > maxy)
            {
                maxyi = i;
                maxy = tri[i].y;
            }
            if (tri[i].y < miny)
            {
                minyi = i;
                miny = tri[i].y;
            }
            if (tri[i].x > maxx)
            {
                maxxi = i;
                maxx = tri[i].x;
            }
            if (tri[i].x < minx)
            {
                minxi = i;
                minx = tri[i].x;
            }
        }

        BYTE mediumyi = 3 - (minyi + maxyi);
        float mediumy = tri[mediumyi].y;
        coordinate d;
        d.y = tri[mediumyi].y;
        d.x = tri[maxyi].x + (tri[minyi].x - tri[maxyi].x) * ((tri[mediumyi].y - tri[maxyi].y) / (tri[minyi].y - tri[maxyi].y));

        coordinate texd;
        float alpha;
        if (inmod(tri[maxyi].x - tri[minyi].x) > inmod(tri[maxyi].y - tri[minyi].y))
        {
            alpha = (d.x - tri[minyi].x) / (tri[maxyi].x - tri[minyi].x);
        }
        else
        {
            alpha = (d.y - tri[minyi].y) / (tri[maxyi].y - tri[minyi].y);
        }
        d.z = alpha * (1 / tri[maxyi].z) + (1 - alpha) * (1 / tri[minyi].z);
        d.z = 1 / d.z;
        coordinate temptexmax = coordinate(tritex[maxyi].x / tri[maxyi].z, tritex[maxyi].y / tri[maxyi].z);
        coordinate temptexmin = coordinate(tritex[minyi].x / tri[minyi].z, tritex[minyi].y / tri[minyi].z);
        texd = coordinate(temptexmin, temptexmax, alpha);
        texd.x = texd.x * d.z;
        texd.y = texd.y * d.z;

        draw2dtriangleuptex(tri[maxyi], tri[mediumyi], d, winwidth, winheight, tritex[maxyi], tritex[mediumyi], texd, normal, lightdir, diffuse, ambient, image, imgwidth, imgheight);

        draw2dtriangledowntex(tri[minyi], tri[mediumyi], d, winwidth, winheight, tritex[minyi], tritex[mediumyi], texd, normal, lightdir, diffuse, ambient, image, imgwidth, imgheight);
    }
}

void draw3dtriangletex(vertex a2, vertex b2, vertex c2, int winwidth, int winheight, coordinate tex1, coordinate tex2, coordinate tex3, vertex tempnormal, vertex templightdir, float diffuse, float ambient, std::vector<unsigned char> &image, float imgwidth, float imgheight, float area)
{

    vertex tempa = transformmat(vertex(a2.x - campos.x, a2.y - campos.y, a2.z - campos.z));
    vertex tempb = transformmat(vertex(b2.x - campos.x, b2.y - campos.y, b2.z - campos.z));
    vertex tempc = transformmat(vertex(c2.x - campos.x, c2.y - campos.y, c2.z - campos.z));
    vertex lightdir = transformmat(templightdir);
    vertex normal = crossprod(vertex(tempc, tempb), vertex(tempa, tempb));
    float invarea = 1 / area;
    normal = vertex(normal.x * invarea, normal.y * invarea, normal.z * invarea);
    /*
    vertex tempa=matrirxmulti(vertex(a2.x-campos.x,a2.y-campos.y,a2.z-campos.z),matrix);
    vertex tempb=matrirxmulti(vertex(b2.x-campos.x,b2.y-campos.y,b2.z-campos.z),matrix);
    vertex tempc=matrirxmulti(vertex(c2.x-campos.x,c2.y-campos.y,c2.z-campos.z),matrix);
    */

    coordinate a2d, b2d, c2d;
    vertex posa, posb, posc;
    float adistx, adisty;
    float bdistx, bdisty;
    float cdistx, cdisty;
    float adist = tempa.z;
    float bdist = tempb.z;
    float cdist = tempc.z;

    if (adist <= uznear || bdist <= uznear || cdist <= uznear)
    {

        if (adist <= uznear && bdist <= uznear && cdist <= uznear || adist <= 0 && bdist <= 0 && cdist <= 0)
        {
        }
        else
        {
            vertex tri[3] = {tempa, tempb, tempc};
            coordinate tritex[3] = {tex1, tex2, tex3};
            BYTE minzi = 0;
            float minz = tempa.z;
            BYTE maxzi = 0;
            float maxz = tempa.z;
            for (BYTE i = 1; i < 3; i++)
            {
                if (tri[i].z < minz)
                {
                    minz = tri[i].z;
                    minzi = i;
                }
                if (tri[i].z > maxz)
                {
                    maxz = tri[i].z;
                    maxzi = i;
                }
            }
            BYTE mediumzi = 3 - minzi - maxzi;
            float mediumz = tri[mediumzi].z;
            if (maxz > uznear)
            {

                if (mediumz > uznear)
                {
                    coordinate min2d, medium2d, max2d;
                    float invmediumz = 1 / mediumz;
                    float invmaxz = 1 / maxz;
                    vertex tempposmed = vertex(tri[mediumzi].x * ud * invmediumz, tri[mediumzi].y * ud * invmediumz, ud);
                    vertex tempposmax = vertex(tri[maxzi].x * ud * invmaxz, tri[maxzi].y * ud * invmaxz, ud);
                    float meddistx = tempposmed.x * invuw;
                    float meddisty = tempposmed.y * invuh;
                    float maxdistx = tempposmax.x * invuw;
                    float maxdisty = tempposmax.y * invuh;
                    medium2d.x = winwidth * (meddistx + 1) * 0.5;
                    medium2d.y = winheight * (meddisty + 1) * 0.5;
                    max2d.x = winwidth * (maxdistx + 1) * 0.5;
                    max2d.y = winheight * (maxdisty + 1) * 0.5;
                    medium2d.z = mediumz;
                    max2d.z = maxz;
                    coordinate texmax = tritex[maxzi];
                    coordinate texmed = tritex[mediumzi];

                    float lamdamax = rayplanelamda(0, 0, 1, vertex(0, 0, uznear), tri[maxzi], tri[minzi]);
                    float lamdamed = rayplanelamda(0, 0, 1, vertex(0, 0, uznear), tri[mediumzi], tri[minzi]);

                    vertex maxnear = vertex(tri[maxzi], tri[minzi], lamdamax);
                    coordinate nearmax2d;
                    float nearmax2ddistx = maxnear.x * invnearuw;
                    float nearmax2ddisty = maxnear.y * invnearuh;
                    nearmax2d.x = winwidth * (nearmax2ddistx + 1) * 0.5;
                    nearmax2d.y = winheight * (nearmax2ddisty + 1) * 0.5;
                    nearmax2d.z = maxnear.z;

                    coordinate texnearmax = coordinate(tritex[maxzi], tritex[minzi], lamdamax);

                    vertex mednear = vertex(tri[mediumzi], tri[minzi], lamdamed);
                    coordinate nearmed2d;
                    float nearmed2ddistx = mednear.x * invnearuw;
                    float nearmed2ddisty = mednear.y * invnearuh;
                    nearmed2d.x = winwidth * (nearmed2ddistx + 1) * 0.5;
                    nearmed2d.y = winheight * (nearmed2ddisty + 1) * 0.5;
                    nearmed2d.z = mednear.z;

                    coordinate texnearmed = coordinate(tritex[mediumzi], tritex[minzi], lamdamed);

                    draw2dtriangletex(nearmed2d, medium2d, max2d, winwidth, winheight, texnearmed, texmed, texmax, normal, lightdir, diffuse, ambient, image, imgwidth, imgheight);
                    draw2dtriangletex(nearmed2d, nearmax2d, max2d, winwidth, winheight, texnearmed, texnearmax, texmax, normal, lightdir, diffuse, ambient, image, imgwidth, imgheight);
                }
                else
                {
                    coordinate max2d;
                    float invmaxz = 1 / maxz;
                    vertex tempposmax = vertex(tri[maxzi].x * ud * invmaxz, tri[maxzi].y * ud * invmaxz, ud);
                    float maxdistx = tempposmax.x * invuw;
                    float maxdisty = tempposmax.y * invuh;
                    max2d.x = winwidth * (maxdistx + 1) * 0.5;
                    max2d.y = winheight * (maxdisty + 1) * 0.5;
                    max2d.z = maxz;
                    coordinate texmax = tritex[maxzi];

                    float lamdamed = rayplanelamda(0, 0, 1, vertex(0, 0, uznear), tri[maxzi], tri[mediumzi]);
                    float lamdamin = rayplanelamda(0, 0, 1, vertex(0, 0, uznear), tri[maxzi], tri[minzi]);

                    vertex mednear = vertex(tri[maxzi], tri[mediumzi], lamdamed);
                    coordinate nearmed2d;
                    float nearmed2ddistx = mednear.x * invnearuw;
                    float nearmed2ddisty = mednear.y * invnearuh;
                    nearmed2d.x = winwidth * (nearmed2ddistx + 1) * 0.5;
                    nearmed2d.y = winheight * (nearmed2ddisty + 1) * 0.5;
                    nearmed2d.z = mednear.z;

                    coordinate texnearmed = coordinate(tritex[maxzi], tritex[mediumzi], lamdamed);

                    vertex minnear = vertex(tri[maxzi], tri[minzi], lamdamin);
                    coordinate nearmin2d;
                    float nearmin2ddistx = minnear.x * invnearuw;
                    float nearmin2ddisty = minnear.y * invnearuh;
                    nearmin2d.x = winwidth * (nearmin2ddistx + 1) * 0.5;
                    nearmin2d.y = winheight * (nearmin2ddisty + 1) * 0.5;
                    nearmin2d.z = minnear.z;

                    coordinate texnearmin = coordinate(tritex[maxzi], tritex[minzi], lamdamin);

                    draw2dtriangletex(max2d, nearmin2d, nearmed2d, winwidth, winheight, texmax, texnearmin, texnearmed, normal, lightdir, diffuse, ambient, image, imgwidth, imgheight);
                }
            }
            else
            {
            }
        }
    }
    else
    {

        /*
        posa=rayplane(lx,ly, lz, camzpoint,campos, vertex(a.x,a.y,a.z));
        posb=rayplane(lx,ly, lz, camzpoint,campos, vertex(b.x,b.y,b.z));
        posc=rayplane(lx,ly, lz, camzpoint,campos, vertex(c.x,c.y,c.z));
        */
        float invadist = 1 / adist;
        float invbdist = 1 / bdist;
        float invcdist = 1 / cdist;
        vertex tempposa = vertex(tempa.x * ud * invadist, tempa.y * ud * invadist, ud);
        vertex tempposb = vertex(tempb.x * ud * invbdist, tempb.y * ud * invbdist, ud);
        vertex tempposc = vertex(tempc.x * ud * invcdist, tempc.y * ud * invcdist, ud);

        adistx = tempposa.x * invuw;
        adisty = tempposa.y * invuh;
        bdistx = tempposb.x * invuw;
        bdisty = tempposb.y * invuh;
        cdistx = tempposc.x * invuw;
        cdisty = tempposc.y * invuh;
        a2d.x = winwidth * (adistx + 1) * 0.5;
        a2d.y = winheight * (adisty + 1) * 0.5;
        b2d.x = winwidth * (bdistx + 1) * 0.5;
        b2d.y = winheight * (bdisty + 1) * 0.5;
        c2d.x = winwidth * (cdistx + 1) * 0.5;
        c2d.y = winheight * (cdisty + 1) * 0.5;
        a2d.z = tempa.z;
        b2d.z = tempb.z;
        c2d.z = tempc.z;

        draw2dtriangletex(a2d, b2d, c2d, winwidth, winheight, tex1, tex2, tex3, normal, lightdir, diffuse, ambient, image, imgwidth, imgheight);
    }
}

class object
{
public:
    string line;
    string trash;
    float pubobjsize;
    string pubtexturename;
    string texture, rectexture;
    string recmonotexture;
    string recfilename;
    float recscale;
    float scale;
    int rendermaterial = 0;
    int materialpos = 0, recpos = 0;
    int materialpos2 = 0, recpos2 = 0;
    int vcount = 0;
    int vncount = 0;
    int vtcount = 0;
    int fcount = 0;
    int mtlcount = 0;
    bool readonce = false;
    bool readonce2 = false;
    vertex objv[50000];
    normal objvn[50000];
    texturecoord objvt[50000];
    faces objf[50000];
    /*
    material objmtl[700];
    */

    float centerx = 0;
    float centery = 0;
    float centerz = 0;

    float scalecenterx = 0;
    float scalecentery = 0;
    float scalecenterz = 0;

    void LoadObj(string filename)
    {
        if (recfilename != filename)
        {
            ifstream myfile(filename.c_str());
            if (myfile.is_open())
            {
                while (getline(myfile, line))
                {
                    istringstream ss(line);
                    if (searchstring(line, "v ", 0))
                    {
                        vcount++;
                        ss >> trash >> objv[vcount].x >> objv[vcount].y >> objv[vcount].z;
                        centerx += objv[vcount].x;
                        centery += objv[vcount].y;
                        centerz += objv[vcount].z;
                    }
                    else if (searchstring(line, "vt ", 0))
                    {
                        vtcount++;
                        ss >> trash >> objvt[vtcount].s >> objvt[vtcount].t;
                    }
                    else if (searchstring(line, "vn ", 0))
                    {
                        vncount++;
                        ss >> trash >> objvn[vncount].x >> objvn[vncount].y >> objvn[vncount].z;
                    }
                    else if (searchstring(line, "f ", 0))
                    {
                        replacecharacter(line, '/', ' ');
                        istringstream ss2(line);
                        ss2 >> trash >> objf[fcount].v1 >> objf[fcount].vt1 >> objf[fcount].vn1 >> objf[fcount].v2 >> objf[fcount].vt2 >> objf[fcount].vn2 >> objf[fcount].v3 >> objf[fcount].vt3 >> objf[fcount].vn3;
                        objf[fcount].area = crossprodtrianglearea(objv[objf[fcount].v1], objv[objf[fcount].v2], objv[objf[fcount].v3]);
                        /**
                        objv[objf[fcount].v1].poly[objv[objf[fcount].v1].polycount]=fcount;
                        objv[objf[fcount].v1].polycount+=1;
                        objv[objf[fcount].v2].poly[objv[objf[fcount].v2].polycount]=fcount;
                        objv[objf[fcount].v2].polycount+=1;
                        objv[objf[fcount].v3].poly[objv[objf[fcount].v3].polycount]=fcount;
                        objv[objf[fcount].v3].polycount+=1;
                        **/
                        fcount++;
                    }
                    /*
                    else if(searchstring(line,"usemtl ",0))
                    {
                      ss >> trash >> objmtl[mtlcount].texturename;
                      objmtl[mtlcount].texturename=addstring(objmtl[mtlcount].texturename,".jpg");
                      objmtl[mtlcount].start=fcount;
                      objf[fcount].changetex=true;
                        mtlcount++;
                    }
                    */
                }
                /*
                for(int i=0; i<mtlcount; i++)
                {
                    objmtl[i].tex_2d = SOIL_load_OGL_texture
                    (
                        objmtl[i].texturename.c_str(),
                        SOIL_LOAD_AUTO,
                        SOIL_CREATE_NEW_ID,
                        SOIL_FLAG_MIPMAPS | SOIL_FLAG_INVERT_Y | SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT
                    );
                }
                */

                centerx = centerx / vcount;
                centery = centery / vcount;
                centerz = centerz / vcount;

                myfile.close();
            }
            else
            {
                cout << "Unable to open " << filename << endl;
            }
            recfilename = filename;
        }
    }

    void RenderObj(float objsize)
    {
        if (scale != objsize)
        {
            scalecenterx = centerx * objsize;
            scalecentery = centery * objsize;
            scalecenterz = centerz * objsize;
            scale = objsize;
        }
        for (int i = 0; i < fcount; i++)
        {
            rectexture = texture;
            recpos = materialpos;
            recpos2 = materialpos2;
            /*
            for(int j=0; j<mtlcount; j++)
            {
                if(objmtl[j].start==i)
                {

                    texture=objmtl[j].texturename;

                    materialpos=j;
                }
                else if(objmtl[j].start==i+1)
                {
                    materialpos2=j;
                }
            }
            */

            /*
                if(recpos!=materialpos || rectexture!=texture)
                {

                    glEnable(GL_TEXTURE_2D);
                    glBindTexture(GL_TEXTURE_2D,objmtl[materialpos].tex_2d);
                    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
                    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );


                }
            */
            /*
            glNormal3f(objvn[objf[i].vn1].x,objvn[objf[i].vn1].y, objvn[objf[i].vn1].z);
            glTexCoord2d(objvt[objf[i].vt1].s,objvt[objf[i].vt1].t);
            glVertex3f(objv[objf[i].v1].x*objsize,objv[objf[i].v1].y*objsize,objv[objf[i].v1].z*objsize);

            glNormal3f(objvn[objf[i].vn2].x,objvn[objf[i].vn2].y, objvn[objf[i].vn2].z);
            glTexCoord2d(objvt[objf[i].vt2].s,objvt[objf[i].vt2].t);
            glVertex3f(objv[objf[i].v2].x*objsize,objv[objf[i].v2].y*objsize,objv[objf[i].v2].z*objsize);

            glNormal3f(objvn[objf[i].vn3].x,objvn[objf[i].vn3].y, objvn[objf[i].vn3].z);
            glTexCoord2d(objvt[objf[i].vt3].s,objvt[objf[i].vt3].t);
            glVertex3f(objv[objf[i].v3].x*objsize,objv[objf[i].v3].y*objsize,objv[objf[i].v3].z*objsize);
            */
            draw3dtriangletex(vertex(objv[objf[i].v1], objsize), vertex(objv[objf[i].v2], objsize), vertex(objv[objf[i].v3], objsize), width, height, coordinate(objvt[objf[i].vt1].s, objvt[objf[i].vt1].t), coordinate(objvt[objf[i].vt2].s, objvt[objf[i].vt2].t), coordinate(objvt[objf[i].vt3].s, objvt[objf[i].vt3].t), vertex(objvn[objf[i].vn1].x, objvn[objf[i].vn1].y, objvn[objf[i].vn1].z), pointdir, 0.8, 0.3, image1.image, image1.imagewidth, image1.imageheight, objf[i].area * objsize);

            /*
            if(recpos2!=materialpos2 || i==fcount-1)
            {
                glEnd();

                glDisable(GL_TEXTURE_2D);

            }
            */
        }
        recpos = 0;
        materialpos = 0;
        recpos2 = 0;
        materialpos2 = 0;
        texture = "";
    }
};
object obj1;

void clear_screen(u32 color)
{
    u32 *pixel = (u32 *)memory;
    for (int i = 0;
         i < client_width * client_height;
         ++i)
    {
        *pixel++ = color;
    }
}

LRESULT CALLBACK
window_proc(HWND window,
            UINT message,
            WPARAM w_param,
            LPARAM l_param)
{
    switch (message)
    {

    case WM_KEYDOWN:
    {
        switch (w_param)
        {
        // "o" exits the program
        case 'O':
        {
            DestroyWindow(window);
        };
        }
    }
    break;

    case WM_DESTROY:
    {
        PostQuitMessage(0);
    }
    break;

    default:
    {
        return DefWindowProc(window,
                             message,
                             w_param,
                             l_param);
    }
    }

    return 0;
}
float mov;
long long last = 0;
long long frequency = 0;
int WINAPI
WinMain(HINSTANCE instance,
        HINSTANCE prev_instance,
        LPSTR cmd_line,
        int cmd_show)
{
    // window creation

    WNDCLASS window_class = {};

    const wchar_t class_name[] = L"MyWindowClass";

    window_class.lpfnWndProc = window_proc;
    window_class.hInstance = instance;
    window_class.lpszClassName = class_name;
    window_class.hCursor = LoadCursor(0, IDC_CROSS);

    if (!RegisterClass(&window_class))
    {
        MessageBox(0, L"RegisterClass failed", 0, 0);
        return GetLastError();
    }

    HWND window = CreateWindowEx(0,
                                 class_name,
                                 L"Window",
                                 WS_OVERLAPPEDWINDOW | WS_VISIBLE,
                                 CW_USEDEFAULT,
                                 CW_USEDEFAULT,
                                 CW_USEDEFAULT,
                                 CW_USEDEFAULT,
                                 0,
                                 0,
                                 instance,
                                 0);

    if (!window)
    {
        MessageBox(0, L"CreateWindowEx failed", 0, 0);
        return GetLastError();
    }

    // allocate memory

    RECT rect;
    GetClientRect(window, &rect);
    client_width = rect.right - rect.left;
    client_height = rect.bottom - rect.top;
    height = client_height;
    width = client_width;
    memory = VirtualAlloc(0,
                          client_width * client_height * 4,
                          MEM_RESERVE | MEM_COMMIT,
                          PAGE_READWRITE);
    // create BITMAPINFO struct for StretchDIBits

    BITMAPINFO bitmap_info;
    bitmap_info.bmiHeader.biSize = sizeof(bitmap_info.bmiHeader);
    bitmap_info.bmiHeader.biWidth = client_width;
    bitmap_info.bmiHeader.biHeight = client_height;
    bitmap_info.bmiHeader.biPlanes = 1;
    bitmap_info.bmiHeader.biBitCount = 32;
    bitmap_info.bmiHeader.biCompression = BI_RGB;

    HDC hdc = GetDC(window);

    // loop

    bool running = true;

    while (running)
    {
        MSG msg;
        while (PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE))
        {
            if (msg.message == WM_QUIT)
                running = false;
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
        RECT rect;
        GetClientRect(window, &rect);
        int temp_client_width = rect.right - rect.left;
        int temp_client_height = rect.bottom - rect.top;
        if (temp_client_height != height || temp_client_width != width)
        {
            running = false;
        }
        // clear screen with gray color
        clear_screen(0x333333);
        depthbuffer = std::vector<std::vector<float>>(height);
        for (int i = 0; i < height; i++)
        {
            depthbuffer[i] = std::vector<float>(width);
            for (int j = 0; j < width; j++)
            {
                depthbuffer[i][j] = /*INFINITE*/ 0;
            }
        }
        // draw white pixel at 100, 100 (from bottom left)
        long long current = PerformanceCounter();
        long long frequency = PerformanceFrequency();
        if (last == 0)
        {
            last = current;
        }
        float dt = float((current - last)) / float(frequency);
        last = current;
        /**
        currentTime = GetTickCount();
        if(init==true)
        {
            lasttime=currentTime;
            init=false;
        }
        float dt=(currentTime-lasttime);
        dt=dt/1000;
        **/
        lasttime = currentTime;
        GetCursorPos(&P);

        angley += (P.x - xOrigin) * sensitivity * dt;
        xOrigin = P.x;
        anglex += (P.y - yOrigin) * sensitivity * dt;
        yOrigin = P.y;
        ShowCursor(false);
        if (angley >= 360 || angley <= -360)
        {
            angley = 0;
        }
        if (anglex > 90)
        {
            anglex = 90;
        }
        if (anglex < -90)
        {
            anglex = -90;
        }

        if (P.x >= 1200)
        {
            SetCursorPos(20, P.y);
            xOrigin = 20;
        }
        else if (P.x <= 10)
        {
            SetCursorPos(1000, P.y);
            xOrigin = 1000;
        }
        if (P.y >= 700)
        {
            SetCursorPos(P.x, 20);
            yOrigin = 20;
        }
        else if (P.y <= 10)
        {
            SetCursorPos(P.x, 600);
            yOrigin = 600;
        }

        /*
        float targetdist=sqrt(campos.x*campos.x+campos.z*campos.z);
        if(targetdist!=0)
        {
        float sinangle=asin((-campos.x/targetdist))*180/PI;
        if(sinangle<0)
        {
        angley=-acos((-campos.z/targetdist))*180/PI;
        }
        else
        {
        angley=acos((-campos.z/targetdist))*180/PI;
        }
        }
        */

        if (moving)
        {
            /*
            campos.x+=movespeed*sin(angley*PI/180)*cos(anglex*PI/180)*dt;
            campos.y+=movespeed*-sin(anglex*PI/180)*dt;
            campos.z+=movespeed*cos(angley*PI/180)*cos(anglex*PI/180)*dt;
            */
            campos.x += movespeed * sin(angley * PI / 180) * dt;
            campos.z += movespeed * cos(angley * PI / 180) * dt;
        }
        if (backmoving)
        {
            /*
            campos.x-=movespeed*sin(angley*PI/180)*cos(anglex*PI/180)*dt;
            campos.y-=movespeed*-sin(anglex*PI/180)*dt;
            campos.z-=movespeed*cos(angley*PI/180)*cos(anglex*PI/180)*dt;
            */
            campos.x -= movespeed * sin(angley * PI / 180) * dt;
            campos.z -= movespeed * cos(angley * PI / 180) * dt;
        }

        if (rightmov)
        {

            campos.x += movespeed * cos(angley * PI / 180) * dt;
            campos.z -= movespeed * sin(angley * PI / 180) * dt;
        }

        if (leftmov)
        {
            campos.x -= movespeed * cos(angley * PI / 180) * dt;
            campos.z += movespeed * sin(angley * PI / 180) * dt;
        }

        if (!readonce)
        {
            image1.filename = "texture2.jpg";
            image1.success = load_image_lin(image1.image, image1.filename, image1.imagewidth, image1.imageheight);
            if (!image1.success)
            {
                std::cout << "Error loading image\n";
                return 1;
            }
            image2.filename = "chess_tex.jpg";
            image2.success = load_image_lin(image2.image, image2.filename, image2.imagewidth, image2.imageheight);
            if (!image2.success)
            {
                std::cout << "Error loading image\n";
                return 1;
            }
            obj1.LoadObj("cube.obj");
            readonce = true;
        }

        setfrustrum(campos, anglex, angley, anglez, fov, 0.11f, float(width) / float(height));
        /*
        draw3dtriangletex(vertex(0,0,5),vertex(0,0,1),vertex(0,4,5),width,height, coordinate(0,1.0f),coordinate(1.0f,1.0f),coordinate(0,0),image2.image,image2.imagewidth, image2.imageheight);
        */
        obj1.RenderObj(2.01f);
        /*
        draw3dtriangle(vertex(-2,0,3),vertex(2,0,3),vertex(-2,4,3),width,height, 0x00ff00);
        */

        /*
        draw3dtriangletex(vertex(4.3*2,0,2.5*2),vertex(4.3*2,3*2,2.5*2),vertex(5.3*2,0,0.8*2),width,height, coordinate(0,0.99f),coordinate(0.99f,0.99f),coordinate(0,0),image1.image,image1.imagewidth, image1.imageheight);
        */

        /*
        draw2dtriangle3(coordinate(100,100),coordinate(300,20),coordinate(200,300),width,height);
        */

        BresenhamInt(width / 2, 10 + (height / 2), width / 2, -10 + (height / 2));
        BresenhamInt(-10 + (width / 2), (height / 2), 10 + (width / 2), (height / 2));

        if (GetKeyState(87) & 0x8000)
        {
            moving = true;
        }
        else
        {
            moving = false;
        }
        if (GetKeyState(83) & 0x8000)
        {
            backmoving = true;
        }
        else
        {
            backmoving = false;
        }
        if (GetKeyState(68) & 0x8000)
        {
            rightmov = true;
        }
        else
        {
            rightmov = false;
        }
        if (GetKeyState(65) & 0x8000)
        {
            leftmov = true;
        }
        else
        {
            leftmov = false;
        }
        if (GetKeyState(VK_SHIFT) & 0x8000)
        {
            movespeed = 10.0f;
        }
        else
        {
            movespeed = 5.0f;
        }
        for (int i = 0; i < height; i++)
        {
            depthbuffer[i].clear();
        }
        depthbuffer.clear();
        StretchDIBits(hdc,
                      0,
                      0,
                      client_width,
                      client_height,
                      0,
                      0,
                      client_width,
                      client_height,
                      memory,
                      &bitmap_info,
                      DIB_RGB_COLORS,
                      SRCCOPY);
    }

    return 0;
}
