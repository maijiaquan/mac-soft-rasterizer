#pragma once

#include "predefine.h"
#include "ds.h"

#include <iostream>
using namespace std;

extern int display_w, display_h;
extern int *framebuffer; 

//绘制像素点
void device_pixel(int *fb, const int &x, const int &y, const int &c);
void device_pixel(int *fb, const int &x, const int &y, const int &r, const int &g, const int &b);

//绘制一条纯色线
void device_draw_line(int *fb, int x1, int y1, int x2, int y2, int c); 

//绘制线框三角形
void DrawTriangleWireframe(int *fb, float x1, float y1, float x2, float y2, float x3, float y3, int color);

//绘制纯色TopTriangle
void DrawTopTrianglePureColor(int *fb, int x1, int y1, int x2, int y2, int x3, int y3, IUINT32 color);

//绘制纯色DownTriangle
void DrawDownTrianglePureColor(int *fb, int x1, int y1, int x2, int y2, int x3, int y3, IUINT32 color);



//绘制纯色三角形
void DrawTrianglePureColor2(int *fb, float x1, float y1, float x2, float y2, float x3, float y3, int color);


//绘制Grouraud三角形
void Draw_Gouraud_Triangle16(int *fb, POLYF4DV2_PTR face);