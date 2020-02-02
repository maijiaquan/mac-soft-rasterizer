#pragma once

#include "predefine.h"

#include <iostream>
using namespace std;

extern int display_w, display_h;
extern int *framebuffer; 

void device_pixel(int *fb, int x, int y, int c);
void device_draw_line(int *fb, int x1, int y1, int x2, int y2, int c);
void DrawTopTrianglePureColor(int *fb, int x1, int y1, int x2, int y2, int x3, int y3, IUINT32 color);
void DrawDownTrianglePureColor(int *fb, int x1, int y1, int x2, int y2, int x3, int y3, IUINT32 color);
void DrawTrianglePureColor2(int *fb, float x1, float y1, float x2, float y2, float x3, float y3, int color);
