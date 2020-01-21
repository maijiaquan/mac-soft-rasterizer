#pragma once

#include <iostream>
using namespace std;

extern int display_w, display_h;
extern int *framebuffer; //帧缓

void device_pixel(int *fb, int x, int y, int c);
void device_draw_line(int *fb, int x1, int y1, int x2, int y2, int c);