#include"draw.h"


// void device_pixel(int x, int y, int c);
void device_pixel(int *fb, int x, int y, int c)
// void device_pixel(int x, int y, int c)
{
    // framebuffer[y * display_w + x] = c;
    fb[y * display_w + x] = c;
}


void device_draw_line(int *fb, int x1, int y1, int x2, int y2, int c)
{
    // std::cout<<"x1 = "<<x1<<"y1 = "<<y1<<std::endl;
    // std::cout<<"x2 = "<<x2<<"y2 = "<<y2<<std::endl;

    int x, y, rem = 0;
    if (x1 == x2 && y1 == y2)
    {
        device_pixel(fb, x1, y1, c);
    }
    else if (x1 == x2)
    {
        int inc = (y1 <= y2) ? 1 : -1;
        for (y = y1; y != y2; y += inc)
            device_pixel(fb, x1, y, c);
        device_pixel(fb, x2, y2, c);
    }
    else if (y1 == y2)
    {
        int inc = (x1 <= x2) ? 1 : -1;
        for (x = x1; x != x2; x += inc)
            device_pixel(fb, x, y1, c);
        device_pixel(fb, x2, y2, c);
    }
    else
    {
        int dx = (x1 < x2) ? x2 - x1 : x1 - x2;
        int dy = (y1 < y2) ? y2 - y1 : y1 - y2;
        if (dx >= dy)
        {
            if (x2 < x1)
                x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
            for (x = x1, y = y1; x <= x2; x++)
            {
                device_pixel(fb, x, y, c);
                rem += dy;
                if (rem >= dx)
                {
                    rem -= dx;
                    y += (y2 >= y1) ? 1 : -1;
                    device_pixel(fb,x, y, c);
                }
            }
            device_pixel(fb,x2, y2, c);
        }
        else
        {
            if (y2 < y1)
                x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
            for (x = x1, y = y1; y <= y2; y++)
            {
                device_pixel(fb,x, y, c);
                rem += dx;
                if (rem >= dy)
                {
                    rem -= dy;
                    x += (x2 >= x1) ? 1 : -1;
                    device_pixel(fb,x, y, c);
                }
            }
            device_pixel(fb,x2, y2, c);
        }
    }
}