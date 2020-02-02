#include "draw.h"

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
                    device_pixel(fb, x, y, c);
                }
            }
            device_pixel(fb, x2, y2, c);
        }
        else
        {
            if (y2 < y1)
                x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
            for (x = x1, y = y1; y <= y2; y++)
            {
                device_pixel(fb, x, y, c);
                rem += dx;
                if (rem >= dy)
                {
                    rem -= dy;
                    x += (x2 >= x1) ? 1 : -1;
                    device_pixel(fb, x, y, c);
                }
            }
            device_pixel(fb, x2, y2, c);
        }
    }
}
void DrawDownTrianglePureColor(int *fb, int x1, int y1, int x2, int y2, int x3, int y3, IUINT32 color)
{
        float dx_right, // the dx/dy ratio of the right edge of line
        dx_left,    // the dx/dy ratio of the left edge of line
        xs, xe,     // the starting and ending points of the edges
        height;     // the height of the triangle

    int temp_x, // used during sorting as temps
        temp_y,
        right, // used by clipping
        left;

    // test order of x1 and x2
    if (x3 < x2)
    {
        temp_x = x2;
        x2 = x3;
        x3 = temp_x;
    } // end if swap

    // compute delta's
    height = y3 - y1;

    dx_left = (x2 - x1) / height;
    dx_right = (x3 - x1) / height;

    // set starting points
    xs = (float)x1;
    xe = (float)x1; // +(float)0.5;

    // perform y clipping
    if (y1 < min_clip_y)
    {
        // compute new xs and ys
        xs = xs + dx_left * (float)(-y1 + min_clip_y);
        xe = xe + dx_right * (float)(-y1 + min_clip_y);

        // reset y1
        y1 = min_clip_y;

    } // end if top is off screen

    if (y3 > max_clip_y)
        y3 = max_clip_y;

    // compute starting address in video memory

    // test if x clipping is needed
    if (x1 >= min_clip_x && x1 <= max_clip_x &&
        x2 >= min_clip_x && x2 <= max_clip_x &&
        x3 >= min_clip_x && x3 <= max_clip_x)
    {
        // draw the triangle
        for (temp_y = y1; temp_y <= y3; temp_y++)
        {
            // draw the line
            device_draw_line(fb, xs, temp_y, xe, temp_y, color);

            // adjust starting point and ending point
            xs += dx_left;
            xe += dx_right;

        } // end for

    } // end if no x clipping needed
    else
    {
        // clip x axis with slower version

        // draw the triangle
        for (temp_y = y1; temp_y <= y3; temp_y++)
        {
            // do x clip
            left = (int)xs;
            right = (int)xe;

            // adjust starting point and ending point
            xs += dx_left;
            xe += dx_right;

            // clip line
            if (left < min_clip_x)
            {
                left = min_clip_x;

                if (right < min_clip_x)
                    continue;
            }

            if (right > max_clip_x)
            {
                right = max_clip_x;

                if (left > max_clip_x)
                    continue;
            }
            // draw the line
            device_draw_line(fb, left, temp_y, right, temp_y, color);
        } // end for

    } // end else x clipping needed
}

void DrawTopTrianglePureColor(int *fb, int x1, int y1, int x2, int y2, int x3, int y3, IUINT32 color)
{
   
    float dx_right, // the dx/dy ratio of the right edge of line
        dx_left,    // the dx/dy ratio of the left edge of line
        xs, xe,     // the starting and ending points of the edges
        height;     // the height of the triangle

    int temp_x, // used during sorting as temps
        temp_y,
        right, // used by clipping
        left;

    // test order of x1 and x2
    //保证 x1 < x2
    if (x2 < x1)
    {
        temp_x = x2;
        x2 = x1;
        x1 = temp_x;
    } // end if swap

    // compute delta's
    height = y3 - y1;

    dx_left = (x3 - x1) / height;
    dx_right = (x3 - x2) / height;

    // set starting points
    xs = (float)x1;
    xe = (float)x2; // +(float)0.5;

    // perform y clipping
    if (y1 < min_clip_y)
    {
        // compute new xs and ys
        xs = xs + dx_left * (float)(-y1 + min_clip_y);
        xe = xe + dx_right * (float)(-y1 + min_clip_y);

        // reset y1
        y1 = min_clip_y;

    } // end if top is off screen

    if (y3 > max_clip_y)
        y3 = max_clip_y;

    // compute starting address in video memory

    // test if x clipping is needed
    if (x1 >= min_clip_x && x1 <= max_clip_x &&
        x2 >= min_clip_x && x2 <= max_clip_x &&
        x3 >= min_clip_x && x3 <= max_clip_x)
    {
        // draw the triangle
        //for (temp_y = y1; temp_y <= y3; temp_y++, dest_addr += mempitch)
        for (temp_y = y1; temp_y <= y3; temp_y++)
        {
            // draw the line

            device_draw_line(fb, xs, temp_y, xe, temp_y, color);

            xs += dx_left;
            xe += dx_right;

        } // end for

    } // end if no x clipping needed
    else
    {
        // clip x axis with slower version

        // draw the triangle
        for (temp_y = y1; temp_y <= y3; temp_y++)
        {
            // do x clip
            left = (int)xs;
            right = (int)xe;

            // adjust starting point and ending point
            xs += dx_left;
            xe += dx_right;

            // clip line
            if (left < min_clip_x)
            {
                left = min_clip_x;

                if (right < min_clip_x)
                    continue;
            }

            if (right > max_clip_x)
            {
                right = max_clip_x;

                if (left > max_clip_x)
                    continue;
            }

            // draw the line
            // IUINT32 c = (0 << 16) | (255 << 8) | 0;
            device_draw_line(fb, left, temp_y, right, temp_y, color);
        } // end for

    } // end else x clipping needed 
}

void DrawTrianglePureColor2(int *fb, float x1, float y1, float x2, float y2, float x3, float y3, int color)
{
    //还原RGB值
    int r_base, g_base, b_base;
    _RGB565FROM16BIT(color, &r_base, &g_base, &b_base);
    // scale to 8 bit
    r_base <<= 3;
    g_base <<= 2;
    b_base <<= 3;
    IUINT32 c = (r_base << 16) | (g_base << 8) | b_base;

    int temp_x, // used for sorting
        temp_y,
        new_x;

    if ((FCMP(x1, x2) && FCMP(x2, x3)) || (FCMP(y1, y2) && FCMP(y2, y3)))
        return;
    // test for h lines and v lines
    // if ((x1 == x2 && x2 == x3) || (y1 == y2 && y2 == y3))
    // 	return;
    // sort p1,p2,p3 in ascending y order
    if (y2 < y1)
    {
        SWAP(x1, x2, temp_x);
        SWAP(y1, y2, temp_y);
    } // end if

    // now we know that p1 and p2 are in order
    if (y3 < y1)
    {
        SWAP(x1, x3, temp_x);
        SWAP(y1, y3, temp_y);
    } // end if

    // finally test y3 against y2
    if (y3 < y2)
    {
        SWAP(x2, x3, temp_x);
        SWAP(y2, y3, temp_y);
    } // end if

    // do trivial rejection tests for clipping
    if (y3 < min_clip_y || y1 > max_clip_y ||
        (x1 < min_clip_x && x2 < min_clip_x && x3 < min_clip_x) ||
        (x1 > max_clip_x && x2 > max_clip_x && x3 > max_clip_x))
        return;

    // test if top of triangle is flat
    if (FCMP(y1, y2))
    {

        DrawTopTrianglePureColor(fb, x1, y1, x2, y2, x3, y3, c);
        // DrawTopTriangle2(fb, x1, y1, x2, y2, x3, y3, c);
    } // end if
    else
    {
        if (FCMP(y2, y3))
        {
            DrawDownTrianglePureColor(fb, x1, y1, x2, y2, x3, y3, c);
            // DrawDownTriangle2(device, x1, y1, x2, y2, x3, y3, c);
        } // end if bottom is flat
        else
        {
            // draw each sub-triangle
            // new_x = x1 + (int)(0.5 + (float)(y2 - y1) * (float)(x3 - x1) / (float)(y3 - y1));
            new_x = x1 + (y2 - y1) * (x3 - x1) / (y3 - y1);
            DrawDownTrianglePureColor(fb, x1, y1, new_x, y2, x2, y2, c);
            // DrawDownTriangle2(device, x1, y1, new_x, y2, x2, y2, c);
            DrawTopTrianglePureColor(fb, x2, y2, new_x, y2, x3, y3, c);
            // DrawTopTriangle2(device, x2, y2, new_x, y2, x3, y3, c);
        } // end else
    }
}
}

// //等价于函数 Draw_Triangle_2D2_16
// void DrawTrianglePureColor2(int *fb, float x1, float y1, float x2, float y2, float x3, float y3, int color)
// {
//     //还原RGB值
//     int r_base, g_base, b_base;
//     _RGB565FROM16BIT(color, &r_base, &g_base, &b_base);
//     // scale to 8 bit
//     r_base <<= 3;
//     g_base <<= 2;
//     b_base <<= 3;
//     IUINT32 c = (r_base << 16) | (g_base << 8) | b_base;

//     int temp_x, // used for sorting
//         temp_y,
//         new_x;

//     if ((FCMP(x1, x2) && FCMP(x2, x3)) || (FCMP(y1, y2) && FCMP(y2, y3)))
//         return;
//     // test for h lines and v lines
//     // if ((x1 == x2 && x2 == x3) || (y1 == y2 && y2 == y3))
//     // 	return;
//     // sort p1,p2,p3 in ascending y order
//     if (y2 < y1)
//     {
//         SWAP(x1, x2, temp_x);
//         SWAP(y1, y2, temp_y);
//     } // end if

//     // now we know that p1 and p2 are in order
//     if (y3 < y1)
//     {
//         SWAP(x1, x3, temp_x);
//         SWAP(y1, y3, temp_y);
//     } // end if

//     // finally test y3 against y2
//     if (y3 < y2)
//     {
//         SWAP(x2, x3, temp_x);
//         SWAP(y2, y3, temp_y);
//     } // end if

//     // do trivial rejection tests for clipping
//     if (y3 < min_clip_y || y1 > max_clip_y ||
//         (x1 < min_clip_x && x2 < min_clip_x && x3 < min_clip_x) ||
//         (x1 > max_clip_x && x2 > max_clip_x && x3 > max_clip_x))
//         return;

//     // test if top of triangle is flat
//     if (FCMP(y1, y2))
//     {

//         DrawTopTrianglePureColor(fb, x1, y1, x2, y2, x3, y3, c);
//         // DrawTopTriangle2(fb, x1, y1, x2, y2, x3, y3, c);
//     } // end if
//     else if (FCMP(y2, y3))
//     {
//         DrawDownTrianglePureColor(fb, x1, y1, x2, y2, x3, y3, c);
//         // DrawDownTriangle2(device, x1, y1, x2, y2, x3, y3, c);
//     } // end if bottom is flat
//     else
//     {
//         // draw each sub-triangle
//         // new_x = x1 + (int)(0.5 + (float)(y2 - y1) * (float)(x3 - x1) / (float)(y3 - y1));
//         new_x = x1 + (y2 - y1) * (x3 - x1) / (y3 - y1);
//         DrawDownTrianglePureColor(fb, x1, y1, new_x, y2, x2, y2, c);
//         // DrawDownTriangle2(device, x1, y1, new_x, y2, x2, y2, c);
//         DrawTopTrianglePureColor(fb, x2, y2, new_x, y2, x3, y3, c);
//         // DrawTopTriangle2(device, x2, y2, new_x, y2, x3, y3, c);
//     } // end else
// }

// void DrawDownTrianglePureColor(int *fb, int x1, int y1, int x2, int y2, int x3, int y3, IUINT32 color)
// {
//     float dx_right, // the dx/dy ratio of the right edge of line
//         dx_left,    // the dx/dy ratio of the left edge of line
//         xs, xe,     // the starting and ending points of the edges
//         height;     // the height of the triangle

//     int temp_x, // used during sorting as temps
//         temp_y,
//         right, // used by clipping
//         left;

//     // test order of x1 and x2
//     if (x3 < x2)
//     {
//         temp_x = x2;
//         x2 = x3;
//         x3 = temp_x;
//     } // end if swap

//     // compute delta's
//     height = y3 - y1;

//     dx_left = (x2 - x1) / height;
//     dx_right = (x3 - x1) / height;

//     // set starting points
//     xs = (float)x1;
//     xe = (float)x1; // +(float)0.5;

//     // perform y clipping
//     if (y1 < min_clip_y)
//     {
//         // compute new xs and ys
//         xs = xs + dx_left * (float)(-y1 + min_clip_y);
//         xe = xe + dx_right * (float)(-y1 + min_clip_y);

//         // reset y1
//         y1 = min_clip_y;

//     } // end if top is off screen

//     if (y3 > max_clip_y)
//         y3 = max_clip_y;

//     // compute starting address in video memory

//     // test if x clipping is needed
//     if (x1 >= min_clip_x && x1 <= max_clip_x &&
//         x2 >= min_clip_x && x2 <= max_clip_x &&
//         x3 >= min_clip_x && x3 <= max_clip_x)
//     {
//         // draw the triangle
//         for (temp_y = y1; temp_y <= y3; temp_y++)
//         {
//             // draw the line
//             device_draw_line(fb, xs, temp_y, xe, temp_y, color);

//             // adjust starting point and ending point
//             xs += dx_left;
//             xe += dx_right;

//         } // end for

//     } // end if no x clipping needed
//     else
//     {
//         // clip x axis with slower version

//         // draw the triangle
//         for (temp_y = y1; temp_y <= y3; temp_y++)
//         {
//             // do x clip
//             left = (int)xs;
//             right = (int)xe;

//             // adjust starting point and ending point
//             xs += dx_left;
//             xe += dx_right;

//             // clip line
//             if (left < min_clip_x)
//             {
//                 left = min_clip_x;

//                 if (right < min_clip_x)
//                     continue;
//             }

//             if (right > max_clip_x)
//             {
//                 right = max_clip_x;

//                 if (left > max_clip_x)
//                     continue;
//             }
//             // draw the line
//             device_draw_line(fb, left, temp_y, right, temp_y, color);
//         } // end for

//     } // end else x clipping needed
// }

// void DrawTopTrianglePureColor(int *fb, int x1, int y1, int x2, int y2, int x3, int y3, IUINT32 color)
// {

//     float dx_right, // the dx/dy ratio of the right edge of line
//         dx_left,    // the dx/dy ratio of the left edge of line
//         xs, xe,     // the starting and ending points of the edges
//         height;     // the height of the triangle

//     int temp_x, // used during sorting as temps
//         temp_y,
//         right, // used by clipping
//         left;

//     // test order of x1 and x2
//     //保证 x1 < x2
//     if (x2 < x1)
//     {
//         temp_x = x2;
//         x2 = x1;
//         x1 = temp_x;
//     } // end if swap

//     // compute delta's
//     height = y3 - y1;

//     dx_left = (x3 - x1) / height;
//     dx_right = (x3 - x2) / height;

//     // set starting points
//     xs = (float)x1;
//     xe = (float)x2; // +(float)0.5;

//     // perform y clipping
//     if (y1 < min_clip_y)
//     {
//         // compute new xs and ys
//         xs = xs + dx_left * (float)(-y1 + min_clip_y);
//         xe = xe + dx_right * (float)(-y1 + min_clip_y);

//         // reset y1
//         y1 = min_clip_y;

//     } // end if top is off screen

//     if (y3 > max_clip_y)
//         y3 = max_clip_y;

//     // compute starting address in video memory

//     // test if x clipping is needed
//     if (x1 >= min_clip_x && x1 <= max_clip_x &&
//         x2 >= min_clip_x && x2 <= max_clip_x &&
//         x3 >= min_clip_x && x3 <= max_clip_x)
//     {
//         // draw the triangle
//         //for (temp_y = y1; temp_y <= y3; temp_y++, dest_addr += mempitch)
//         for (temp_y = y1; temp_y <= y3; temp_y++)
//         {
//             // draw the line

//             device_draw_line(fb, xs, temp_y, xe, temp_y, color);

//             xs += dx_left;
//             xe += dx_right;

//         } // end for

//     } // end if no x clipping needed
//     else
//     {
//         // clip x axis with slower version

//         // draw the triangle
//         for (temp_y = y1; temp_y <= y3; temp_y++)
//         {
//             // do x clip
//             left = (int)xs;
//             right = (int)xe;

//             // adjust starting point and ending point
//             xs += dx_left;
//             xe += dx_right;

//             // clip line
//             if (left < min_clip_x)
//             {
//                 left = min_clip_x;

//                 if (right < min_clip_x)
//                     continue;
//             }

//             if (right > max_clip_x)
//             {
//                 right = max_clip_x;

//                 if (left > max_clip_x)
//                     continue;
//             }

//             // draw the line
//             // IUINT32 c = (0 << 16) | (255 << 8) | 0;
//             device_draw_line(fb, left, temp_y, right, temp_y, color);
//         } // end for

//     } // end else x clipping needed
// }
