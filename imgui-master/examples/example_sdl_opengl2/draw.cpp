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

void DrawTriangleWireframe(int *fb, float x1, float y1, float x2, float y2, float x3, float y3, int color)
{
    device_draw_line(fb, x1, y1, x2, y2, color); //3 1
    device_draw_line(fb, x1, y1, x3, y3, color); //3 1
    device_draw_line(fb, x2, y2, x3, y3, color); //3 1
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

void Draw_Gouraud_Triangle16(int *fb, POLYF4DV2_PTR face)
{

    int v0 = 0,
        v1 = 1,
        v2 = 2,
        temp = 0,
        tri_type = TRI_TYPE_NONE,
        irestart = INTERP_LHS;

    int dx, dy, dyl, dyr, // general deltas
        // u, v, w,
        du, dv, dw,
        xi, yi,           // the current interpolated x,y
        ui, vi, wi,       // the current interpolated u,v
        index_x, index_y, // looping vars
        x, y,             // hold general x,y
        xstart,
        xend,
        ystart,
        yrestart,
        yend,
        xl,
        dxdyl,
        xr,
        dxdyr,
        dudyl,
        ul,
        dvdyl,
        vl,
        dwdyl,
        wl,
        dudyr,
        ur,
        dvdyr,
        vr,
        dwdyr,
        wr;

    int x0, y0, tu0, tv0, tw0, // cached vertices
        x1, y1, tu1, tv1, tw1,
        x2, y2, tu2, tv2, tw2;

    int r_base0, g_base0, b_base0,
        r_base1, g_base1, b_base1,
        r_base2, g_base2, b_base2;

    //USHORT *screen_ptr = NULL,
    //       *screen_line = NULL,
    //       *textmap = NULL,
    //       *dest_buffer = (USHORT *)_dest_buffer;

    // adjust memory pitch to words, divide by 2

    // first trivial clipping rejection tests
    if (((face->tvlist[0].y < min_clip_y) &&
         (face->tvlist[1].y < min_clip_y) &&
         (face->tvlist[2].y < min_clip_y)) ||

        ((face->tvlist[0].y > max_clip_y) &&
         (face->tvlist[1].y > max_clip_y) &&
         (face->tvlist[2].y > max_clip_y)) ||

        ((face->tvlist[0].x < min_clip_x) &&
         (face->tvlist[1].x < min_clip_x) &&
         (face->tvlist[2].x < min_clip_x)) ||

        ((face->tvlist[0].x > max_clip_x) &&
         (face->tvlist[1].x > max_clip_x) &&
         (face->tvlist[2].x > max_clip_x)))
        return;

    // degenerate triangle
    if (((face->tvlist[0].x == face->tvlist[1].x) && (face->tvlist[1].x == face->tvlist[2].x)) ||
        ((face->tvlist[0].y == face->tvlist[1].y) && (face->tvlist[1].y == face->tvlist[2].y)))
        return;

    // sort vertices
    if (face->tvlist[v1].y < face->tvlist[v0].y)
    {
        SWAP(v0, v1, temp);
    }

    if (face->tvlist[v2].y < face->tvlist[v0].y)
    {
        SWAP(v0, v2, temp);
    }

    if (face->tvlist[v2].y < face->tvlist[v1].y)
    {
        SWAP(v1, v2, temp);
    }

    // // now test for trivial flat sided cases
    if (face->tvlist[v0].y == face->tvlist[v1].y)
    {
        // set triangle type
        tri_type = TRI_TYPE_FLAT_TOP;

        // sort vertices left to right
        if (face->tvlist[v1].x < face->tvlist[v0].x)
        {
            SWAP(v0, v1, temp);
        }

    } // end if
    else if (face->tvlist[v1].y == face->tvlist[v2].y)
    {
        // set triangle type
        tri_type = TRI_TYPE_FLAT_BOTTOM;

        // sort vertices left to right
        if (face->tvlist[v2].x < face->tvlist[v1].x)
        {
            SWAP(v1, v2, temp);
        }

    } // end if
    else
    {
        // must be a general triangle
        tri_type = TRI_TYPE_GENERAL;

    } // end else

    // // assume 5.6.5 format -- sorry!
    // // we can't afford a function call in the inner loops, so we must write
    // // two hard coded versions, if we want support for both 5.6.5, and 5.5.5
    // cout<<endl;
    _RGB565FROM16BIT(face->lit_color[v0], &r_base0, &g_base0, &b_base0);
    _RGB565FROM16BIT(face->lit_color[v1], &r_base1, &g_base1, &b_base1);
    _RGB565FROM16BIT(face->lit_color[v2], &r_base2, &g_base2, &b_base2);

    // // scale to 8 bit
    r_base0 <<= 3;
    g_base0 <<= 2;
    b_base0 <<= 3;

    // scale to 8 bit
    r_base1 <<= 3;
    g_base1 <<= 2;
    b_base1 <<= 3;

    // scale to 8 bit
    r_base2 <<= 3;
    g_base2 <<= 2;
    b_base2 <<= 3;

    // extract vertices for processing, now that we have order
    x0 = (int)(face->tvlist[v0].x + 0.5);
    y0 = (int)(face->tvlist[v0].y + 0.5);

    tu0 = r_base0;
    tv0 = g_base0;
    tw0 = b_base0;

    x1 = (int)(face->tvlist[v1].x + 0.5);
    y1 = (int)(face->tvlist[v1].y + 0.5);

    tu1 = r_base1;
    tv1 = g_base1;
    tw1 = b_base1;

    x2 = (int)(face->tvlist[v2].x + 0.5);
    y2 = (int)(face->tvlist[v2].y + 0.5);

    tu2 = r_base2;
    tv2 = g_base2;
    tw2 = b_base2;

    // set interpolation restart value
    yrestart = y1;

    // what kind of triangle
    if (tri_type & TRI_TYPE_FLAT_MASK)
    {
        // std::cout<<"special "<<std::endl;

        if (tri_type == TRI_TYPE_FLAT_TOP)
        {
            // compute all deltas
            dy = (y2 - y0);

            dxdyl = ((x2 - x0) << FIXP16_SHIFT) / dy;
            dudyl = ((tu2 - tu0) << FIXP16_SHIFT) / dy;
            dvdyl = ((tv2 - tv0) << FIXP16_SHIFT) / dy;
            dwdyl = ((tw2 - tw0) << FIXP16_SHIFT) / dy;

            dxdyr = ((x2 - x1) << FIXP16_SHIFT) / dy;
            dudyr = ((tu2 - tu1) << FIXP16_SHIFT) / dy;
            dvdyr = ((tv2 - tv1) << FIXP16_SHIFT) / dy;
            dwdyr = ((tw2 - tw1) << FIXP16_SHIFT) / dy;

            // test for y clipping
            if (y0 < min_clip_y)
            {
                // compute overclip
                dy = (min_clip_y - y0);

                // computer new LHS starting values
                xl = dxdyl * dy + (x0 << FIXP16_SHIFT);
                ul = dudyl * dy + (tu0 << FIXP16_SHIFT);
                vl = dvdyl * dy + (tv0 << FIXP16_SHIFT);
                wl = dwdyl * dy + (tw0 << FIXP16_SHIFT);

                // compute new RHS starting values
                xr = dxdyr * dy + (x1 << FIXP16_SHIFT);
                ur = dudyr * dy + (tu1 << FIXP16_SHIFT);
                vr = dvdyr * dy + (tv1 << FIXP16_SHIFT);
                wr = dwdyr * dy + (tw1 << FIXP16_SHIFT);

                // compute new starting y
                ystart = min_clip_y;

            } // end if
            else
            {
                // no clipping

                // set starting values
                xl = (x0 << FIXP16_SHIFT);
                xr = (x1 << FIXP16_SHIFT);

                ul = (tu0 << FIXP16_SHIFT);
                vl = (tv0 << FIXP16_SHIFT);
                wl = (tw0 << FIXP16_SHIFT);

                ur = (tu1 << FIXP16_SHIFT);
                vr = (tv1 << FIXP16_SHIFT);
                wr = (tw1 << FIXP16_SHIFT);

                // set starting y
                ystart = y0;

            } // end else

        } // end if flat top
        else
        {
            // must be flat bottom

            // compute all deltas
            dy = (y1 - y0);

            dxdyl = ((x1 - x0) << FIXP16_SHIFT) / dy;
            dudyl = ((tu1 - tu0) << FIXP16_SHIFT) / dy;
            dvdyl = ((tv1 - tv0) << FIXP16_SHIFT) / dy;
            dwdyl = ((tw1 - tw0) << FIXP16_SHIFT) / dy;

            dxdyr = ((x2 - x0) << FIXP16_SHIFT) / dy;
            dudyr = ((tu2 - tu0) << FIXP16_SHIFT) / dy;
            dvdyr = ((tv2 - tv0) << FIXP16_SHIFT) / dy;
            dwdyr = ((tw2 - tw0) << FIXP16_SHIFT) / dy;

            // test for y clipping
            if (y0 < min_clip_y)
            {
                // compute overclip
                dy = (min_clip_y - y0);

                // computer new LHS starting values
                xl = dxdyl * dy + (x0 << FIXP16_SHIFT);
                ul = dudyl * dy + (tu0 << FIXP16_SHIFT);
                vl = dvdyl * dy + (tv0 << FIXP16_SHIFT);
                wl = dwdyl * dy + (tw0 << FIXP16_SHIFT);

                // compute new RHS starting values
                xr = dxdyr * dy + (x0 << FIXP16_SHIFT);
                ur = dudyr * dy + (tu0 << FIXP16_SHIFT);
                vr = dvdyr * dy + (tv0 << FIXP16_SHIFT);
                wr = dwdyr * dy + (tw0 << FIXP16_SHIFT);

                // compute new starting y
                ystart = min_clip_y;

            } // end if
            else
            {
                // no clipping

                // set starting values
                xl = (x0 << FIXP16_SHIFT);
                xr = (x0 << FIXP16_SHIFT);

                ul = (tu0 << FIXP16_SHIFT);
                vl = (tv0 << FIXP16_SHIFT);
                wl = (tw0 << FIXP16_SHIFT);

                ur = (tu0 << FIXP16_SHIFT);
                vr = (tv0 << FIXP16_SHIFT);
                wr = (tw0 << FIXP16_SHIFT);

                // set starting y
                ystart = y0;

            } // end else

        } // end else flat bottom

        // // test for bottom clip, always
        if ((yend = y2) > max_clip_y)
            yend = max_clip_y;

        // test for horizontal clipping
        if ((x0 < min_clip_x) || (x0 > max_clip_x) ||
            (x1 < min_clip_x) || (x1 > max_clip_x) ||
            (x2 < min_clip_x) || (x2 > max_clip_x))
        {
            // clip version

            // point screen ptr to starting line
            // screen_ptr = dest_buffer + (ystart * mem_pitch);

            for (yi = ystart; yi <= yend; yi++)
            {
                // compute span endpoints
                xstart = ((xl + FIXP16_ROUND_UP) >> FIXP16_SHIFT);
                xend = ((xr + FIXP16_ROUND_UP) >> FIXP16_SHIFT);

                // compute starting points for u,v,w interpolants
                ui = ul + FIXP16_ROUND_UP;
                vi = vl + FIXP16_ROUND_UP;
                wi = wl + FIXP16_ROUND_UP;

                // compute u,v interpolants
                if ((dx = (xend - xstart)) > 0)
                {
                    du = (ur - ul) / dx;
                    dv = (vr - vl) / dx;
                    dw = (wr - wl) / dx;
                } // end if
                else
                {
                    du = (ur - ul);
                    dv = (vr - vl);
                    dw = (wr - wl);
                } // end else

                ///////////////////////////////////////////////////////////////////////

                // test for x clipping, LHS
                if (xstart < min_clip_x)
                {
                    // compute x overlap
                    dx = min_clip_x - xstart;

                    // slide interpolants over
                    ui += dx * du;
                    vi += dx * dv;
                    wi += dx * dw;

                    // reset vars
                    xstart = min_clip_x;

                } // end if

                // test for x clipping RHS
                if (xend > max_clip_x)
                    xend = max_clip_x;

                ///////////////////////////////////////////////////////////////////////

                // draw span
                for (xi = xstart; xi <= xend; xi++)
                {
                    // write textel assume 5.6.5
                    // screen_ptr[xi] = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));

                    int color = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));
                    IUINT32 c;
                    RGBFrom565(color, c);
                    device_pixel(fb, xi, yi, c);
                    // interpolate u,v
                    // IUINT32 c = (ui<<16) | (vi<<8) | wi;
                    // IUINT32 c = (red << 16) | (green << 8) | blue;
                    // IUINT32 c = ((ui >> (FIXP16_SHIFT + 3))<<16) | ((vi >> (FIXP16_SHIFT + 2))<<8) | (wi >> (FIXP16_SHIFT + 3));
                    ui += du;
                    vi += dv;
                    wi += dw;
                } // end for xi

                // interpolate u,v,w,x along right and left edge
                xl += dxdyl;
                ul += dudyl;
                vl += dvdyl;
                wl += dwdyl;

                xr += dxdyr;
                ur += dudyr;
                vr += dvdyr;
                wr += dwdyr;

                // advance screen ptr

            } // end for y

        } // end if clip
        else
        {
            // non-clip version

            // point screen ptr to starting line
            // screen_ptr = dest_buffer + (ystart * mem_pitch);

            for (yi = ystart; yi <= yend; yi++)
            {
                // compute span endpoints
                xstart = ((xl + FIXP16_ROUND_UP) >> FIXP16_SHIFT);
                xend = ((xr + FIXP16_ROUND_UP) >> FIXP16_SHIFT);

                // compute starting points for u,v,w interpolants
                ui = ul + FIXP16_ROUND_UP;
                vi = vl + FIXP16_ROUND_UP;
                wi = wl + FIXP16_ROUND_UP;

                // compute u,v interpolants
                if ((dx = (xend - xstart)) > 0)
                {
                    du = (ur - ul) / dx;
                    dv = (vr - vl) / dx;
                    dw = (wr - wl) / dx;
                } // end if
                else
                {
                    du = (ur - ul);
                    dv = (vr - vl);
                    dw = (wr - wl);
                } // end else

                // draw span
                for (xi = xstart; xi <= xend; xi++)
                {
                    // write textel 5.6.5
                    // screen_ptr[xi] = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));

                    // IUINT32 c = (ui<<16) | (vi<<8) | wi;
                    // IUINT32 c = ((ui >> (FIXP16_SHIFT + 3))<<16) | ((vi >> (FIXP16_SHIFT + 2))<<8) | (wi >> (FIXP16_SHIFT + 3));
                    // device_pixel(fb, xi,  yi,  c);

                    int color = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));
                    IUINT32 c;
                    RGBFrom565(color, c);
                    device_pixel(fb, xi, yi, c);
                    // interpolate u,v
                    ui += du;
                    vi += dv;
                    wi += dw;
                } // end for xi

                // interpolate u,v,w,x along right and left edge
                xl += dxdyl;
                ul += dudyl;
                vl += dvdyl;
                wl += dwdyl;

                xr += dxdyr;
                ur += dudyr;
                vr += dvdyr;
                wr += dwdyr;

                // advance screen ptr
                // screen_ptr += mem_pitch;

            } // end for y

        } // end if non-clipped

    } // end if
    else if (tri_type == TRI_TYPE_GENERAL)
    {

        // first test for bottom clip, always
        if ((yend = y2) > max_clip_y)
            yend = max_clip_y;

        // pre-test y clipping status
        if (y1 < min_clip_y)
        {
            std::cout << "y1 < min_clip_y" << std::endl;
            // compute all deltas
            // LHS
            dyl = (y2 - y1);

            dxdyl = ((x2 - x1) << FIXP16_SHIFT) / dyl;
            dudyl = ((tu2 - tu1) << FIXP16_SHIFT) / dyl;
            dvdyl = ((tv2 - tv1) << FIXP16_SHIFT) / dyl;
            dwdyl = ((tw2 - tw1) << FIXP16_SHIFT) / dyl;

            // RHS
            dyr = (y2 - y0);

            dxdyr = ((x2 - x0) << FIXP16_SHIFT) / dyr;
            dudyr = ((tu2 - tu0) << FIXP16_SHIFT) / dyr;
            dvdyr = ((tv2 - tv0) << FIXP16_SHIFT) / dyr;
            dwdyr = ((tw2 - tw0) << FIXP16_SHIFT) / dyr;

            // compute overclip
            dyr = (min_clip_y - y0);
            dyl = (min_clip_y - y1);

            // computer new LHS starting values
            xl = dxdyl * dyl + (x1 << FIXP16_SHIFT);

            ul = dudyl * dyl + (tu1 << FIXP16_SHIFT);
            vl = dvdyl * dyl + (tv1 << FIXP16_SHIFT);
            wl = dwdyl * dyl + (tw1 << FIXP16_SHIFT);

            // compute new RHS starting values
            xr = dxdyr * dyr + (x0 << FIXP16_SHIFT);

            ur = dudyr * dyr + (tu0 << FIXP16_SHIFT);
            vr = dvdyr * dyr + (tv0 << FIXP16_SHIFT);
            wr = dwdyr * dyr + (tw0 << FIXP16_SHIFT);

            // compute new starting y
            ystart = min_clip_y;

            // test if we need swap to keep rendering left to right
            if (dxdyr > dxdyl)
            {
                SWAP(dxdyl, dxdyr, temp);
                SWAP(dudyl, dudyr, temp);
                SWAP(dvdyl, dvdyr, temp);
                SWAP(dwdyl, dwdyr, temp);
                SWAP(xl, xr, temp);
                SWAP(ul, ur, temp);
                SWAP(vl, vr, temp);
                SWAP(wl, wr, temp);
                SWAP(x1, x2, temp);
                SWAP(y1, y2, temp);
                SWAP(tu1, tu2, temp);
                SWAP(tv1, tv2, temp);
                SWAP(tw1, tw2, temp);

                // set interpolation restart
                irestart = INTERP_RHS;

            } // end if

        } // end if
        else if (y0 < min_clip_y)
        {
            std::cout << "y0 < min_clip_y" << std::endl;
            // compute all deltas
            // LHS
            dyl = (y1 - y0);

            dxdyl = ((x1 - x0) << FIXP16_SHIFT) / dyl;
            dudyl = ((tu1 - tu0) << FIXP16_SHIFT) / dyl;
            dvdyl = ((tv1 - tv0) << FIXP16_SHIFT) / dyl;
            dwdyl = ((tw1 - tw0) << FIXP16_SHIFT) / dyl;

            // RHS
            dyr = (y2 - y0);

            dxdyr = ((x2 - x0) << FIXP16_SHIFT) / dyr;
            dudyr = ((tu2 - tu0) << FIXP16_SHIFT) / dyr;
            dvdyr = ((tv2 - tv0) << FIXP16_SHIFT) / dyr;
            dwdyr = ((tw2 - tw0) << FIXP16_SHIFT) / dyr;

            // compute overclip
            dy = (min_clip_y - y0);

            // computer new LHS starting values
            xl = dxdyl * dy + (x0 << FIXP16_SHIFT);
            ul = dudyl * dy + (tu0 << FIXP16_SHIFT);
            vl = dvdyl * dy + (tv0 << FIXP16_SHIFT);
            wl = dwdyl * dy + (tw0 << FIXP16_SHIFT);

            // compute new RHS starting values
            xr = dxdyr * dy + (x0 << FIXP16_SHIFT);
            ur = dudyr * dy + (tu0 << FIXP16_SHIFT);
            vr = dvdyr * dy + (tv0 << FIXP16_SHIFT);
            wr = dwdyr * dy + (tw0 << FIXP16_SHIFT);

            // compute new starting y
            ystart = min_clip_y;

            // test if we need swap to keep rendering left to right
            if (dxdyr < dxdyl)
            {
                SWAP(dxdyl, dxdyr, temp);
                SWAP(dudyl, dudyr, temp);
                SWAP(dvdyl, dvdyr, temp);
                SWAP(dwdyl, dwdyr, temp);
                SWAP(xl, xr, temp);
                SWAP(ul, ur, temp);
                SWAP(vl, vr, temp);
                SWAP(wl, wr, temp);
                SWAP(x1, x2, temp);
                SWAP(y1, y2, temp);
                SWAP(tu1, tu2, temp);
                SWAP(tv1, tv2, temp);
                SWAP(tw1, tw2, temp);

                // set interpolation restart
                irestart = INTERP_RHS;

            } // end if

        } // end if
        else
        {

            // no initial y clipping

            // compute all deltas
            // LHS
            dyl = (y1 - y0);

            dxdyl = ((x1 - x0) << FIXP16_SHIFT) / dyl;
            dudyl = ((tu1 - tu0) << FIXP16_SHIFT) / dyl;
            dvdyl = ((tv1 - tv0) << FIXP16_SHIFT) / dyl;
            dwdyl = ((tw1 - tw0) << FIXP16_SHIFT) / dyl;

            // RHS
            dyr = (y2 - y0);

            dxdyr = ((x2 - x0) << FIXP16_SHIFT) / dyr;
            dudyr = ((tu2 - tu0) << FIXP16_SHIFT) / dyr;
            dvdyr = ((tv2 - tv0) << FIXP16_SHIFT) / dyr;
            dwdyr = ((tw2 - tw0) << FIXP16_SHIFT) / dyr;

            // no clipping y

            // set starting values
            xl = (x0 << FIXP16_SHIFT);
            xr = (x0 << FIXP16_SHIFT);

            ul = (tu0 << FIXP16_SHIFT);
            vl = (tv0 << FIXP16_SHIFT);
            wl = (tw0 << FIXP16_SHIFT);

            ur = (tu0 << FIXP16_SHIFT);
            vr = (tv0 << FIXP16_SHIFT);
            wr = (tw0 << FIXP16_SHIFT);

            // set starting y
            ystart = y0;

            // test if we need swap to keep rendering left to right
            if (dxdyr < dxdyl)
            {
                SWAP(dxdyl, dxdyr, temp);
                SWAP(dudyl, dudyr, temp);
                SWAP(dvdyl, dvdyr, temp);
                SWAP(dwdyl, dwdyr, temp);
                SWAP(xl, xr, temp);
                SWAP(ul, ur, temp);
                SWAP(vl, vr, temp);
                SWAP(wl, wr, temp);
                SWAP(x1, x2, temp);
                SWAP(y1, y2, temp);
                SWAP(tu1, tu2, temp);
                SWAP(tv1, tv2, temp);
                SWAP(tw1, tw2, temp);

                // set interpolation restart
                irestart = INTERP_RHS;

            } // end if

        } // end else

        // test for horizontal clipping
        if ((x0 < min_clip_x) || (x0 > max_clip_x) ||
            (x1 < min_clip_x) || (x1 > max_clip_x) ||
            (x2 < min_clip_x) || (x2 > max_clip_x))
        {
            // clip version
            // x clipping

            // point screen ptr to starting line
            // screen_ptr = dest_buffer + (ystart * mem_pitch);

            for (yi = ystart; yi <= yend; yi++)
            {
                // compute span endpoints
                xstart = ((xl + FIXP16_ROUND_UP) >> FIXP16_SHIFT);
                xend = ((xr + FIXP16_ROUND_UP) >> FIXP16_SHIFT);

                // compute starting points for u,v,w interpolants
                ui = ul + FIXP16_ROUND_UP;
                vi = vl + FIXP16_ROUND_UP;
                wi = wl + FIXP16_ROUND_UP;

                // compute u,v interpolants
                if ((dx = (xend - xstart)) > 0)
                {
                    du = (ur - ul) / dx;
                    dv = (vr - vl) / dx;
                    dw = (wr - wl) / dx;
                } // end if
                else
                {
                    du = (ur - ul);
                    dv = (vr - vl);
                    dw = (wr - wl);
                } // end else

                ///////////////////////////////////////////////////////////////////////

                // test for x clipping, LHS
                if (xstart < min_clip_x)
                {
                    // compute x overlap
                    dx = min_clip_x - xstart;

                    // slide interpolants over
                    ui += dx * du;
                    vi += dx * dv;
                    wi += dx * dw;

                    // set x to left clip edge
                    xstart = min_clip_x;

                } // end if

                // test for x clipping RHS
                if (xend > max_clip_x)
                    xend = max_clip_x;

                ///////////////////////////////////////////////////////////////////////

                // draw span
                for (xi = xstart; xi <= xend; xi++)
                {
                    // write textel assume 5.6.5
                    // screen_ptr[xi] = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));

                    // IUINT32 c = (ui<<16) | (vi<<8) | wi;
                    // IUINT32 c = ((ui >> (FIXP16_SHIFT + 3))<<16) | ((vi >> (FIXP16_SHIFT + 2))<<8) | (wi >> (FIXP16_SHIFT + 3));
                    // device_pixel(fb, xi,  yi,  c);

                    int color = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));
                    IUINT32 c;
                    RGBFrom565(color, c);
                    device_pixel(fb, xi, yi, c);
                    // interpolate u,v
                    ui += du;
                    vi += dv;
                    wi += dw;
                } // end for xi

                // interpolate u,v,w,x along right and left edge
                xl += dxdyl;
                ul += dudyl;
                vl += dvdyl;
                wl += dwdyl;

                xr += dxdyr;
                ur += dudyr;
                vr += dvdyr;
                wr += dwdyr;

                // advance screen ptr
                // screen_ptr += mem_pitch;

                // test for yi hitting second region, if so change interpolant
                if (yi == yrestart)
                {
                    // test interpolation side change flag

                    if (irestart == INTERP_LHS)
                    {
                        // LHS
                        dyl = (y2 - y1);

                        dxdyl = ((x2 - x1) << FIXP16_SHIFT) / dyl;
                        dudyl = ((tu2 - tu1) << FIXP16_SHIFT) / dyl;
                        dvdyl = ((tv2 - tv1) << FIXP16_SHIFT) / dyl;
                        dwdyl = ((tw2 - tw1) << FIXP16_SHIFT) / dyl;

                        // set starting values
                        xl = (x1 << FIXP16_SHIFT);
                        ul = (tu1 << FIXP16_SHIFT);
                        vl = (tv1 << FIXP16_SHIFT);
                        wl = (tw1 << FIXP16_SHIFT);

                        // interpolate down on LHS to even up
                        xl += dxdyl;
                        ul += dudyl;
                        vl += dvdyl;
                        wl += dwdyl;
                    } // end if
                    else
                    {
                        // RHS
                        dyr = (y1 - y2);

                        dxdyr = ((x1 - x2) << FIXP16_SHIFT) / dyr;
                        dudyr = ((tu1 - tu2) << FIXP16_SHIFT) / dyr;
                        dvdyr = ((tv1 - tv2) << FIXP16_SHIFT) / dyr;
                        dwdyr = ((tw1 - tw2) << FIXP16_SHIFT) / dyr;

                        // set starting values
                        xr = (x2 << FIXP16_SHIFT);
                        ur = (tu2 << FIXP16_SHIFT);
                        vr = (tv2 << FIXP16_SHIFT);
                        wr = (tw2 << FIXP16_SHIFT);

                        // interpolate down on RHS to even up
                        xr += dxdyr;
                        ur += dudyr;
                        vr += dvdyr;
                        wr += dwdyr;

                    } // end else

                } // end if

            } // end for y

        } // end if
        else
        {
            // no x clipping
            // point screen ptr to starting line
            // screen_ptr = dest_buffer + (ystart * mem_pitch);

            for (yi = ystart; yi <= yend; yi++)
            {
                // compute span endpoints
                xstart = ((xl + FIXP16_ROUND_UP) >> FIXP16_SHIFT);
                xend = ((xr + FIXP16_ROUND_UP) >> FIXP16_SHIFT);

                // compute starting points for u,v,w interpolants
                ui = ul + FIXP16_ROUND_UP;
                vi = vl + FIXP16_ROUND_UP;
                wi = wl + FIXP16_ROUND_UP;

                // compute u,v interpolants
                if ((dx = (xend - xstart)) > 0)
                {
                    du = (ur - ul) / dx;
                    dv = (vr - vl) / dx;
                    dw = (wr - wl) / dx;
                } // end if
                else
                {
                    du = (ur - ul);
                    dv = (vr - vl);
                    dw = (wr - wl);
                } // end else

                // draw span
                for (xi = xstart; xi <= xend; xi++)
                {
                    // write textel assume 5.6.5
                    // screen_ptr[xi] = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));

                    // IUINT32 c = (ui<<16) | (vi<<8) | wi;
                    // IUINT32 c = ((ui >> (FIXP16_SHIFT + 3))<<16) | ((vi >> (FIXP16_SHIFT + 2))<<8) | (wi >> (FIXP16_SHIFT + 3));
                    // device_pixel(fb, xi,  yi,  c);

                    int color = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));
                    IUINT32 c;
                    RGBFrom565(color, c);
                    device_pixel(fb, xi, yi, c);
                    // interpolate u,v
                    ui += du;
                    vi += dv;
                    wi += dw;
                } // end for xi

                // interpolate u,v,w,x along right and left edge
                xl += dxdyl;
                ul += dudyl;
                vl += dvdyl;
                wl += dwdyl;

                xr += dxdyr;
                ur += dudyr;
                vr += dvdyr;
                wr += dwdyr;

                // advance screen ptr
                // screen_ptr += mem_pitch;

                // test for yi hitting second region, if so change interpolant
                if (yi == yrestart)
                {
                    // test interpolation side change flag

                    if (irestart == INTERP_LHS)
                    {
                        // LHS
                        dyl = (y2 - y1);

                        dxdyl = ((x2 - x1) << FIXP16_SHIFT) / dyl;
                        dudyl = ((tu2 - tu1) << FIXP16_SHIFT) / dyl;
                        dvdyl = ((tv2 - tv1) << FIXP16_SHIFT) / dyl;
                        dwdyl = ((tw2 - tw1) << FIXP16_SHIFT) / dyl;

                        // set starting values
                        xl = (x1 << FIXP16_SHIFT);
                        ul = (tu1 << FIXP16_SHIFT);
                        vl = (tv1 << FIXP16_SHIFT);
                        wl = (tw1 << FIXP16_SHIFT);

                        // interpolate down on LHS to even up
                        xl += dxdyl;
                        ul += dudyl;
                        vl += dvdyl;
                        wl += dwdyl;
                    } // end if
                    else
                    {
                        // RHS
                        dyr = (y1 - y2);

                        dxdyr = ((x1 - x2) << FIXP16_SHIFT) / dyr;
                        dudyr = ((tu1 - tu2) << FIXP16_SHIFT) / dyr;
                        dvdyr = ((tv1 - tv2) << FIXP16_SHIFT) / dyr;
                        dwdyr = ((tw1 - tw2) << FIXP16_SHIFT) / dyr;

                        // set starting values
                        xr = (x2 << FIXP16_SHIFT);
                        ur = (tu2 << FIXP16_SHIFT);
                        vr = (tv2 << FIXP16_SHIFT);
                        wr = (tw2 << FIXP16_SHIFT);

                        // interpolate down on RHS to even up
                        xr += dxdyr;
                        ur += dudyr;
                        vr += dvdyr;
                        wr += dwdyr;
                    } // end else

                } // end if

            } // end for y

        } // end else

    } // end if
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
// }