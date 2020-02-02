// dear imgui: standalone example application for SDL2 + OpenGL
// If you are new to dear imgui, see examples/README.txt and documentation at the top of imgui.cpp.
// (SDL is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan graphics context creation, etc.)

// **DO NOT USE THIS CODE IF YOUR CODE/ENGINE IS USING MODERN OPENGL (SHADERS, VBO, VAO, etc.)**
// **Prefer using the code in the example_sdl_opengl3/ folder**
// See imgui_impl_sdl.cpp for details.

#include "imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl2.h"
#include <stdio.h>
#include <SDL.h>
#include <SDL_opengl.h>
//#include "string.h"
#include "ds.h"
#include "draw.h"
#include "cobloader.h"
//using namespace std;

//全局变量
int *framebuffer; //帧缓冲
ImVec4 clear_color = ImVec4(0.0f, 0.0f, 0.00f, 1.00f);
int display_w, display_h;

POINT4D cam_pos = {0, 0, -100, 1};
VECTOR4D cam_dir = {0, 0, 0, 1};

// all your initialization code goes here...
VECTOR4D vscale = {.5, .5, .5, 1},
         vpos = {0, 0, 0, 1},
         vrot = {0, 0, 0, 1};

RENDERLIST4DV1 rend_list;           // the single renderlist
RENDERLIST4DV2 rend_list2;          // the render list
POLYF4DV1 poly1;                    // our lonely polygon
CAM4DV1 cam;                        // the single camera
POINT4D poly1_pos = {0, 0, 100, 1}; // world position of polygon
OBJECT4DV1 obj;                     // used to hold our cube mesh
    static float x_ang = 0, y_ang = 0, z_ang = 0;
USHORT(*RGB16Bit)
(int r, int g, int b) = nullptr;
void ClearFrame();

void Build_Sin_Cos_Tables(void);

void Build_Sin_Cos_Tables(void)
{
    for (int ang = 0; ang <= 360; ang++)
    {
        float theta = (float)ang * PI / (float)180;
        cos_look[ang] = cos(theta);
        sin_look[ang] = sin(theta);
    }
}

void ClearFrame()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //        glOrtho(0.0, 500.0, 0.0, 500.0, -1, 1); //设置正射投影的剪裁空间
    //        glOrtho2D(0.0, 500.0, 50 0.0, 0.0);

    //        glOrtho2D(0.0, 500.0, 500.0, 0.0);
    //        glOrtho2D(0.0, 500.0, 0.0, 500.0);
    glOrtho(0.0, display_w, 0.0, display_h, 0.0, 1.0);
    //    glOrtho(0.0,display_w,0.0,display_h,0.0,1.0);

    glBegin(GL_POINTS);
    float delta = (float)1.0 / 255;

    for (int i = 0; i < display_w; i++)
    {
        glColor3f(0, 0, 0);
        for (int j = 0; j < display_h; j++)
        {
            glVertex3f(i, j, 0);
        }
    }
    glEnd();
}

void GameMain();
void GameInit();

void GameInit()
{
    RGB16Bit = RGB16Bit565;

    Build_Sin_Cos_Tables();
    poly1.state = POLY4DV1_STATE_ACTIVE;
    poly1.attr = 0;
    poly1.color = RGB16Bit(0, 255, 0);

    poly1.vlist[0].x = 0;
    poly1.vlist[0].y = 50;
    poly1.vlist[0].z = 0;
    poly1.vlist[0].w = 1;

    poly1.vlist[1].x = 50;
    poly1.vlist[1].y = -50;
    poly1.vlist[1].z = 0;
    poly1.vlist[1].w = 1;

    poly1.vlist[2].x = -50;
    poly1.vlist[2].y = -50;
    poly1.vlist[2].z = 0;
    poly1.vlist[2].w = 1;

    poly1.next = poly1.prev = NULL;

    // initialize the camera with 90 FOV, normalized coordinates
    Init_CAM4DV1(&cam,            // the camera object
                 CAM_MODEL_EULER, // euler camera model
                 &cam_pos,        // initial camera position
                 &cam_dir,        // initial camera angles
                 NULL,            // no initial target
                 50.0,            // near and far clipping planes
                 500.0,
                 90.0,         // field of view in degrees
                 WINDOW_WIDTH, // size of final screen viewport
                 WINDOW_HEIGHT);
}
void GameMain()
{
    int greyLevel = 255;
    // for (int i = 0; i < display_w; i++)
    for (int i = 0; i < display_h; i++)
    {
        // for (int j = 0; j < display_h; j++)
        for (int j = 0; j < display_w; j++)
        {
            // RGB2Color(framebuffer[i * display_w + j], greyLevel * clear_color.x, greyLevel * clear_color.y, greyLevel * clear_color.z);
            RGB2Color(framebuffer[i * display_w + j], greyLevel * clear_color.x, greyLevel * clear_color.y, greyLevel * clear_color.z);
        }
    }

    for (int i = 0; i < 100; i++)
    {
        for (int j = 0; j < 200; j++)
        {
            RGB2Color(framebuffer[i * display_w + j], 255 * clear_color.x, 255 * clear_color.y, 255 * clear_color.z);
        }
    }
    static MATRIX4X4 gMatrixRotate; // general rotation matrix

    static float ang_y = 0; // rotation angle

    Reset_RENDERLIST4DV1(&rend_list);
    Insert_POLYF4DV1_RENDERLIST4DV1(&rend_list, &poly1);
    Build_XYZ_Rotation_MATRIX4X4(0, ang_y, 0, &gMatrixRotate);

    ang_y += 10;
    if (ang_y >= 360.0)
        ang_y = 0;

    Transform_RENDERLIST4DV1(&rend_list, &gMatrixRotate, TRANSFORM_LOCAL_ONLY);
    Model_To_World_RENDERLIST4DV1(&rend_list, &poly1_pos);
    Build_CAM4DV1_Matrix_Euler(&cam, CAM_ROT_SEQ_ZYX);
    World_To_Camera_RENDERLIST4DV1(&rend_list, &cam);
    Camera_To_Perspective_RENDERLIST4DV1(&rend_list, &cam);
    Perspective_To_Screen_RENDERLIST4DV1(&rend_list, &cam);

    RENDERLIST4DV1_PTR rend_list_ptr = &rend_list;

    for (int idx_poly = 0; idx_poly < rend_list_ptr->num_polys; idx_poly++)
    {
        // std::cout << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[0].x << std::endl;
        // std::cout << "x1 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[0].x << std::endl;
        // std::cout << "y1 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[0].y << std::endl;
        // std::cout << "x2 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[1].x << std::endl;
        // std::cout << "y2 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[1].y << std::endl;
        // std::cout << "x3 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[2].x << std::endl;
        // std::cout << "y3 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[2].y << std::endl;
        // std::cout << "-----" << std::endl;

        float x1 = rend_list_ptr->poly_ptrs[idx_poly]->tvlist[0].x;
        float y1 = rend_list_ptr->poly_ptrs[idx_poly]->tvlist[0].y;
        float x2 = rend_list_ptr->poly_ptrs[idx_poly]->tvlist[1].x;
        float y2 = rend_list_ptr->poly_ptrs[idx_poly]->tvlist[1].y;
        float x3 = rend_list_ptr->poly_ptrs[idx_poly]->tvlist[2].x;
        float y3 = rend_list_ptr->poly_ptrs[idx_poly]->tvlist[2].y;

        int c;
        RGB2Color(c, 255, 255, 255);

        device_draw_line(framebuffer,x1, y1, x2, y2, c); //3 1
        device_draw_line(framebuffer,x1, y1, x3, y3, c); //3 1
        device_draw_line(framebuffer,x2, y2, x3, y3, c); //3 1

    } // end for poly
}

void DrawFrame();
void DrawFrame()
{

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //        glOrtho(0.0, 500.0, 0.0, 500.0, -1, 1); //设置正射投影的剪裁空间
    //        glOrtho2D(0.0, 500.0, 50 0.0, 0.0);

    //        glOrtho2D(0.0, 500.0, 500.0, 0.0);
    //        glOrtho2D(0.0, 500.0, 0.0, 500.0);
    glOrtho(0.0, display_w, 0.0, display_h, 0.0, 1.0);
    //    glOrtho(0.0,display_w,0.0,display_h,0.0,1.0);

    glBegin(GL_POINTS);
    float delta = (float)1.0 / 255;

    for (int i = 0; i < display_h; i++)
    {
        // int c = 0;
        // int R = 255*clear_color.x;
        // int G = 255*clear_color.y;
        // int B = 255*clear_color.z;
        //
        // cout<<"R = "<<R<<"G = "<<G<<"b = "<<B;
        // c = (R<<16) | (G<<8) | B;
        // RGB2Color(c, R, G, B);

        // int r = (0xff << 16 & c) >> 16;
        // int g = (0xff << 8 & c) >> 8;
        // int b = 0xff & c;

        int r, g, b;
        // Color2RGB(c, r, g, b);
        for (int j = 0; j < display_w; j++)
        {
            Color2RGB(framebuffer[i * display_w + j], r, g, b);
            // cout<<"r = "<<r<<"g = "<<g<<"b = "<<b;
            // glColor3f((float)r/255,(float)g/255,(float)b/255);
            glColor3f((float)r / 255, (float)g / 255, (float)b / 255);
            // glVertex3f(i,j,0);
            glVertex3f(j, display_h - i, 0);
        }
    }

    glEnd();
}

void InitDemo9_2();
void DrawDemo9_2();

OBJECT4DV2 obj_constant_water,
    obj_flat_water,
    obj_gouraud_water,
    obj_constant_light;

RGBAV1 white, gray, black, red, green, blue;
void InitDemo9_2()
{
    Build_Sin_Cos_Tables();
    // POINT4D cam_pos = {0, 0, 0, 1};
    POINT4D cam_pos = {30, 0, 0, 1};
    POINT4D cam_target = {0, 0, 0, 1};
    VECTOR4D cam_dir = {0, -42, 0, 1};

    // all your initialization code goes here...
    VECTOR4D vscale = {1.0, 1.0, 1.0, 1},
             vpos = {0, 0, 0, 1},
             vrot = {0, 0, 0, 1};

    Init_CAM4DV1(&cam,            // the camera object
                 CAM_MODEL_EULER, // the euler model
                 &cam_pos,        // initial camera position
                 &cam_dir,        // initial camera angles
                 &cam_target,     // no target
                 200.0,           // near and far clipping planes
                 12000.0,
                 120.0,        // field of view in degrees
                 WINDOW_WIDTH, // size of final screen viewport
                 WINDOW_HEIGHT);

    VECTOR4D_INITXYZ(&vscale, 20.00, 20.00, 20.00);

    // Load_OBJECT4DV2_COB(&obj_constant_water, "./cob/water_constant_01.cob",
    // Load_OBJECT4DV2_COB(&obj_constant_water, "./cob/cube_constant_01.cob",						&vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);
    Load_OBJECT4DV2_COB(&obj_constant_water, "./cob/water_flat_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);

    // load flat shaded water
    VECTOR4D_INITXYZ(&vscale, 20.00, 20.00, 20.00);
    // Load_OBJECT4DV2_COB(&obj_flat_water, "./cob/water_gouraud_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD); //修改颜色后的水分子
    Load_OBJECT4DV2_COB(&obj_flat_water, "./cob/water_flat_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);

    // load gouraud shaded water
    VECTOR4D_INITXYZ(&vscale, 20.00, 20.00, 20.00);
    // Load_OBJECT4DV2_COB(&obj_gouraud_water, "./cob/water_flat_01_gouraud.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD); //修改颜色后的水分子
    Load_OBJECT4DV2_COB(&obj_gouraud_water, "./cob/water_flat_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);
                                                                                                                                                                                                 //    Load_OBJECT4DV2_COB(&obj_gouraud_water, "./cob/water_gouraud_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD); //修改颜色后的水分子

    VECTOR4D_INITXYZ(&vscale, 5.00, 5.00, 5.00);

    Load_OBJECT4DV2_COB(&obj_constant_light, "./cob/cube_flat_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);
    Reset_Lights_LIGHTV1();

    // create some working colors
    white.rgba = _RGBA32BIT(255, 255, 255, 0);
    gray.rgba = _RGBA32BIT(100, 100, 100, 0);
    black.rgba = _RGBA32BIT(0, 0, 0, 0);
    red.rgba = _RGBA32BIT(255, 0, 0, 0);
    green.rgba = _RGBA32BIT(0, 255, 0, 0);
    blue.rgba = _RGBA32BIT(0, 0, 255, 0);

    // ambient light
    Init_Light_LIGHTV1(AMBIENT_LIGHT_INDEX,
                       LIGHTV1_STATE_ON,     // turn the light on
                       LIGHTV1_ATTR_AMBIENT, // ambient light type
                                             //   gray, black, black,   // color for ambient term only
                       black, black, black,  // color for ambient term only
                       NULL, NULL,           // no need for pos or dir
                       0, 0, 0,              // no need for attenuation
                       0, 0, 0);             // spotlight info NA

    VECTOR4D dlight_dir = {-1, 0, -1, 0};

    // directional light
    Init_Light_LIGHTV1(INFINITE_LIGHT_INDEX,
                       LIGHTV1_STATE_ON,      // turn the light on
                       LIGHTV1_ATTR_INFINITE, // infinite light type
                                              //   black, gray, black,	// color for diffuse term only
                       gray, gray, gray,      // color for diffuse term only
                       NULL, &dlight_dir,     // need direction only
                       0, 0, 0,               // no need for attenuation
                       0, 0, 0);              // spotlight info NA

    VECTOR4D plight_pos = {0, 200, 0, 0};

    // point light
    Init_Light_LIGHTV1(POINT_LIGHT_INDEX,
                       LIGHTV1_STATE_ON,    // turn the light on
                       LIGHTV1_ATTR_POINT,  // pointlight type
                                            //   black, green, black, // color for diffuse term only
                       black, white, white, // color for diffuse term only
                       &plight_pos, NULL,   // need pos only
                       0, .007, 0,          // linear attenuation only
                       0, 0, 1);            // spotlight info NA

    VECTOR4D slight2_pos = {0, 200, 0, 0};
    VECTOR4D slight2_dir = {-1, 0, -1, 0};

    // spot light2
    //    Init_Light_LIGHTV1(SPOT_LIGHT2_INDEX,
    // 					  LIGHTV1_STATE_ON,			  // turn the light on
    // 					  LIGHTV1_ATTR_SPOTLIGHT2,	// spot light type 2
    // 					//   black, red, black,		  // color for diffuse term only
    // 					  black, white, black,		  // color for diffuse term only
    // 					  &slight2_pos, &slight2_dir, // need pos only
    // 					  0, .001, 0,				  // linear attenuation only
    // 					  0, 0, 1);
}
void DrawDemo9_2()
{
    // Sleep(20);

    static MATRIX4X4 mrot; // general rotation matrix

    static float view_angle = 0;
    static float camera_distance = 6000;
    static VECTOR4D pos = {0, 0, 0, 0};
    static float tank_speed = 10;
    static float turning = 0;

    char work_string[256]; // temp string

    int index; // looping var

    static float plight_ang = 0, slight_ang = 0; // angles for light motion

    // move point light source in ellipse around game world
    // lights[POINT_LIGHT_INDEX].pos.x = 1000 * Fast_Cos(plight_ang);
    // lights[POINT_LIGHT_INDEX].pos.x = 200;
    // lights[POINT_LIGHT_INDEX].pos.y = 10;
    // lights[POINT_LIGHT_INDEX].pos.z = 1000 * Fast_Sin(plight_ang);

    // if ((plight_ang += 3) > 360)
    // 	plight_ang = 0;

    // move spot light source in ellipse around game world
    lights[SPOT_LIGHT2_INDEX].pos.x = 1000 * Fast_Cos(slight_ang);
    lights[SPOT_LIGHT2_INDEX].pos.y = 200;
    lights[SPOT_LIGHT2_INDEX].pos.z = 1000 * Fast_Sin(slight_ang);

    // if ((slight_ang -= 5) < 0)
    // 	slight_ang = 360;
    //Reset_RENDERLIST4DV1(&rend_list);
    Reset_RENDERLIST4DV2(&rend_list2);

    // if (KEY_DOWN(VK_NUMPAD8))
    // {
    //     // move forward
    //     cam.pos.y += 1;
    // } // end if

    // if (KEY_DOWN(VK_NUMPAD2))
    // {
    //     cam.pos.y -= 1;
    // } // end if

    // if (KEY_DOWN(VK_NUMPAD9))
    // {
    //     lights[POINT_LIGHT_INDEX].pos.z += 10;
    //     // move forward
    // } // end if

    // if (KEY_DOWN(VK_NUMPAD3))
    // {
    //     lights[POINT_LIGHT_INDEX].pos.z -= 10;
    // } // end if

    // if (KEY_DOWN(VK_NUMPAD7))
    // {
    //     lights[POINT_LIGHT_INDEX].pos.y += 10;
    //     // move forward
    // } // end if

    // if (KEY_DOWN(VK_NUMPAD1))
    // {
    //     lights[POINT_LIGHT_INDEX].pos.y -= 10;
    // } // end if

    // if (KEY_DOWN(VK_NUMPAD4))
    // {
    //     lights[POINT_LIGHT_INDEX].pos.x -= 10;
    //     // move forward
    // } // end if

    // if (KEY_DOWN(VK_NUMPAD6))
    // {
    //     lights[POINT_LIGHT_INDEX].pos.x += 10;
    // } // end if

    // if (KEY_DOWN(VK_UP))
    // {
    //     // move forward
    //     cam.pos.x += tank_speed * Fast_Sin(cam.dir.y);

    //     cam.pos.z += tank_speed * Fast_Cos(cam.dir.y);
    // } // end if

    // if (KEY_DOWN(VK_DOWN))
    // {
    //     // move backward
    //     cam.pos.x -= tank_speed * Fast_Sin(cam.dir.y);
    //     cam.pos.z -= tank_speed * Fast_Cos(cam.dir.y);
    // } // end if

    // // rotate
    // if (KEY_DOWN(VK_RIGHT))
    // {
    //     cam.dir.y += 3;

    //     // add a little turn to object
    //     if ((turning += 2) > 15)
    //         turning = 15;

    // } // end if

    static bool gF1 = false;
    static bool gF2 = false;

    // if (KEY_DOWN(VK_F1))
    // {
    //     gF1 = !gF1;
    // } // end if

    // if (KEY_DOWN(VK_F2))
    // {
    //     gF2 = !gF2;
    // } // end if

    // if (KEY_DOWN(VK_LEFT))
    // {
    //     cam.dir.y -= 3;

    //     // add a little turn to object
    //     if ((turning -= 2) < -15)
    //         turning = -15;

    // }    // end if
    // else // center heading again
    // {
    //     if (turning > 0)
    //         turning -= 1;
    //     else if (turning < 0)
    //         turning += 1;

    // } // end else
    Build_CAM4DV1_Matrix_Euler(&cam, CAM_ROT_SEQ_ZYX);

    Reset_OBJECT4DV2(&obj_constant_water);
    obj_constant_water.world_pos.x = -50;
    obj_constant_water.world_pos.y = 0;
    obj_constant_water.world_pos.z = 120;
    Build_XYZ_Rotation_MATRIX4X4(x_ang, y_ang, z_ang, &mrot);
    Transform_OBJECT4DV2(&obj_constant_water, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);
    Model_To_World_OBJECT4DV2(&obj_constant_water, TRANSFORM_TRANS_ONLY);
    Insert_OBJECT4DV2_RENDERLIST4DV2(&rend_list2, &obj_constant_water, 0);

    Reset_OBJECT4DV2(&obj_constant_light);
    obj_constant_light.world_pos.x = lights[POINT_LIGHT_INDEX].pos.x;
    obj_constant_light.world_pos.y = lights[POINT_LIGHT_INDEX].pos.y;
    obj_constant_light.world_pos.z = lights[POINT_LIGHT_INDEX].pos.z;
    Build_XYZ_Rotation_MATRIX4X4(x_ang, y_ang, z_ang, &mrot);
    Transform_OBJECT4DV2(&obj_constant_light, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);
    Model_To_World_OBJECT4DV2(&obj_constant_light, TRANSFORM_TRANS_ONLY);
    // Insert_OBJECT4DV2_RENDERLIST4DV2(&rend_list2, &obj_constant_light, 0);

    Reset_OBJECT4DV2(&obj_flat_water);
    obj_flat_water.world_pos.x = 0;
    obj_flat_water.world_pos.y = 0;
    obj_flat_water.world_pos.z = 120;
    Build_XYZ_Rotation_MATRIX4X4(x_ang, y_ang, z_ang, &mrot);
    Transform_OBJECT4DV2(&obj_flat_water, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);
    Model_To_World_OBJECT4DV2(&obj_flat_water, TRANSFORM_TRANS_ONLY);
    // Insert_OBJECT4DV2_RENDERLIST4DV2(&rend_list2, &obj_flat_water, 0);

    Reset_OBJECT4DV2(&obj_gouraud_water);
    obj_gouraud_water.world_pos.x = 50;
    obj_gouraud_water.world_pos.y = 0;
    obj_gouraud_water.world_pos.z = 120;
    Build_XYZ_Rotation_MATRIX4X4(x_ang, y_ang, z_ang, &mrot);
    Transform_OBJECT4DV2(&obj_gouraud_water, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);
    Model_To_World_OBJECT4DV2(&obj_gouraud_water, TRANSFORM_TRANS_ONLY);
    // Insert_OBJECT4DV2_RENDERLIST4DV2(&rend_list2, &obj_gouraud_water, 0);

    RENDERLIST4DV2_PTR rend_list_ptr = &rend_list2;
    // update rotation angles
    // if ((x_ang += 1) > 360)
    // 	x_ang = 0;
    // if (gF1)
    if (true)
    {
        if ((y_ang += 2) > 360)
            y_ang = 0;
    }
    // if ((z_ang += 3) > 360)
    // 	z_ang = 0;

    // remove backfaces
    Remove_Backfaces_RENDERLIST4DV2(&rend_list2, &cam);

    // light scene all at once
    Light_RENDERLIST4DV2_World16(&rend_list2, &cam, lights, 4);

    // apply world to camera transform
    World_To_Camera_RENDERLIST4DV2(&rend_list2, &cam);

    // sort the polygon list (hurry up!)
    Sort_RENDERLIST4DV2(&rend_list2, SORT_POLYLIST_AVGZ);

    // apply camera to perspective transformation
    Camera_To_Perspective_RENDERLIST4DV2(&rend_list2, &cam);

    // apply screen transform
    Perspective_To_Screen_RENDERLIST4DV2(&rend_list2, &cam);

    //  RENDERLIST4DV1_PTR rend_list_ptr = &rend_list;

    POLYF4DV2 face; // temp face used to render polygon

    for (int poly = 0; poly < rend_list_ptr->num_polys; poly++)
    {
        if (!(rend_list_ptr->poly_ptrs[poly]->state & POLY4DV1_STATE_ACTIVE) || (rend_list_ptr->poly_ptrs[poly]->state & POLY4DV1_STATE_CLIPPED) || (rend_list_ptr->poly_ptrs[poly]->state & POLY4DV1_STATE_BACKFACE))
            continue;
        float x1 = rend_list_ptr->poly_ptrs[poly]->tvlist[0].x;
        float y1 = rend_list_ptr->poly_ptrs[poly]->tvlist[0].y;
        float x2 = rend_list_ptr->poly_ptrs[poly]->tvlist[1].x;
        float y2 = rend_list_ptr->poly_ptrs[poly]->tvlist[1].y;
        float x3 = rend_list_ptr->poly_ptrs[poly]->tvlist[2].x;
        float y3 = rend_list_ptr->poly_ptrs[poly]->tvlist[2].y;

        int color = rend_list_ptr->poly_ptrs[poly]->lit_color[0];
        if ((rend_list_ptr->poly_ptrs[poly]->attr & POLY4DV2_ATTR_SHADE_MODE_FLAT) ||
            (rend_list_ptr->poly_ptrs[poly]->attr & POLY4DV2_ATTR_SHADE_MODE_CONSTANT))
        {

            IUINT32 c = (255 << 16) | (255 << 8) | 255;

            //  std::cout<<"color = "<<color<<std::endl;
            //  int color = rend_list_ptr->poly_ptrs[poly]->color;
            //  rend_list->poly_ptrs[poly]->lit_color[0]
            // DrawTrianglePureColor(&device, x1, y1, x2, y2, x3, y3, color);
            // DrawTrianglePureColor2(&device, x1, y1, x2, y2, x3, y3, color);
            DrawTrianglePureColor2(framebuffer, x1, y1, x2, y2, x3, y3, color);
            //  device_draw_line(&device, x1, y1, x2, y2, c);
            //  device_draw_line(&device, x1, y1, x3, y3, c);
            //  device_draw_line(&device, x2, y2, x3, y3, c);
        }
        else if (rend_list_ptr->poly_ptrs[poly]->attr & POLY4DV2_ATTR_SHADE_MODE_GOURAUD)
        {
            // {andre take advantage of the data structures later..}
            // set the vertices
            face.tvlist[0].x = (int)rend_list_ptr->poly_ptrs[poly]->tvlist[0].x;
            face.tvlist[0].y = (int)rend_list_ptr->poly_ptrs[poly]->tvlist[0].y;
            face.lit_color[0] = rend_list_ptr->poly_ptrs[poly]->lit_color[0];

            face.tvlist[1].x = (int)rend_list_ptr->poly_ptrs[poly]->tvlist[1].x;
            face.tvlist[1].y = (int)rend_list_ptr->poly_ptrs[poly]->tvlist[1].y;
            face.lit_color[1] = rend_list_ptr->poly_ptrs[poly]->lit_color[1];

            face.tvlist[2].x = (int)rend_list_ptr->poly_ptrs[poly]->tvlist[2].x;
            face.tvlist[2].y = (int)rend_list_ptr->poly_ptrs[poly]->tvlist[2].y;
            face.lit_color[2] = rend_list_ptr->poly_ptrs[poly]->lit_color[2];

            // draw the gouraud shaded triangle
            // Draw_Gouraud_Triangle16(&face, video_buffer, lpitch);

            // Draw_Gouraud_Triangle16(&device, &face);

            // DrawPhongTriangle(&device, &cam, &face, rend_list_ptr->poly_ptrs[poly], lights, gF2);

            //     std::cout<<"...tu0 = "<< tu0    <<"tv0 = "<< tv0     <<"tw0 = "<< tw0
            //  <<"tu1 = "<< tu1     <<"tv1 = "<< tv1     <<"tw1 = "<< tw1
            //  <<"tu2 = "<< tu2     <<"tv2 = "<< tv2     <<"tw2 = "<< tw2<< std::endl;

            //          std::cout<<"...tpu0 = "<< tpu0    <<"tpv0 = "<< tpv0     <<"tpw0 = "<< tpw0
            //  <<"tpu1 = "<< tpu1     <<"tpv1 = "<< tpv1     <<"tpw1 = "<< tpw1
            //  <<"tpu2 = "<< tpu2     <<"tpv2 = "<< tpv2     <<"tpw2 = "<< tpw2<< std::endl;

            // DrawTrianglePureColor2(&device, x1, y1, x2, y2, x3, y3, color);
        } // end if gouraud
    }
}

int main(int, char **)
{


    // Setup SDL
    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) != 0)
    {
        printf("Error: %s\n", SDL_GetError());
        return -1;
    }

    // Setup window
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
    SDL_DisplayMode current;
    SDL_GetCurrentDisplayMode(0, &current);
    SDL_Window *window = SDL_CreateWindow("Dear ImGui SDL2+OpenGL example", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1280, 720, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
    SDL_GLContext gl_context = SDL_GL_CreateContext(window);
    SDL_GL_SetSwapInterval(1); // Enable vsync

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    (void)io;

    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsClassic();

    // Setup Platform/Renderer bindings
    ImGui_ImplSDL2_InitForOpenGL(window, gl_context);
    ImGui_ImplOpenGL2_Init();

    // Load Fonts
    // - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
    // - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
    // - If the file cannot be loaded, the function will return NULL. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
    // - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which ImGui_ImplXXXX_NewFrame below will call.
    // - Read 'misc/fonts/README.txt' for more instructions and details.
    // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
    //io.Fonts->AddFontDefault();
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/ProggyTiny.ttf", 10.0f); //ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, NULL, io.Fonts->GetGlyphRangesJapanese());
    //IM_ASSERT(font != NULL);

    bool show_demo_window = false;
    bool show_another_window = false;

    //    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    // Main loop
    bool done = false;
    bool gInited = false;
    while (!done)
    {

        // Poll and handle events (inputs, window resize, etc.)
        // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
        // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application.
        // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application.
        // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            ImGui_ImplSDL2_ProcessEvent(&event);
            if (event.type == SDL_QUIT)
                done = true;
        }

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplSDL2_NewFrame(window);
        ImGui::NewFrame();

        // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
        if (show_demo_window)
            ImGui::ShowDemoWindow(&show_demo_window);

        // 2. Show a simple window that we create ourselves. We use a Begin/End pair to created a named window.
        {
            static float f = 0.0f;
            static int counter = 0;

            ImGui::Begin("Hello, world!"); // Create a window called "Hello, world!" and append into it.

            ImGui::Text("This is some useful text.");          // Display some text (you can use a format strings too)
            ImGui::Checkbox("Demo Window", &show_demo_window); // Edit bools storing our window open/close state
            ImGui::Checkbox("Another Window", &show_another_window);

            ImGui::SliderFloat("cam.pos.z", &cam.pos.z, -100.0f, 100.0f);             
            ImGui::SliderFloat("cam.pos.y", &cam.pos.y, -100.0f, 100.0f);             
            ImGui::SliderFloat("cam.pos.x", &cam.pos.x, -100.0f, 100.0f);             
            ImGui::SliderFloat("y_ang", &y_ang, 0.0f, 360.0f);             
            ImGui::ColorEdit3("clear color", (float *)&clear_color); // Edit 3 floats representing a color

            if (ImGui::Button("Button")) // Buttons return true when clicked (most widgets return true when edited/activated)
                counter++;
            ImGui::SameLine();
            ImGui::Text("counter = %d", counter);

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            ImGui::End();
        }

        // 3. Show another simple window.
        if (show_another_window)
        {
            ImGui::Begin("Another Window", &show_another_window); // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
            ImGui::Text("Hello from another window!");
            if (ImGui::Button("Close Me"))
                show_another_window = false;
            ImGui::End();
        }

        // Rendering
        ImGui::Render();
        glViewport(0, 0, (int)io.DisplaySize.x, (int)io.DisplaySize.y);
        if (!gInited)
        {
            display_w = (int)io.DisplaySize.x;
            display_h = (int)io.DisplaySize.y;
            framebuffer = new int[display_w * display_h];
            InitDemo9_2();
            gInited = true;
        }
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);

        //清屏
        int greyLevel = 255;
        // for (int i = 0; i < display_w; i++)
        for (int i = 0; i < display_h; i++)
        {
            // for (int j = 0; j < display_h; j++)
            for (int j = 0; j < display_w; j++)
            {
                // RGB2Color(framebuffer[i * display_w + j], greyLevel * clear_color.x, greyLevel * clear_color.y, greyLevel * clear_color.z);
                RGB2Color(framebuffer[i * display_w + j], greyLevel * clear_color.x, greyLevel * clear_color.y, greyLevel * clear_color.z);
            }
        }
        DrawDemo9_2();


        DrawFrame();

        //glUseProgram(0); // You may want this if using this code in an OpenGL 3+ context where shaders may be bound
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
        SDL_GL_SwapWindow(window);
    }

    // Cleanup
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImGui::DestroyContext();

    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}

