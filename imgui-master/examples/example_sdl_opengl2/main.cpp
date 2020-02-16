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


// #include "sphere.h"
// #include "hitable_list.h"
// #include "float.h"
// #include "camera.h"
// #include "material.h"
// #include "random.h"
#include <thread>
//全局变量
enum RenderType
{
    Wireframe,
    PureTriangle
};

static int gRenderType = 0;
int *framebuffer; //帧缓冲
ImVec4 clear_color = ImVec4(0.0f, 0.0f, 0.00f, 1.00f);
int display_w, display_h = 2000;

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
// OBJECT4DV1 obj;                     // used to hold our cube mesh
static float x_ang = 0, y_ang = 0, z_ang = 0;
USHORT(*RGB16Bit)
(int r, int g, int b) = nullptr;
void ClearFrame();

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

    // glMatrixMode(GL_PROJECTION);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //        glOrtho(0.0, 500.0, 0.0, 500.0, -1, 1); //设置正射投影的剪裁空间
    //        glOrtho2D(0.0, 500.0, 50 0.0, 0.0);

    //        glOrtho2D(0.0, 500.0, 500.0, 0.0);
    //        glOrtho2D(0.0, 500.0, 0.0, 500.0);
    // glOrtho(0.0, display_w, 0.0, display_h, 0.0, 1.0);
    // glOrtho(0.0, display_h, 0.0, display_w, 0.0, 1.0);
    glOrtho(0.0, display_w, 0.0, display_h, 0.0, 1.0);

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

        int r, g, b = 255;
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

    //water_flat_01.cob ok

    Load_OBJECT4DV2_COB(&obj_constant_water, "./cob/water_flat_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);
    // Load_OBJECT4DV2_COB(&obj_constant_water, "./cob/water_gouraud_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);
    
                    // load flat shaded water
                    VECTOR4D_INITXYZ(&vscale, 20.00, 20.00, 20.00);
    // Load_OBJECT4DV2_COB(&obj_flat_water, "./cob/water_gouraud_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD); //修改颜色后的水分子
    Load_OBJECT4DV2_COB(&obj_flat_water, "./cob/water_flat_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);

    // load gouraud shaded water
    VECTOR4D_INITXYZ(&vscale, 20.00, 20.00, 20.00);
    // Load_OBJECT4DV2_COB(&obj_gouraud_water, "./cob/water_flat_01_gouraud.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD); //修改颜色后的水分子
    Load_OBJECT4DV2_COB(&obj_gouraud_water, "./cob/water_flat_01_gouraud.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);
                                                                                                                                                                                                 //    Load_OBJECT4DV2_COB(&obj_gouraud_water, "./cob/water_gouraud_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD); //修改颜色后的水分子

    VECTOR4D_INITXYZ(&vscale, 5.00, 5.00, 5.00);

    Load_OBJECT4DV2_COB(&obj_constant_light, "./cob/water_gouraud_01.cob", &vscale, &vpos, &vrot, VERTEX_FLAGS_SWAP_YZ | VERTEX_FLAGS_TRANSFORM_LOCAL | VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD);
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
    obj_constant_water.world_pos.x = -150;
    obj_constant_water.world_pos.y = 20;
    obj_constant_water.world_pos.z = 120;
    Build_XYZ_Rotation_MATRIX4X4(x_ang, y_ang, z_ang, &mrot);
    Transform_OBJECT4DV2(&obj_constant_water, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);
    Model_To_World_OBJECT4DV2(&obj_constant_water, TRANSFORM_TRANS_ONLY);
    Insert_OBJECT4DV2_RENDERLIST4DV2(&rend_list2, &obj_constant_water, 0);

    // Reset_OBJECT4DV2(&obj_constant_light);
    // obj_constant_light.world_pos.x = lights[POINT_LIGHT_INDEX].pos.x;
    // obj_constant_light.world_pos.y = lights[POINT_LIGHT_INDEX].pos.y;
    // obj_constant_light.world_pos.z = lights[POINT_LIGHT_INDEX].pos.z;
    // Build_XYZ_Rotation_MATRIX4X4(x_ang, y_ang, z_ang, &mrot);
    // Transform_OBJECT4DV2(&obj_constant_light, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);
    // Model_To_World_OBJECT4DV2(&obj_constant_light, TRANSFORM_TRANS_ONLY);
    // Insert_OBJECT4DV2_RENDERLIST4DV2(&rend_list2, &obj_constant_light, 0);

    // Reset_OBJECT4DV2(&obj_flat_water);
    // obj_flat_water.world_pos.x = 0;
    // obj_flat_water.world_pos.y = 0;
    // obj_flat_water.world_pos.z = 120;
    // Build_XYZ_Rotation_MATRIX4X4(x_ang, y_ang, z_ang, &mrot);
    // Transform_OBJECT4DV2(&obj_flat_water, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);
    // Model_To_World_OBJECT4DV2(&obj_flat_water, TRANSFORM_TRANS_ONLY);
    // Insert_OBJECT4DV2_RENDERLIST4DV2(&rend_list2, &obj_flat_water, 0);

    Reset_OBJECT4DV2(&obj_gouraud_water);
    obj_gouraud_water.world_pos.x = -100;
    obj_gouraud_water.world_pos.y = 0;
    obj_gouraud_water.world_pos.z = 120;
    Build_XYZ_Rotation_MATRIX4X4(x_ang, y_ang, z_ang, &mrot);
    Transform_OBJECT4DV2(&obj_gouraud_water, &mrot, TRANSFORM_LOCAL_TO_TRANS, 1);
    Model_To_World_OBJECT4DV2(&obj_gouraud_water, TRANSFORM_TRANS_ONLY);
    Insert_OBJECT4DV2_RENDERLIST4DV2(&rend_list2, &obj_gouraud_water, 0);

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

    POLYF4DV2 face; // temp face used to render polygon

    // cout<<rend_list_ptr->num_polys<<endl;
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
        if (gRenderType == Wireframe)
        {
            int c;
            RGB2Color(c, 255, 255, 255);
            DrawTriangleWireframe(framebuffer, x1, y1, x2, y2, x3, y3, c);
        }
        else if (gRenderType == PureTriangle)
        {

            if ((rend_list_ptr->poly_ptrs[poly]->attr & POLY4DV2_ATTR_SHADE_MODE_FLAT) || (rend_list_ptr->poly_ptrs[poly]->attr & POLY4DV2_ATTR_SHADE_MODE_CONSTANT))
            {
                DrawTrianglePureColor2(framebuffer, x1, y1, x2, y2, x3, y3, color);
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

                Draw_Gouraud_Triangle16(framebuffer, &face);
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
}

// vec3 color(const ray &r, hitable *world, int depth);
// hitable *random_scene();

// vec3 color(const ray& r, hitable *world, int depth) {
    
//     hit_record rec;
//     if (world->hit(r, 0.001, MAXFLOAT, rec)) {
//         ray scattered;
//         vec3 attenuation;
//         if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//              return attenuation*color(scattered, world, depth+1);
//         }
//         else {
//             return vec3(0,0,0);
//         }
//     }
//     else {
//         vec3 unit_direction = unit_vector(r.direction());
//         float t = 0.5*(unit_direction.y() + 1.0);
//         return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
//     }
// }


// hitable *random_scene() {
//     int n = 500;
//     // int n = 50;
//     hitable **list = new hitable*[n+1];
//     list[0] =  new sphere(vec3(0,-1000,0), 1000, new lambertian(vec3(0.5, 0.5, 0.5)));
//     int i = 1;
//     // for (int a = -11; a < 11; a++) {
//     for (int a = -2; a < 2; a++) {
//         // for (int b = -11; b < 11; b++) {
//         for (int b = -2; b < 2; b++) {
//             float choose_mat = random_double();
//             vec3 center(a+0.9*random_double(),0.2,b+0.9*random_double());
//             if ((center-vec3(4,0.2,0)).length() > 0.9) {
//                 if (choose_mat < 0.8) {  // diffuse
//                     list[i++] = new sphere(
//                         center, 0.2,
//                         new lambertian(vec3(random_double()*random_double(),
//                                             random_double()*random_double(),
//                                             random_double()*random_double()))
//                     );
//                 }
//                 else if (choose_mat < 0.95) { // metal
//                     list[i++] = new sphere(
//                         center, 0.2,
//                         new metal(vec3(0.5*(1 + random_double()),
//                                        0.5*(1 + random_double()),
//                                        0.5*(1 + random_double())),
//                                   0.5*random_double())
//                     );
//                 }
//                 else {  // glass
//                     list[i++] = new sphere(center, 0.2, new dielectric(1.5));
//                 }
//             }
//         }
//     }

//     list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
//     list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
//     list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));
//     return new hitable_list(list,i);

// }

// vec3 color(const ray& r) {
//     vec3 unit_direction = unit_vector(r.direction());
//     float t = 0.5*(unit_direction.y() + 1.0);
//     return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
//     // return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.0, 0.0, 0.0);
//     // return (1.0-t)*vec3(1.0, 0.0, 0.0) + t*vec3(0.0, 1.0, 0.0);
// }

/*
本次改动
1.增加材质抽象类 material
2.增加兰伯特材质
3.增加mental材质
4.修改color函数
*/
#include <unistd.h>
#include "vec3.h"
#include <cstdlib>
inline double random_double() {
    return rand() / (RAND_MAX + 1.0);
}

vec3 random_in_unit_sphere()
{
    vec3 p;
    do
    {
        p = 2.0 * vec3(random_double(), random_double(), random_double()) - vec3(1, 1, 1);
    } while (p.squared_length() >= 1.0);
    return p;
}

float progressIdx = 0.0f, progressDir = 1.0f;

class ray
{
public:
    vec3 A; //起点
    vec3 B; //方向

    ray() {}
    ray(const vec3 &a, const vec3 &b)
    {
        A = a;
        B = b;
    }
    vec3 origin() const { return A; }
    vec3 direction() const { return B; }
    vec3 point_at_parameter(float t) const { return A + t * B; } //终点的坐标

};
struct hit_record;

//attenuation 衰减
class material  {
    public:
        virtual bool scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation, ray &scattered) const = 0;
};


struct hit_record
{
    float t;   //命中射线的长度
    vec3 p;    //命中终点坐标
    vec3 normal; //命中点的法线
    material *mat_ptr;
};


class lambertian : public material
{
public:
    vec3 albedo; //反射率
    lambertian(const vec3 &a) : albedo(a) {}
    virtual bool scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation, ray &scattered) const
    {
        vec3 s_world = rec.p + rec.normal + random_in_unit_sphere();
        scattered = ray(rec.p, s_world - rec.p);
        attenuation = albedo;
        return true;
    }

};

class metal : public material
{
public:
    vec3 albedo;
    float fuzz;

    metal(const vec3 &a, float f) : albedo(a)
    {
        if (f < 1)
            fuzz = f;
        else
            fuzz = 1;
    }
    virtual bool scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation, ray &scattered) const
    {
        vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
        // scattered = ray(rec.p, reflected);
        scattered = ray(rec.p, reflected + fuzz * random_in_unit_sphere());
        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > 0);
    }
};



class hittable
{
public:
    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const = 0;
};

class sphere : public hittable
{
public:
    sphere() {}
    sphere(vec3 cen, float r, material *m) : center(cen), radius(r), mat_ptr(m){};
    vec3 center;
    float radius;
    material *mat_ptr; /* NEW */

    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const
    {
        vec3 oc = r.origin() - center;
        float a = dot(r.direction(), r.direction());
        float b = dot(oc, r.direction());
        float c = dot(oc, oc) - radius * radius;
        float discriminant = b * b - a * c;
        if (discriminant > 0)
        {
            float temp = (-b - sqrt(discriminant)) / a; //小根
            if (temp < t_max && temp > t_min)
            {
                rec.t = temp;
                rec.p = r.point_at_parameter(rec.t);
                rec.normal = (rec.p - center) / radius;
                rec.mat_ptr = mat_ptr; /* NEW */

                return true;
            }
            temp = (-b + sqrt(discriminant)) / a; //大根
            if (temp < t_max && temp > t_min)
            {
                rec.t = temp;
                rec.p = r.point_at_parameter(rec.t);
                rec.normal = (rec.p - center) / radius;
                rec.mat_ptr = mat_ptr; /* NEW */

                return true;
            }
        }
        return false;
    }
};

class hittable_list : public hittable
{
public:
    hittable **list;
    int list_size;

    hittable_list() {}
    hittable_list(hittable **l, int n)
    {
        list = l;
        list_size = n;
    }
    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const
    {

        hit_record temp_rec;
        bool hit_anything = false;
        double closest_so_far = t_max;
        for (int i = 0; i < list_size; i++)
        {
            if (list[i]->hit(r, t_min, closest_so_far, temp_rec))
            {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec; //只记录打到的最近的球
            }
        }
        return hit_anything;
    }
};



vec3 color(const ray &r, hittable *world, int depth)
{
    hit_record rec;
    if (world->hit(r, 0.001, MAXFLOAT, rec)) //射线命中物体
    {
        ray scattered;
        vec3 attenuation;
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            return attenuation*color(scattered, world, depth+1);
        }
        else {
            return vec3(0,0,0);
        }
    }
    else
    {
        vec3 unit_direction = unit_vector(r.direction());
        float t = 0.5 * (unit_direction.y() + 1.0);
        return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
    }
}

class camera
{
public:
    vec3 origin;
    vec3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;

    camera()
    {
        lower_left_corner = vec3(-2.0, -1.0, -1.0);
        horizontal = vec3(4.0, 0.0, 0.0);
        vertical = vec3(0.0, 2.0, 0.0);
        origin = vec3(0.0, 0.0, 0.0);
    }

    ray get_ray(float u, float v)
    {
        return ray(origin, lower_left_corner + u * horizontal + v * vertical - origin);
    }

};

void testThread(int i)
{
            // std::cout << "Hello from lamda thread " << std::this_thread::get_id() << std::endl;
            cout<<"hello -- "<<i<<endl;
}

void foo(const int  &x,char *mychar)
{
	std::cout << &x << "   " << &mychar << std::endl;
	std::cout << "正在运行的线程为：" << std::this_thread::get_id() << "线程的参数为： " << x <<"  "<<mychar<< std::endl;
	return;
}

#include <fstream>
#include <vector>
#include <ctime>
using namespace std;
time_t curTime;
float totTime;
void RenderInThread(int numThread)
{
    // for (int j = (ny - 1) - i; j >= 0; j -= numThread)
    // {
    //     progressIdx = float(ny - 1 - j) / (ny - 1);
    //     for (int i = 0; i < nx; i++)
    //     {

    //         vec3 col(0, 0, 0);
    //         for (int s = 0; s < ns; s++)
    //         {
    //             float u = float(i + random_double()) / float(nx);
    //             float v = float(j + random_double()) / float(ny);
    //             ray r = cam.get_ray(u, v);
    //             col += color(r, world, 0);
    //             // col += color(r, world);
    //         }
    //         col /= float(ns);
    //         col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));

    //         int ir = int(255.99 * col[0]);
    //         int ig = int(255.99 * col[1]);
    //         int ib = int(255.99 * col[2]);

    //         device_pixel(framebuffer, i, ny - 1 - j, ir, ig, ib);
    //         // std::cout << ir << " " << ig << " " << ib << "\n";
    //         outFile << ir << " " << ig << " " << ib << endl;
    //         // outFile << icolor.r << " " << icolor.g << " " << icolor.b << endl;
    //     }
    // }
}
void RayTracing()
{
    usleep(1000);    // will sleep for 1 ms
    usleep(1);       // will sleep for 0.001 ms
    usleep(1000000); // will sleep for 1 s
                     // usleep(1000000); // will sleep for 1 s
                     // int nx = 600;
                     // int ny = 300;
    int ns = 100;

    int nx = 1200;
    int ny = 600;
    // int nx = 800;
    // int ny = 400;
    // int nx = 400;
    // int ny = 200;
    // int nx = 200;
    // int ny = 100;
    std::cout << "P3\n"
              << nx << " " << ny << "\n255\n";

    hittable *list[4];
    list[0] = new sphere(vec3(0, 0, -1), 0.5, new lambertian(vec3(0.8, 0.3, 0.3)));
    list[1] = new sphere(vec3(0, -100.5, -1), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
    // list[2] = new sphere(vec3(1, 0, -1), 0.5, new metal(vec3(0.8, 0.6, 0.2)));
    list[2] = new sphere(vec3(1, 0, -1), 0.5, new metal(vec3(0.8, 0.6, 0.2), 0.3));
    // list[3] = new sphere(vec3(-1, 0, -1), 0.5, new metal(vec3(0.8, 0.8, 0.8)));
    list[3] = new sphere(vec3(-1, 0, -1), 0.5, new metal(vec3(0.8, 0.8, 0.8), 1.0));
    hittable *world = new hittable_list(list, 4);

     camera cam;
    ofstream outFile("output_" + to_string(nx) + "x" + to_string(ny) + ".ppm");

      vector<thread> threads;


//     for (int i = 0; i < 5; ++i)
//     {

//                 threads.push_back(thread([i]() {

//             cout << "Hello from lamda thread " <<i<< this_thread::get_id() << endl;

//         }));

//            
//     }


//         for (auto &thread : threads)
//     {

//                 thread.join();

//            
//     }

    int num = 100;

    int numThread = 4;
    for (int k = 0; k < numThread; k++)
    {
        // threads.push_back(thread([=]() {

            camera cam;
//             cout << "Hello from lamda thread " <<i<< this_thread::get_id() << endl;
             for (int j = (ny - 1) - k; j >= 0; j -= numThread)
            {
                progressIdx = float(ny - 1 - j) / (ny - 1);
                for (int i = 0; i < nx; i++)
                {

                    vec3 col(0, 0, 0);
                    for (int s = 0; s < ns; s++)
                    {
                        float u = float(i + random_double()) / float(nx);
                        float v = float(j + random_double()) / float(ny);
                        ray r = cam.get_ray(u, v);
                        col += color(r, world, 0);
                        // col += color(r, world);
                    }
                    col /= float(ns);
                    col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));

                    int ir = int(255.99 * col[0]);
                    int ig = int(255.99 * col[1]);
                    int ib = int(255.99 * col[2]);

                    device_pixel(framebuffer, i, ny - 1 - j, ir, ig, ib);
                    // std::cout << ir << " " << ig << " " << ib << "\n";
                    // outFile << ir << " " << ig << " " << ib << endl;
                    // outFile << icolor.r << " " << icolor.g << " " << icolor.b << endl;
        totTime = (float)((int)time(nullptr) - (int)curTime);
                }
            }

//         }));
    }

            for (auto &thread : threads)
    {

                thread.join();

           
    }
        cout << "totTime = " << time(nullptr) - curTime << endl;


        return;
}

int main(int, char **)
{
    totTime = 0.0;
       curTime = time(nullptr);
    // framebuffer = new int[display_w * display_h];
    framebuffer = new int[1280 * 720]; //todo

    thread t(RayTracing);
    // RayTracing();

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

            ImGui::RadioButton("Wireframe", &gRenderType, Wireframe);
            ImGui::SameLine();
            ImGui::RadioButton("PureTriangle", &gRenderType, PureTriangle);
            ImGui::SameLine();
            ImGui::RadioButton("radio c", &gRenderType, 2);

            ImGui::SliderFloat("cam.pos.z", &cam.pos.z, -100.0f, 100.0f);             
            ImGui::SliderFloat("cam.pos.y", &cam.pos.y, -100.0f, 100.0f);             
            ImGui::SliderFloat("cam.pos.x", &cam.pos.x, -100.0f, 100.0f);             
            ImGui::SliderFloat("y_ang", &y_ang, 0.0f, 360.0f);             
            ImGui::ColorEdit3("clear color", (float *)&clear_color); // Edit 3 floats representing a color

            if (ImGui::Button("Button")) // Buttons return true when clicked (most widgets return true when edited/activated)
                counter++;
            ImGui::SameLine();
            ImGui::Text("counter = %d", counter);

                    // Animate a simple progress bar
            // if (true)
            // {
            //     progressIdx += progressDir * 0.4f * ImGui::GetIO().DeltaTime;
            //     if (progressIdx >= +1.1f)
            //     {
            //         progressIdx = +1.1f;
            //         progressDir *= -1.0f;
            //     }
            //     if (progressIdx <= -0.1f)
            //     {
            //         progressIdx = -0.1f;
            //         progressDir *= -1.0f;
            //     }
            // }

            // Typically we would use ImVec2(-1.0f,0.0f) to use all available width, or ImVec2(width,0.0f) for a specified width. ImVec2(0.0f,0.0f) uses ItemWidth.
            ImGui::ProgressBar(progressIdx, ImVec2(0.0f, 0.0f));
            ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
            ImGui::Text("Progress Bar");
            ImGui::Text("Total time %.1f s", totTime);
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
            InitDemo9_2();
            gInited = true;
                    // RayTracing();
            // thread t(RayTracing);
            // cout<< display_w <<display_h<<endl;

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
                // RGB2Color(framebuffer[i * display_w + j], greyLevel * clear_color.x, greyLevel * clear_color.y, greyLevel * clear_color.z);
            }
        }
        // DrawDemo9_2();
        // RayTracing();
        DrawFrame();


        // glUseProgram(0); // You may want this if using this code in an OpenGL 3+ context where shaders may be bound
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
        SDL_GL_SwapWindow(window);
    }

    // t.join();


    // Cleanup
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImGui::DestroyContext();

    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}

