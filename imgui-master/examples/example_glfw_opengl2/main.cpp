// dear imgui: standalone example application for GLFW + OpenGL2, using legacy fixed pipeline
// If you are new to dear imgui, see examples/README.txt and documentation at the top of imgui.cpp.
// (GLFW is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan graphics context creation, etc.)

// **DO NOT USE THIS CODE IF YOUR CODE/ENGINE IS USING MODERN OPENGL (SHADERS, VBO, VAO, etc.)**
// **Prefer using the code in the example_glfw_opengl2/ folder**
// See imgui_impl_glfw.cpp for details.

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include <stdio.h>
#include <iostream>
#include <sys/time.h>                // for gettimeofday()
//#include <opencv2/core/core.hpp>
//#include <opencv2/highgui/highgui.hpp>
#include <ctime>
#include "ds.h"
//using namespace cv;

using namespace std;

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GLFW/glfw3.h>

// [Win32] Our example includes a copy of glfw3.lib pre-compiled with VS2010 to maximize ease of testing and compatibility with old VS compilers.
// To link with VS2010-era libraries, VS2015+ requires linking with legacy_stdio_definitions.lib, which we do using this pragma.
// Your own project should not be affected, as you are likely to link with a newer binary of GLFW that is adequate for your version of Visual Studio.
#if defined(_MSC_VER) && (_MSC_VER >= 1900) && !defined(IMGUI_DISABLE_WIN32_FUNCTIONS)
#pragma comment(lib, "legacy_stdio_definitions")
#endif

static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

//全局变量
int *framebuffer;
ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

int display_w, display_h;
float colorR = 0;

POINT4D cam_pos = {0, 0, -100, 1};
VECTOR4D cam_dir = {0, 0, 0, 1};

// all your initialization code goes here...
VECTOR4D vscale = {.5, .5, .5, 1},
         vpos = {0, 0, 0, 1},
         vrot = {0, 0, 0, 1};

RENDERLIST4DV1 rend_list;            // the single renderlist
POLYF4DV1 poly1;                    // our lonely polygon
CAM4DV1 cam;                        // the single camera
POINT4D poly1_pos = {0, 0, 100, 1}; // world position of polygon
USHORT(*RGB16Bit)
(int r, int g, int b) = nullptr;
void DrawFrame();
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
void device_pixel(int x, int y, int c)
{
  framebuffer[y*display_w + x] = c;
}
void device_draw_line(int x1, int y1, int x2, int y2, int c)
{
    // std::cout<<"x1 = "<<x1<<"y1 = "<<y1<<std::endl;
    // std::cout<<"x2 = "<<x2<<"y2 = "<<y2<<std::endl;

    int x, y, rem = 0;
    if (x1 == x2 && y1 == y2)
    {
      device_pixel(x1, y1, c);
    }
    else if (x1 == x2)
    {
        int inc = (y1 <= y2) ? 1 : -1;
        for (y = y1; y != y2; y += inc)
          device_pixel(x1, y, c);
        device_pixel(x2, y2, c);
    }
    else if (y1 == y2)
    {
        int inc = (x1 <= x2) ? 1 : -1;
        for (x = x1; x != x2; x += inc)
            device_pixel(x, y1, c);
        device_pixel(x2, y2, c);
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
                device_pixel(x, y, c);
                rem += dy;
                if (rem >= dx)
                {
                    rem -= dx;
                    y += (y2 >= y1) ? 1 : -1;
                    device_pixel(x, y, c);
                }
            }
            device_pixel(x2, y2, c);
        }
        else
        {
            if (y2 < y1)
                x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
            for (x = x1, y = y1; y <= y2; y++)
            {
                device_pixel(x, y, c);
                rem += dx;
                if (rem >= dy)
                {
                    rem -= dy;
                    x += (x2 >= x1) ? 1 : -1;
                    device_pixel(x, y, c);
                }
            }
            device_pixel(x2, y2, c);
        }
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
    glOrtho(0.0,display_w,0.0,display_h,0.0,1.0);
    //    glOrtho(0.0,display_w,0.0,display_h,0.0,1.0);
    
    
    glBegin(GL_POINTS);
    float delta = (float)1.0/255;
    
    for(int i = 0; i < display_w; i++)
    {
        glColor3f(0,0,0);
        for(int j = 0; j < display_h; j++)
        {
            glVertex3f(i,j,0);
        }
    }
    glEnd();
}
void Color2RGB(const int &hex, int &r, int &g, int &b);

void RGB2Color(int &hex, const int &r, const int &g, const int &b);

void GameMain();
void GameInit();

void Color2RGB(const int &c, int &r, int &g, int &b)
{
  r = (0xff << 16 & c) >> 16;
  g = (0xff << 8 & c) >> 8;
  b = 0xff & c;
}

void RGB2Color(int &c, const int &r, const int &g, const int &b)
{
  c = (r<<16) | (g<<8) | b;
}


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
    Init_CAM4DV1(&cam,              // the camera object
                 CAM_MODEL_EULER, // euler camera model
                 &cam_pos,          // initial camera position
                 &cam_dir,          // initial camera angles
                 NULL,              // no initial target
                 50.0,              // near and far clipping planes
                 500.0,
                 90.0,           // field of view in degrees
                 WINDOW_WIDTH, // size of final screen viewport
                 WINDOW_HEIGHT);
}
void GameMain()
{
  // for(int i = 0; i < display_h; i++)
  // {
    // for(int j = 0; j < display_w; j++)
    // {
        // RGB2Color(framebuffer[i*display_w + j], 255*clear_color.x, 255*clear_color.y, 255*clear_color.z);
    // }
  // }


  for(int i = 0; i < 100; i++)
  {
    for(int j = 0; j < 200; j++)
    {
        RGB2Color(framebuffer[i*display_w + j], 255*clear_color.x, 255*clear_color.y, 255*clear_color.z);
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
        RGB2Color(c,255,255,255);

       device_draw_line(x1, y1, x2, y2,c); //3 1
       device_draw_line(x1, y1, x3, y3,c); //3 1
       device_draw_line(x2, y2, x3, y3,c); //3 1

    } // end for poly
}
void DrawFrame()
{

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  //        glOrtho(0.0, 500.0, 0.0, 500.0, -1, 1); //设置正射投影的剪裁空间
    //        glOrtho2D(0.0, 500.0, 50 0.0, 0.0);
    
    //        glOrtho2D(0.0, 500.0, 500.0, 0.0);
    //        glOrtho2D(0.0, 500.0, 0.0, 500.0);
    glOrtho(0.0,display_w,0.0,display_h,0.0,1.0);
    //    glOrtho(0.0,display_w,0.0,display_h,0.0,1.0);
    
    
    glBegin(GL_POINTS);
    float delta = (float)1.0/255;
    
    for(int i = 0; i < display_h; i++)
    {
      // glColor3f(colorR,0,0);
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

        int r,g,b;
        // Color2RGB(c, r, g, b);
        for(int j = 0; j < display_w; j++)
        {
          Color2RGB(framebuffer[i*display_w + j], r, g, b);
          // cout<<"r = "<<r<<"g = "<<g<<"b = "<<b;
          // glColor3f((float)r/255,(float)g/255,(float)b/255);
          glColor3f((float)r/255,(float)g/255,(float)b/255);
          // glVertex3f(i,j,0);
          glVertex3f(j,display_h-i,0);
        }
    }
    
    glEnd();
}

int main(int, char**)
{
//    #ifdef __APPLE__
//        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
//    #endif
//    Mat image;
//    image = imread("0.jpg", CV_LOAD_IMAGE_COLOR);
//    if(! image.data )                              // Check for invalid input
//    {
//        cout <<  "Could not open or find the image" << std::endl ;
//        return -1;
//    }
//
//    namedWindow( "Display window", CV_WINDOW_AUTOSIZE );// Create a window for display.
//    imshow( "Display window", image );
    
    
    // Setup window
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
        return 1;
    GLFWwindow* window = glfwCreateWindow(600, 600, "Dear ImGui GLFW+OpenGL2 example", NULL, NULL);
    if (window == NULL)
        return 1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync
    

    
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
    
    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsClassic();
    
    // Setup Platform/Renderer bindings
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();
    
    // Load Fonts
    // - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
    // - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
    // - If the file cannot be loaded, the function will return NULL. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
    // - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which ImGui_ImplXXXX_NewFrame below will call.
    // - Read 'docs/FONTS.txt' for more instructions and details.
    // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
    //io.Fonts->AddFontDefault();
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/ProggyTiny.ttf", 10.0f);
    //ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, NULL, io.Fonts->GetGlyphRangesJapanese());
    //IM_ASSERT(font != NULL);
    
    // Our state
    bool show_demo_window = true;
    bool show_another_window = false;
    // ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
    GameInit();
    
    // Main loop
    while (!glfwWindowShouldClose(window))
    {
        // Poll and handle events (inputs, window resize, etc.)
        // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
        // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application.
        // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application.
        // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
        glfwPollEvents();
        
        // Start the Dear ImGui frame
        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        
        // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
        if (show_demo_window)
            ImGui::ShowDemoWindow(&show_demo_window);
        
        // 2. Show a simple window that we create ourselves. We use a Begin/End pair to created a named window.
        {
            static float f = 0.0f;
            static int counter = 0;
            
            ImGui::Begin("Hello, world!");                          // Create a window called "Hello, world!" and append into it.
            
            ImGui::Text("This is some useful text.");               // Display some text (you can use a format strings too)
            ImGui::Checkbox("Demo Window", &show_demo_window);      // Edit bools storing our window open/close state
            ImGui::Checkbox("Another Window", &show_another_window);
            
            ImGui::SliderFloat("float", &f, 0.0f, 1.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
            ImGui::ColorEdit3("clear color", (float*)&clear_color); // Edit 3 floats representing a color
            
            if (ImGui::Button("Button"))                            // Buttons return true when clicked (most widgets return true when edited/activated)
                counter++;
            ImGui::SameLine();
            ImGui::Text("counter = %d", counter);
            
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            ImGui::End();
        }
        
        // 3. Show another simple window.
        if (show_another_window)
        {
            ImGui::Begin("Another Window", &show_another_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
            ImGui::Text("Hello from another window!");
            if (ImGui::Button("Close Me"))
                show_another_window = false;
            ImGui::End();
        }
        
        // Rendering
        
        ImGui::Render();
        
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        framebuffer = new int[display_w*display_h];
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);
        
        timeval t1, t2;
        double elapsedTime;
        gettimeofday(&t1, NULL);
        

        GameMain();
        DrawFrame();
        
        gettimeofday(&t2, NULL);
        
        // compute and print the elapsed time in millisec
        elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
        elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
//        cout << elapsedTime << " ms.\n";
//        cout << 1000/elapsedTime << " fps.\n";
        
        
        // If you are using this code with non-legacy OpenGL header/contexts (which you should not, prefer using imgui_impl_opengl3.cpp!!),
        // you may need to backup/reset/restore current shader using the commented lines below.
        //GLint last_program;
        //glGetIntegerv(GL_CURRENT_PROGRAM, &last_program);
        //glUseProgram(0);
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
        //glUseProgram(last_program);
        
        
        glfwMakeContextCurrent(window);
        
        
        
        glfwSwapBuffers(window);
    }
    
    // Cleanup
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    
    glfwDestroyWindow(window);
    glfwTerminate();
    
    return 0;
}
