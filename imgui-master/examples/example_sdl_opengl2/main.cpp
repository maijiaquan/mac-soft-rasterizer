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

//using namespace std;


//全局变量
int *framebuffer; //帧缓冲
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
void device_pixel(int x, int y, int c);
void device_pixel(int x, int y, int c)
{
    framebuffer[y*display_w + x] = c;
}
void device_draw_line(int x1, int y1, int x2, int y2, int c);
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
    for (int i = 0; i < display_h; i++)
    {
        for (int j = 0; j < display_w; j++)
        {
            RGB2Color(framebuffer[i * display_w + j], 255 * clear_color.x, 255 * clear_color.y, 255 * clear_color.z);
        }
    }


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
         std::cout << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[0].x << std::endl;
         std::cout << "x1 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[0].x << std::endl;
         std::cout << "y1 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[0].y << std::endl;
         std::cout << "x2 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[1].x << std::endl;
         std::cout << "y2 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[1].y << std::endl;
         std::cout << "x3 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[2].x << std::endl;
         std::cout << "y3 = " << rend_list_ptr->poly_ptrs[idx_poly]->tvlist[2].y << std::endl;
         std::cout << "-----" << std::endl;



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

void DrawFrame();
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

    cout<<testPrint()<<"hello"<<endl;
//    string s = "fff";
//    testPrint(s);

    // Setup SDL
    if (SDL_Init(SDL_INIT_VIDEO|SDL_INIT_TIMER) != 0)
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
    SDL_Window* window = SDL_CreateWindow("Dear ImGui SDL2+OpenGL example", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1280, 720, SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);
    SDL_GLContext gl_context = SDL_GL_CreateContext(window);
    SDL_GL_SetSwapInterval(1); // Enable vsync

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;

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

    bool show_demo_window = true;
    bool show_another_window = false;

//    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);


    // Main loop
    bool done = false;
    int tmpCount = 0;
    cout<<"(int)io.DisplaySize.x = "<<(int)io.DisplaySize.x<<"(int)io.DisplaySize.y = "<<(int)io.DisplaySize.y<<endl;
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
        glViewport(0, 0, (int)io.DisplaySize.x, (int)io.DisplaySize.y);
        if(tmpCount < 1)
        {
            display_w = (int)io.DisplaySize.x;
            display_h = (int)io.DisplaySize.y;
            framebuffer = new int[display_w*display_h];
            GameInit();
            tmpCount++;
            // cout<<"while (int)io.DisplaySize.x = "<<(int)io.DisplaySize.x<<"while (int)io.DisplaySize.y = "<<(int)io.DisplaySize.y<<endl;
            cout<<"h = "<<display_h<<"w = "<<display_w<<endl;
            // cout<<"h = "<<(int)io.DisplaySize.x<<"w = "<<(int)io.DisplaySize.y<<endl;
        }
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);

        GameMain();
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
