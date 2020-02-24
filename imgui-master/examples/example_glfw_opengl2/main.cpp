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

#include <thread>
#include <fstream>
#include <vector>
#include <ctime>

static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}


//全局变量
int *framebuffer;
ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
int display_w, display_h;
float colorR = 0;

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

void RGB2Color(int &hex, const int &r, const int &g, const int &b);
void RGB2Color(int &c, const int &r, const int &g, const int &b)
{
  c = (r<<16) | (g<<8) | b;
}

void device_pixel(int x, int y, int c)
{
  framebuffer[y*display_w + x] = c;
}

void device_pixel( const int &x, const int &y, const int &r, const int &g, const int &b)
{
    int c;
    RGB2Color(c, r, g, b);
    device_pixel( x, y, c);
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


void GameMain();
void GameInit();

void Color2RGB(const int &c, int &r, int &g, int &b)
{
  r = (0xff << 16 & c) >> 16;
  g = (0xff << 8 & c) >> 8;
  b = 0xff & c;
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


//   for(int i = 0; i < 100; i++)
  for(int i = 0; i < display_h; i++)
  {
    // for(int j = 0; j < 200; j++)
    for(int j = 0; j < display_w; j++)
    {
        RGB2Color(framebuffer[i*display_w + j], 255*clear_color.x, 255*clear_color.y, 255*clear_color.z);
    }
  }
    
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

float progressDone = 0.0f, progressDir = 1.0f;

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

class material  {
    public:
        //r_in为入射光线, scattered为散射光线, attenuation 为衰减量
        virtual bool scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation, ray &scattered) const = 0;
};


struct hit_record
{
    float t;   //命中射线的长度
    vec3 p;    //命中终点坐标
    vec3 normal; //命中点的法向量
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
        scattered = ray(rec.p, s_world - rec.p); //scattered为散射光线
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
        vec3 v = unit_vector(r_in.direction());
        vec3 n = rec.normal;
        vec3 p = rec.p;
        vec3 r = reflect(v, n);
        vec3 offset = fuzz * random_in_unit_sphere();
        scattered = ray(p, r+offset);

        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > 0);
    }
};

// class dielectric : public material {
//     public:
//         dielectric(float ri) : ref_idx(ri) {}
//         virtual bool scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation, ray &scattered) const
//         {
//             vec3 outward_normal;
//             vec3 reflected = reflect(r_in.direction(), rec.normal);
//             float ni_over_nt;
//             attenuation = vec3(1.0, 1.0, 1.0); //玻璃全透明不吸收光
//             vec3 refracted;

//             if (dot(r_in.direction(), rec.normal) > 0) {
//                 outward_normal = -rec.normal;
//                 ni_over_nt = ref_idx;
//             }
//             else {
//                 outward_normal = rec.normal;
//                 ni_over_nt = 1.0 / ref_idx;
//             }

//             if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted)) {
//                 scattered = ray(rec.p, refracted);
//             }
//             else {
//                 scattered = ray(rec.p, reflected);
//                 return false;
//             }

//             return true;
//         }

//         float ref_idx;
// };


class dielectric : public material {
    public:
        dielectric(float ri) : ref_idx(ri) {}
        virtual bool scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation, ray &scattered) const
        {
            vec3 outward_normal;
            vec3 reflected = reflect(r_in.direction(), rec.normal);
            float ni_over_nt;
            attenuation = vec3(1.0, 1.0, 1.0);
            vec3 refracted;

            float reflect_prob;
            float cosine;

            if (dot(r_in.direction(), rec.normal) > 0) {
                 outward_normal = -rec.normal;
                 ni_over_nt = ref_idx;
                 cosine = ref_idx * dot(r_in.direction(), rec.normal)
                        / r_in.direction().length();
            }
            else {
                 outward_normal = rec.normal;
                 ni_over_nt = 1.0 / ref_idx;
                 cosine = -dot(r_in.direction(), rec.normal)
                        / r_in.direction().length();
            }

            if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted)) {
               reflect_prob = schlick(cosine, ref_idx);
            }
            else {
               reflect_prob = 1.0;
            }

            if (random_double() < reflect_prob) {
               scattered = ray(rec.p, reflected);
            }
            else {
               scattered = ray(rec.p, refracted);
            }

            return true;
        }

        float ref_idx;
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

    //如果命中了，命中记录保存到rec
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
    //如果命中了，命中记录保存到rec
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


//在场景中发射一条射线，并采样该射线最终输出到屏幕的颜色值
vec3 color(const ray &r, hittable *world, int depth)
{
    hit_record rec;
    if (world->hit(r, 0.001, MAXFLOAT, rec)) //射线命中物体
    {
        ray scattered;  //散射光线
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


void Framebuffer2File(int nx, int ny, int ns, int *fb, ofstream &outFile)
{
    for (int j = (ny - 1) ; j >= 0; j --)
    {
        progressDone = float(ny - 1 - j) / (ny - 1);
        for (int i = 0; i < nx; i++)
        {
            // device_pixel(framebuffer, i, ny - 1 - j, ir, ig, ib);
            // device_pixel(fb, x, y, c);
            int x = i;
            int y = ny - 1 - j;
            int c = fb[y * display_w + x];
            int r,g,b = 0;

            Color2RGB(c, r, g, b);

            // std::cout << r << " " << g << " " << b << "\n";
            outFile << r << " " << g << " " << b << endl;
            // outFile << icolor.r << " " << icolor.g << " " << icolor.b << endl;
            // totTime = (float)((int)time(nullptr) - (int)curTime);
        }
    }
}


using namespace std;
time_t curTime;
float totTime;
float timeRemaining;

int nx, ny, ns;
/*

done
输出到ppm文件
进度条除颤
剩余时间（不精确）

todo
剩余时间
精确到微秒
换种方式多线程

*/
void RayTracing()
{

    ns = 100;

    // nx = 1200;
    // ny = 600;
    // nx = 800;
    // ny = 400;
    nx = 400;
    ny = 200;
    // nx = 200;
    // ny = 100;
    usleep(1000);    // will sleep for 1 ms
    usleep(1);       // will sleep for 0.001 ms
    usleep(1000000); // will sleep for 1 s
    usleep(1000000); // will sleep for 1 s
                     // usleep(1000000); // will sleep for 1 s
                     // int nx = 600;
                     // int ny = 300;

    std::cout << "P3\n"
              << nx << " " << ny << "\n255\n";

    hittable *list[4];

    list[0] = new sphere(vec3(0, 0, -1), 0.5, new lambertian(vec3(0.1, 0.2, 0.5)));
    list[1] = new sphere(vec3(0, -100.5, -1), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
    list[2] = new sphere(vec3(1, 0, -1), 0.5, new metal(vec3(0.8, 0.6, 0.2), 0.0));
    list[3] = new sphere(vec3(-1, 0, -1), 0.5, new dielectric(1.5));

    hittable *world = new hittable_list(list, 4);

    camera cam;
    ofstream outFile("output_" + to_string(nx) + "x" + to_string(ny) + ".ppm");
    outFile << "P3\n"
            << nx << " " << ny << "\n255\n";
    vector<thread> threads;


    int num = 100;

    int numThread = 20;
    int numPixelRendered = 0;
    int doneRecord = 0;
    int lastTime = -1;
    for (int k = 0; k < numThread; k++)
    {
        threads.push_back(thread([=, &numPixelRendered, &lastTime, &doneRecord]() {

            camera cam;
            int numPixelTotal = nx * ny;
            for (int j = (ny - 1) - k; j >= 0; j -= numThread)
            {
                progressDone = float(numPixelRendered) / (numPixelTotal);
                // cout<<"progress Done = "<<progressDone<<endl;
                totTime = (float)((int)time(nullptr) - (int)curTime);
                if(lastTime < 0)
                {
                    lastTime = totTime;
                    doneRecord = numPixelRendered;
                }
                else
                {
                    if(totTime - lastTime >= 1.0)
                    {
                        int pixelPerSecond = numPixelRendered - doneRecord;
                        cout<<"pixelPerSecond = "<<pixelPerSecond<<endl;

                        timeRemaining = (float) (numPixelTotal - numPixelRendered) / pixelPerSecond; //法一：根据当前的瞬时速度来更新显示
                        lastTime = totTime;
                        doneRecord = numPixelRendered;
                    }
                }

                // timeRemaining = (float)totTime * (1 - progressDone) / progressDone; //法二：根据过去的平均时间来显示 x/t=pTodo/pDone -> x= t*pTodo/pDone = t*(1-pDone) /pDone

                // cout<<"timeRemaining = "<<timeRemaining<<endl;
                for (int i = 0; i < nx; i++)
                {
                    // usleep(1000); // will sleep for 1 ms
                    numPixelRendered++;

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

                    device_pixel( i, ny - 1 - j, ir, ig, ib);
                    // std::cout << ir << " " << ig << " " << ib << "\n";
                    // outFile << ir << " " << ig << " " << ib << endl;
                    // outFile << icolor.r << " " << icolor.g << " " << icolor.b << endl;
                }
            }

        }));
    }

    for (auto &thread : threads)
    {
        thread.join();
    }
    cout << "totTime = " << time(nullptr) - curTime << endl;
    Framebuffer2File(nx, ny, ns, framebuffer, outFile);

    return;
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
