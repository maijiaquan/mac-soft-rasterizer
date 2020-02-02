#pragma once

#include <math.h>
#include <string.h>
#include "struct.h"
#include "parser.h"

#include <iostream>
using namespace std;

// storage for our lookup tables
extern float cos_look[361]; // 1 extra so we can store 0-360 inclusive
extern float sin_look[361]; // 1 extra so we can store 0-360 inclusive

extern MATV1 materials[MAX_MATERIALS]; // materials in system
extern int num_materials;			   // current number of materials

extern LIGHTV1 lights[MAX_LIGHTS]; // lights in system
extern int num_lights;			   // current number of light

extern char texture_path[80]; // root path to ALL textures, make current directory for now
// extern BITMAP_FILE bitmap16bit; // a 16 bit bitmap file


void Build_Sin_Cos_Tables(void);

//颜色转换
void Color2RGB(const int &hex, int &r, int &g, int &b);
void RGB2Color(int &hex, const int &r, const int &g, const int &b);


//内联函数
inline void VECTOR3D_INITXYZ(VECTOR3D_PTR v, float x, float y, float z)
{
    (v)->x = (x);
    (v)->y = (y);
    (v)->z = (z);
}
inline void VECTOR4D_COPY(VECTOR4D_PTR vdst, VECTOR4D_PTR vsrc)
{
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
    (vdst)->w = (vsrc)->w;
}

inline void VECTOR4D_INITXYZ(VECTOR4D_PTR v, float x, float y, float z)
{
    (v)->x = (x);
    (v)->y = (y);
    (v)->z = (z);
    (v)->w = 1.0;
}

inline void POINT3D_COPY(POINT3D_PTR vdst, POINT3D_PTR vsrc)
{
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
}

inline void VECTOR4D_ZERO(VECTOR4D_PTR v)
{
    (v)->x = (v)->y = (v)->z = 0.0;
    (v)->w = 1.0;
}

inline void VECTOR3D_COPY(VECTOR3D_PTR vdst, VECTOR3D_PTR vsrc)
{
    (vdst)->x = (vsrc)->x;
    (vdst)->y = (vsrc)->y;
    (vdst)->z = (vsrc)->z;
}

inline void VECTOR3D_ZERO(VECTOR3D_PTR v)
{
    (v)->x = (v)->y = (v)->z = 0.0;
}

USHORT RGB16Bit565(int r, int g, int b);

void Reset_RENDERLIST4DV1(RENDERLIST4DV1_PTR renderList);
void Init_CAM4DV1(CAM4DV1_PTR cam, int attr, POINT4D_PTR cam_pos,
                  VECTOR4D_PTR cam_dir, VECTOR4D_PTR cam_target,
                  float near_clip_z, float far_clip_z, float fov,
                  float viewport_width, float viewport_height);

void MAT_IDENTITY_4X4(MATRIX4X4_PTR m);
void PLANE3D_Init(PLANE3D_PTR plane, POINT3D_PTR p0,
                  VECTOR3D_PTR normal, int normalize);

void VECTOR3D_Normalize(VECTOR3D_PTR va);
void VECTOR3D_Normalize(VECTOR3D_PTR va, VECTOR3D_PTR vn);
float VECTOR3D_Length(VECTOR3D_PTR va);

int Insert_POLYF4DV1_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, POLYF4DV1_PTR poly);

void Build_XYZ_Rotation_MATRIX4X4(float theta_x, // euler angles
                                  float theta_y,
                                  float theta_z,
                                  MATRIX4X4_PTR mrot);

float Fast_Sin(float theta);
float Fast_Cos(float theta);
void Mat_Mul_4X4(MATRIX4X4_PTR ma, MATRIX4X4_PTR mb, MATRIX4X4_PTR mprod);

void Mat_Init_4X4(MATRIX4X4_PTR ma,
                  float m00, float m01, float m02, float m03,
                  float m10, float m11, float m12, float m13,
                  float m20, float m21, float m22, float m23,
                  float m30, float m31, float m32, float m33);



void Transform_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, // render list to transform
                              MATRIX4X4_PTR mt,             // transformation matrix
                              int coord_select);

void Model_To_World_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, POINT4D_PTR world_pos,
                                   int coord_select = TRANSFORM_LOCAL_TO_TRANS);

void Build_CAM4DV1_Matrix_Euler(CAM4DV1_PTR cam, int cam_rot_seq);

void World_To_Camera_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, CAM4DV1_PTR cam);

void Camera_To_Perspective_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, CAM4DV1_PTR cam);

void Perspective_To_Screen_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, CAM4DV1_PTR cam);


void Mat_Mul_VECTOR4D_4X4(VECTOR4D_PTR  va, MATRIX4X4_PTR mb, VECTOR4D_PTR  vprod);

void VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vsum);
VECTOR4D VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb);


//new


void Reset_OBJECT4DV1(OBJECT4DV1_PTR obj);


void Transform_OBJECT4DV1(OBJECT4DV1_PTR obj, MATRIX4X4_PTR mt, int coord_select, int transform_basis);
void Model_To_World_OBJECT4DV1(OBJECT4DV1_PTR obj, int coord_select = TRANSFORM_LOCAL_TO_TRANS);
void World_To_Camera_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam);
void Camera_To_Perspective_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam);
void Perspective_To_Screen_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam);


void Remove_Backfaces_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam);  //消除背面
void Remove_Backfaces_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, CAM4DV1_PTR cam);

int Cull_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam, int cull_flags);

int Insert_POLY4DV1_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, POLY4DV1_PTR poly);

int Insert_OBJECT4DV1_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, OBJECT4DV1_PTR obj, int insert_local = 0);

void VECTOR4D_Build(VECTOR4D_PTR init, VECTOR4D_PTR term, VECTOR4D_PTR result);
void VECTOR4D_Cross(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vn);
VECTOR4D VECTOR4D_Cross(VECTOR4D_PTR va, VECTOR4D_PTR vb);


float VECTOR4D_Dot(VECTOR4D_PTR va, VECTOR4D_PTR vb);





                        // lighting system
int Init_Light_LIGHTV1(int index,		   // index of light to create (0..MAX_LIGHTS-1)
					   int _state,		   // state of light
					   int _attr,		   // type of light, and extra qualifiers
					   RGBAV1 _c_ambient,  // ambient light intensity
					   RGBAV1 _c_diffuse,  // diffuse light intensity
					   RGBAV1 _c_specular, // specular light intensity
					   POINT4D_PTR _pos,   // position of light
					   VECTOR4D_PTR _dir,  // direction of light
					   float _kc,		   // attenuation factors
					   float _kl,
					   float _kq,
					   float _spot_inner, // inner angle for spot light
					   float _spot_outer, // outer angle for spot light
					   float _pf);		  // power factor/falloff for spot lights

int Reset_Lights_LIGHTV1(void);
void Reset_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list);
void Reset_OBJECT4DV2(OBJECT4DV2_PTR obj);


void Transform_OBJECT4DV2(OBJECT4DV2_PTR obj, MATRIX4X4_PTR mt,
						  int coord_select, int transform_basis, int all_frames = 0);
void Model_To_World_OBJECT4DV2(OBJECT4DV2_PTR obj, int coord_select = TRANSFORM_LOCAL_TO_TRANS, int all_frames = 0);
int Insert_OBJECT4DV2_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
									 OBJECT4DV2_PTR obj,
									 int insert_local);


void Remove_Backfaces_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list, CAM4DV1_PTR cam);




int Light_RENDERLIST4DV2_World16(RENDERLIST4DV2_PTR rend_list, // list to process
								 CAM4DV1_PTR cam,			   // camera position
								 LIGHTV1_PTR lights,		   // light list (might have more than one)
								 int max_lights);			   // maximum lights in list

void World_To_Camera_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
									CAM4DV1_PTR cam);





void Sort_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list, int sort_method);

void Camera_To_Perspective_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
										  CAM4DV1_PTR cam);
void Perspective_To_Screen_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
										  CAM4DV1_PTR cam);

int Init_OBJECT4DV2(OBJECT4DV2_PTR obj, int _num_vertices, int _num_polys, int _num_frames, int destroy = 0);

float Compute_OBJECT4DV2_Radius(OBJECT4DV2_PTR obj);

int Compute_OBJECT4DV2_Poly_Normals(OBJECT4DV2_PTR obj);
int Compute_OBJECT4DV2_Vertex_Normals(OBJECT4DV2_PTR obj);

void VECTOR4D_Normalize(VECTOR4D_PTR va);
float VECTOR4D_Length_Fast(VECTOR4D_PTR va);
float Fast_Distance_3D(float x, float y, float z);

int Insert_POLY4DV2_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list, POLY4DV2_PTR poly);

inline float VECTOR4D_Length_Fast2(VECTOR4D_PTR va)
{
	// this function computes the distance from the origin to x,y,z

	int temp;	// used for swaping
	int x, y, z; // used for algorithm

	// make sure values are all positive
	x = fabs(va->x) * 1024;
	y = fabs(va->y) * 1024;
	z = fabs(va->z) * 1024;

	// sort values
	if (y < x)
		SWAP(x, y, temp)
	if (z < y)
		SWAP(y, z, temp)
	if (y < x)
		SWAP(x, y, temp)

	int dist = (z + 11 * (y >> 5) + (x >> 2));

	// compute distance with 8% error
	return ((float)(dist >> 10));

} // end VECTOR4D_Length_Fast

// avg z-compare
int Compare_AvgZ_POLYF4DV2(const void *arg1, const void *arg2);

// near z-compare
int Compare_NearZ_POLYF4DV2(const void *arg1, const void *arg2);

// far z-compare
int Compare_FarZ_POLYF4DV2(const void *arg1, const void *arg2);


int Destroy_OBJECT4DV2(OBJECT4DV2_PTR obj);

float VECTOR4D_Length(VECTOR4D_PTR va);

inline void VERTEX4DTV1_COPY(VERTEX4DTV1_PTR vdst, VERTEX4DTV1_PTR vsrc)
{
	*vdst = *vsrc;
}


