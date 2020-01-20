#pragma once

#include <math.h>
#include <string.h>
typedef unsigned short USHORT;
typedef unsigned char UCHAR;

#define RENDERLIST4DV1_MAX_POLYS 32768 // 16384
#define POLY4DV1_STATE_ACTIVE 0x0001
#define CAM_MODEL_EULER 0x0008
#define CAM_MODEL_UVN 0x0010
#define PI ((float)3.141592654f)
#define _RGB16BIT565(r, g, b) ((b & 31) + ((g & 63) << 5) + ((r & 31) << 11))
#define DEG_TO_RAD(ang) ((ang)*PI / 180.0)

// defines for small numbers
#define EPSILON_E3 (float)(1E-3)
#define EPSILON_E4 (float)(1E-4)
#define EPSILON_E5 (float)(1E-5)
#define EPSILON_E6 (float)(1E-6)

#define WINDOW_WIDTH 400 // size of window
#define WINDOW_HEIGHT 400

// transformation control flags
#define TRANSFORM_LOCAL_ONLY 0 // perform the transformation in place on the \
// local/world vertex list
#define TRANSFORM_TRANS_ONLY 1 // perfrom the transformation in place on the \
// "transformed" vertex list

#define TRANSFORM_LOCAL_TO_TRANS 2

// defines for camera rotation sequences
#define CAM_ROT_SEQ_XYZ 0
#define CAM_ROT_SEQ_YXZ 1
#define CAM_ROT_SEQ_XZY 2
#define CAM_ROT_SEQ_YZX 3
#define CAM_ROT_SEQ_ZYX 4
#define CAM_ROT_SEQ_ZXY 5


// states of polygons and faces
#define POLY4DV1_STATE_ACTIVE 0x0001
#define POLY4DV1_STATE_CLIPPED 0x0002
#define POLY4DV1_STATE_BACKFACE 0x0004


// used to compute the min and max of two expresions
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define SWAP(a, b, t) \
	{                 \
		t = a;        \
		a = b;        \
		b = t;        \
	}

// transformation control flags
#define TRANSFORM_LOCAL_ONLY 0
#define TRANSFORM_TRANS_ONLY 1
#define TRANSFORM_LOCAL_TO_TRANS 2

#define OBJECT4DV1_MAX_VERTICES 1024   // 64
#define OBJECT4DV1_MAX_POLYS 1024	  // 128
#define RENDERLIST4DV1_MAX_POLYS 32768 // 16384
#define PI ((float)3.141592654f)
#define DEG_TO_RAD(ang) ((ang)*PI / 180.0)

// states of polygons and faces
#define POLY4DV1_STATE_ACTIVE 0x0001
#define POLY4DV1_STATE_CLIPPED 0x0002
#define POLY4DV1_STATE_BACKFACE 0x0004

// states for objects
#define OBJECT4DV1_STATE_ACTIVE 0x0001
#define OBJECT4DV1_STATE_VISIBLE 0x0002
#define OBJECT4DV1_STATE_CULLED 0x0004

// shading mode of polygon
#define PLX_SHADE_MODE_PURE_FLAG 0x0000		 // this poly is a constant color
#define PLX_SHADE_MODE_CONSTANT_FLAG 0x0000  // alias
#define PLX_SHADE_MODE_FLAT_FLAG 0x2000		 // this poly uses flat shading
#define PLX_SHADE_MODE_GOURAUD_FLAG 0x4000   // this poly used gouraud shading
#define PLX_SHADE_MODE_PHONG_FLAG 0x6000	 // this poly uses phong shading
#define PLX_SHADE_MODE_FASTPHONG_FLAG 0x6000 // this poly uses phong shading (alias)

// attributes of polygons and polygon faces
#define POLY4DV1_ATTR_2SIDED 0x0001
#define POLY4DV1_ATTR_TRANSPARENT 0x0002
#define POLY4DV1_ATTR_8BITCOLOR 0x0004
#define POLY4DV1_ATTR_RGB16 0x0008
#define POLY4DV1_ATTR_RGB24 0x0010

#define POLY4DV1_ATTR_SHADE_MODE_PURE 0x0020
#define POLY4DV1_ATTR_SHADE_MODE_CONSTANT 0x0020 // (alias)
#define POLY4DV1_ATTR_SHADE_MODE_FLAT 0x0040
#define POLY4DV1_ATTR_SHADE_MODE_GOURAUD 0x0080
#define POLY4DV1_ATTR_SHADE_MODE_PHONG 0x0100
#define POLY4DV1_ATTR_SHADE_MODE_FASTPHONG 0x0100 // (alias)
#define POLY4DV1_ATTR_SHADE_MODE_TEXTURE 0x0200

#define PLX_2SIDED_FLAG 0x1000 // this poly is double sided
#define PLX_1SIDED_FLAG 0x0000 // this poly is single sided

#define PLX_COLOR_MODE_RGB_FLAG 0x8000	 // this poly uses RGB color
#define PLX_COLOR_MODE_INDEXED_FLAG 0x0000 // this poly uses an indexed 8-bit color

#define PLX_RGB_MASK 0x8000		   // mask to extract RGB or indexed color
#define PLX_SHADE_MODE_MASK 0x6000 // mask to extract shading mode
#define PLX_2SIDED_MASK 0x1000	 // mask for double sided
#define PLX_COLOR_MASK 0x0fff	  // xxxxrrrrggggbbbb, 4-bits per channel RGB \
                                   // xxxxxxxxiiiiiiii, indexed mode 8-bit index

#define RESET_BIT(word, bit_flag) ((word) = ((word) & (~bit_flag)))

// 相机旋转顺序
#define CAM_ROT_SEQ_XYZ 0
#define CAM_ROT_SEQ_YXZ 1
#define CAM_ROT_SEQ_XZY 2
#define CAM_ROT_SEQ_YZX 3
#define CAM_ROT_SEQ_ZYX 4
#define CAM_ROT_SEQ_ZXY 5


typedef unsigned short USHORT;

#define POLY4DV1_STATE_ACTIVE 0x0001
#define CAM_MODEL_EULER 0x0008
#define CAM_MODEL_UVN 0x0010

// general culling flags
#define CULL_OBJECT_X_PLANE 0x0001 // cull on the x clipping planes
#define CULL_OBJECT_Y_PLANE 0x0002 // cull on the y clipping planes
#define CULL_OBJECT_Z_PLANE 0x0004 // cull on the z clipping planes
#define CULL_OBJECT_XYZ_PLANES (CULL_OBJECT_X_PLANE | CULL_OBJECT_Y_PLANE | CULL_OBJECT_Z_PLANE)

#define RAND_RANGE(x, y) ((x) + (rand() % ((y) - (x) + 1)))


#define SET_BIT(word, bit_flag) ((word) = ((word) | (bit_flag)))

//光照宏定义
// create some constants for ease of access
#define AMBIENT_LIGHT_INDEX 0  // ambient light index
#define INFINITE_LIGHT_INDEX 1 // infinite light index
#define POINT_LIGHT_INDEX 2	// point light index
#define SPOT_LIGHT_INDEX 3	 // spot light index
#define SPOT_LIGHT1_INDEX 4	// point light index
#define SPOT_LIGHT2_INDEX 3	// spot light index

#define LIGHTV1_STATE_ON 1  // light on
#define LIGHTV1_STATE_OFF 0 // light off

// defines for light types
#define LIGHTV1_ATTR_AMBIENT 0x0001		// basic ambient light
#define LIGHTV1_ATTR_INFINITE 0x0002	// infinite light source
#define LIGHTV1_ATTR_DIRECTIONAL 0x0002 // infinite light source (alias)
#define LIGHTV1_ATTR_POINT 0x0004		// point light source
#define LIGHTV1_ATTR_SPOTLIGHT1 0x0008  // spotlight type 1 (simple)
#define LIGHTV1_ATTR_SPOTLIGHT2 0x0010  // spotlight type 2 (complex


#define MAX_LIGHTS 8 // good luck with 1!
#define _RGBA32BIT(r, g, b, a) ((a) + ((b) << 8) + ((g) << 16) + ((r) << 24))

#define SORT_POLYLIST_AVGZ 0  // sorts on average of all vertices
#define SORT_POLYLIST_NEARZ 1 // sorts on closest z vertex of each poly
#define SORT_POLYLIST_FARZ 2  // sorts on farthest z vertex of each poly

#define PARSER_DEBUG_OFF // enables/disables conditional compilation

#define PARSER_STRIP_EMPTY_LINES 1 // strips all blank lines
#define PARSER_LEAVE_EMPTY_LINES 2 // leaves empty lines
#define PARSER_STRIP_WS_ENDS 4	 // strips ws space at ends of line
#define PARSER_LEAVE_WS_ENDS 8	 // leaves it
#define PARSER_STRIP_COMMENTS 16   // strips comments out
#define PARSER_LEAVE_COMMENTS 32   // leaves comments in

#define PARSER_BUFFER_SIZE 256 // size of parser line buffer
#define PARSER_MAX_COMMENT 16  // maximum size of comment delimeter string

#define PARSER_DEFAULT_COMMENT "#" // default comment string for parser

// pattern language
#define PATTERN_TOKEN_FLOAT 'f'
#define PATTERN_TOKEN_INT 'i'
#define PATTERN_TOKEN_STRING 's'
#define PATTERN_TOKEN_LITERAL '\''

// state machine defines for pattern matching
#define PATTERN_STATE_INIT 0

#define PATTERN_STATE_RESTART 1
#define PATTERN_STATE_FLOAT 2
#define PATTERN_STATE_INT 3
#define PATTERN_STATE_LITERAL 4
#define PATTERN_STATE_STRING 5
#define PATTERN_STATE_NEXT 6

#define PATTERN_STATE_MATCH 7
#define PATTERN_STATE_END 8

#define PATTERN_MAX_ARGS 16
#define PATTERN_BUFFER_SIZE 80

// defines for objects version 2
// objects use dynamic allocation now, but keep as max values
#define OBJECT4DV2_MAX_VERTICES 4096 // 64
#define OBJECT4DV2_MAX_POLYS 8192	// 128

// states for objects
#define OBJECT4DV2_STATE_NULL 0x0000
#define OBJECT4DV2_STATE_ACTIVE 0x0001
#define OBJECT4DV2_STATE_VISIBLE 0x0002
#define OBJECT4DV2_STATE_CULLED 0x0004

// new
#define OBJECT4DV2_ATTR_SINGLE_FRAME 0x0001 // single frame object (emulates ver 1.0)
#define OBJECT4DV2_ATTR_MULTI_FRAME 0x0002  // multi frame object for .md2 support etc.
#define OBJECT4DV2_ATTR_TEXTURES 0x0004		// flags if object contains textured polys?

// render list defines ver 2.0
#define RENDERLIST4DV2_MAX_POLYS 32768

// defines for vertices, these are "hints" to the transform and
// lighting systems to help determine if a particular vertex has
// a valid normal that must be rotated, or a texture coordinate
// that must be clipped etc., this helps us minmize load during lighting
// and rendering since we can determine exactly what kind of vertex we
// are dealing with, something like a (direct3d) flexible vertex format in
// as much as it can hold:
// point
// point + normal
// point + normal + texture coordinates
#define VERTEX4DTV1_ATTR_NULL 0x0000 // this vertex is empty
#define VERTEX4DTV1_ATTR_POINT 0x0001
#define VERTEX4DTV1_ATTR_NORMAL 0x0002
#define VERTEX4DTV1_ATTR_TEXTURE 0x0004

// these are some defines for conditional compilation of the new rasterizers
// I don't want 80 million different functions, so I have decided to
// use some conditionals to change some of the logic in each
// these names aren't necessarily the most accurate, but 3 should be enough
#define RASTERIZER_ACCURATE 0 // sub-pixel accurate with fill convention
#define RASTERIZER_FAST 1	 //
#define RASTERIZER_FASTEST 2

// set this to the mode you want the engine to use
#define RASTERIZER_MODE RASTERIZER_ACCURATE

#define VERTEX_FLAGS_INVERT_X 0x0001 // inverts the Z-coordinates
#define VERTEX_FLAGS_INVERT_Y 0x0002 // inverts the Z-coordinates
#define VERTEX_FLAGS_INVERT_Z 0x0004 // inverts the Z-coordinates
#define VERTEX_FLAGS_SWAP_YZ 0x0008  // transforms a RHS model to a LHS model
#define VERTEX_FLAGS_SWAP_XZ 0x0010
#define VERTEX_FLAGS_SWAP_XY 0x0020
#define VERTEX_FLAGS_INVERT_WINDING_ORDER 0x0040 // invert winding order from cw to ccw or ccw to cc

#define VERTEX_FLAGS_TRANSFORM_LOCAL 0x0200		  // if file format has local transform then do it!
#define VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD 0x0400 // if file format has local to world then do it!


// states of polygons and faces
#define POLY4DV2_STATE_NULL 0x0000
#define POLY4DV2_STATE_ACTIVE 0x0001
#define POLY4DV2_STATE_CLIPPED 0x0002
#define POLY4DV2_STATE_BACKFACE 0x0004
#define POLY4DV2_STATE_LIT 0x0008

#define POLY4DV2_ATTR_2SIDED 0x0001
#define POLY4DV2_ATTR_TRANSPARENT 0x0002
#define POLY4DV2_ATTR_8BITCOLOR 0x0004
#define POLY4DV2_ATTR_RGB16 0x0008
#define POLY4DV2_ATTR_RGB24 0x0010

#define POLY4DV2_ATTR_SHADE_MODE_PURE 0x0020
#define POLY4DV2_ATTR_SHADE_MODE_CONSTANT 0x0020 // (alias)
#define POLY4DV2_ATTR_SHADE_MODE_EMISSIVE 0x0020 // (alias)

#define POLY4DV2_ATTR_SHADE_MODE_FLAT 0x0040
#define POLY4DV2_ATTR_SHADE_MODE_GOURAUD 0x0080
#define POLY4DV2_ATTR_SHADE_MODE_PHONG 0x0100
#define POLY4DV2_ATTR_SHADE_MODE_FASTPHONG 0x0100 // (alias)
#define POLY4DV2_ATTR_SHADE_MODE_TEXTURE 0x0200

// new
#define POLY4DV2_ATTR_ENABLE_MATERIAL 0x0800  // use a real material for lighting
#define POLY4DV2_ATTR_DISABLE_MATERIAL 0x1000 // use basic color only for lighting (emulate version 1.0)



// defines for materials, follow our polygon attributes as much as possible
#define MATV1_ATTR_2SIDED 0x0001
#define MATV1_ATTR_TRANSPARENT 0x0002
#define MATV1_ATTR_8BITCOLOR 0x0004
#define MATV1_ATTR_RGB16 0x0008
#define MATV1_ATTR_RGB24 0x0010

#define MATV1_ATTR_SHADE_MODE_CONSTANT 0x0020
#define MATV1_ATTR_SHADE_MODE_EMMISIVE 0x0020 // alias
#define MATV1_ATTR_SHADE_MODE_FLAT 0x0040
#define MATV1_ATTR_SHADE_MODE_GOURAUD 0x0080
#define MATV1_ATTR_SHADE_MODE_FASTPHONG 0x0100
#define MATV1_ATTR_SHADE_MODE_TEXTURE 0x0200

#define MAX_MATERIALS 256

// storage for our lookup tables
extern float cos_look[361]; // 1 extra so we can store 0-360 inclusive
extern float sin_look[361]; // 1 extra so we can store 0-360 inclusive

#include <iostream>
using namespace std;
//std::string testPrint(std::string inputStr);

//颜色转换
void Color2RGB(const int &hex, int &r, int &g, int &b);
void RGB2Color(int &hex, const int &r, const int &g, const int &b);

int testPrint();

// RGB+alpha color
typedef struct RGBAV1_TYP
{
	union {
		int rgba;		 // compressed format
		UCHAR rgba_M[4]; // array format
		struct
		{
			UCHAR a, b, g, r;
		}; // explict name format
	};	 // end union

} RGBAV1, *RGBAV1_PTR;

typedef struct BITMAP_IMAGE_TYP
{
	int state;		   // state of bitmap
	int attr;		   // attributes of bitmap
	int x, y;		   // position of bitmap
	int width, height; // size of bitmap
	int num_bytes;	 // total bytes of bitmap
	int bpp;		   // bits per pixel
	UCHAR *buffer;	 // pixels of bitmap

} BITMAP_IMAGE, *BITMAP_IMAGE_PTR;

typedef struct MATV1_TYP
{
	int state;	 // state of material
	int id;		   // id of this material, index into material array
	char name[64]; // name of material
	int attr;	  // attributes, the modes for shading, constant, flat,
				   // gouraud, fast phong, environment, textured etc.
				   // and other special flags...

	RGBAV1 color;			 // color of material
	float ka, kd, ks, power; // ambient, diffuse, specular,
							 // coefficients, note they are
							 // separate and scalars since many
							 // modelers use this format
							 // along with specular power

	RGBAV1 ra, rd, rs; // the reflectivities/colors pre-
					   // multiplied, to more match our
					   // definitions, each is basically
					   // computed by multiplying the
					   // color by the k's, eg:
					   // rd = color*kd etc.

	char texture_file[80]; // file location of texture
	BITMAP_IMAGE texture;  // actual texture map (if any)

	int iaux1, iaux2; // auxiliary vars for future expansion
	float faux1, faux2;
	void *ptr;

} MATV1, *MATV1_PTR;

extern MATV1 materials[MAX_MATERIALS]; // materials in system
extern int num_materials;			   // current number of materials

// 3D vector, point without the w ////////////////////////
typedef struct VECTOR3D_TYP
{
    union {
        float M[3]; // array indexed storage
        
        // explicit names
        struct
        {
            float x, y, z;
        }; // end struct
        
    }; // end union
    
} VECTOR3D, POINT3D, *VECTOR3D_PTR, *POINT3D_PTR;

typedef struct VECTOR4D_TYP
{
    union {
        float M[4];
        
        struct
        {
            float x, y, z, w;
        };
    };
} VECTOR4D, POINT4D, *VECTOR4D_PTR, *POINT4D_PTR;

//Polygon多边形，其实就是一个三角形
typedef struct POLYF4DV1_TYP
{
    int state; // state information
    int attr;  // physical attributes of polygon
    int color; // color of polygon
    
    POINT4D vlist[3];  // the vertices of this triangle
    POINT4D tvlist[3]; // the vertices after transformation if needed
    
    POLYF4DV1_TYP *next; // pointer to next polygon in list??
    POLYF4DV1_TYP *prev; // pointer to previous polygon in list??
    
} POLYF4DV1, *POLYF4DV1_PTR;

//渲染列表，其实就是三角形数组
typedef struct RENDERLIST4DV1_TYP
{
    int state; // state of renderlist ???
    int attr;  // attributes of renderlist ???
    
    POLYF4DV1_PTR poly_ptrs[RENDERLIST4DV1_MAX_POLYS];
    POLYF4DV1 poly_data[RENDERLIST4DV1_MAX_POLYS];
    
    int num_polys; // number of polys in render list
    
} RENDERLIST4DV1, *RENDERLIST4DV1_PTR;

// 3D 平面 ///////////////////////////////////////////////////
typedef struct PLANE3D_TYP
{
    POINT3D p0; // point on the plane
    VECTOR3D n; // normal to the plane (not necessarily a unit vector)
} PLANE3D, *PLANE3D_PTR;

//矩阵
// 4x4 matrix /////////////////////////////////////////////
typedef struct MATRIX4X4
{
    union {
        float M[4][4]; // array indexed data storage
        
        // storage in row major form with explicit names
        struct
        {
            float M00, M01, M02, M03;
            float M10, M11, M12, M13;
            float M20, M21, M22, M23;
            float M30, M31, M32, M33;
        }; // end explicit names
        
    }; // end union
    
} MATRIX4X4, *MATRIX4X4_PTR;

//矩阵

// 4x4 identity matrix
const MATRIX4X4 IMAT_4X4 = {1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1};

// camera version 1
typedef struct CAM4DV1_TYP
{
    int state; // state of camera
    int attr;  // camera attributes
    
    POINT4D pos; // world position of camera used by both camera models
    
    VECTOR4D dir; // angles or look at direction of camera for simple
    // euler camera models, elevation and heading for
    // uvn model
    
    VECTOR4D u; // extra vectors to track the camera orientation
    VECTOR4D v; // for more complex UVN camera model
    VECTOR4D n;
    
    VECTOR4D target; // look at target
    
    float view_dist; // focal length
    
    float fov; // field of view for both horizontal and vertical axes
    
    // 3d clipping planes
    // if view volume is NOT 90 degree then general 3d clipping
    // must be employed
    float near_clip_z; // near z=constant clipping plane
    float far_clip_z;  // far z=constant clipping plane
    
    PLANE3D rt_clip_plane; // the right clipping plane
    PLANE3D lt_clip_plane; // the left clipping plane
    PLANE3D tp_clip_plane; // the top clipping plane
    PLANE3D bt_clip_plane; // the bottom clipping plane
    
    float viewplane_width;  // width and height of view plane to project onto
    float viewplane_height; // usually 2x2 for normalized projection or
    // the exact same size as the viewport or screen window
    
    // remember screen and viewport are synonomous
    float viewport_width; // size of screen/viewport
    float viewport_height;
    float viewport_center_x; // center of view port (final image destination)
    float viewport_center_y;
    
    // aspect ratio
    float aspect_ratio;
    
    // these matrices are not necessarily needed based on the method of
    // transformation, for example, a manual perspective or screen transform
    // and or a concatenated perspective/screen, however, having these
    // matrices give us more flexibility
    
    MATRIX4X4 mcam; // storage for the world to camera transform matrix
    MATRIX4X4 mper; // storage for the camera to perspective transform matrix
    MATRIX4X4 mscr; // storage for the perspective to screen transform matrix
    
} CAM4DV1, *CAM4DV1_PTR;




// 基于外部顶点列表的多边形
typedef struct POLY4DV1_TYP
{
	int state; // 状态
	int attr;  // 属性
	int color; // 颜色

	POINT4D_PTR vlist; // 顶点指针
	int vert[3];	   // 三个顶点的索引
} POLY4DV1, *POLY4DV1_PTR;

//物体
typedef struct OBJECT4DV1_TYP
{
	int id;
	char name[64];
	int state;
	int attr;
	float avg_radius; // average radius of object used for collision detection
	float max_radius; // maximum radius of object

	POINT4D world_pos;

	VECTOR4D dir;

	VECTOR4D ux, uy, uz;

	int num_vertices;

	POINT4D vlist_local[OBJECT4DV1_MAX_VERTICES]; //局部顶点列表
	POINT4D vlist_trans[OBJECT4DV1_MAX_VERTICES]; //变换后的顶点列表

	int num_polys;
	POLY4DV1 plist[OBJECT4DV1_MAX_POLYS]; //多边形数组

} OBJECT4DV1, *OBJECT4DV1_PTR;





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

#define MAT_COPY_4X4(src_mat, dest_mat)                                   \
{                                                                     \
memcpy((void *)(dest_mat), (void *)(src_mat), sizeof(MATRIX4X4)); \
}

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
