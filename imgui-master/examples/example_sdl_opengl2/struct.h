#pragma once

#include "predefine.h"



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

typedef struct VECTOR2D_TYP
{
	union {
		float M[2]; // array indexed storage

		// explicit names
		struct
		{
			float x, y;
		}; // end struct

	}; // end union

} VECTOR2D, POINT2D, *VECTOR2D_PTR, *POINT2D_PTR;

typedef struct VERTEX4DTV1_TYP
{
	union {
		float M[12]; // array indexed storage

		// explicit names
		struct
		{
			float x, y, z, w;	 // point
			float nx, ny, nz, nw; // normal (vector or point)
			float u0, v0;		  // texture coordinates

			float i;  // final vertex intensity after lighting
			int attr; // attributes/ extra texture coordinates
		};			  // end struct

		// high level types
		struct
		{
			POINT4D v;  // the vertex
			VECTOR4D n; // the normal
			POINT2D t;  // texture coordinates
		};

	}; // end union

} VERTEX4DTV1, *VERTEX4DTV1_PTR;

// a polygon ver 2.0 based on an external vertex list  //////////////////////////////////
typedef struct POLY4DV2_TYP
{
	int state;		  // state information
	int attr;		  // physical attributes of polygon
	int color;		  // color of polygon
	int lit_color[3]; // holds colors after lighting, 0 for flat shading
					  // 0,1,2 for vertex colors after vertex lighting

	BITMAP_IMAGE_PTR texture; // pointer to the texture information for simple texture mapping

	int mati; // material index (-1) no material (new)

	VERTEX4DTV1_PTR vlist; // the vertex list itself
	POINT2D_PTR tlist;	 // the texture list itself (new)
	int vert[3];		   // the indices into the vertex list
	int text[3];		   // the indices into the texture coordinate list (new)
	float nlength;		   // length of normal (new)

} POLY4DV2, *POLY4DV2_PTR;


typedef struct OBJECT4DV2_TYP
{
	int id;			   // numeric id of this object
	char name[64];	 // ASCII name of object just for kicks
	int state;		   // state of object
	int attr;		   // attributes of object
	int mati;		   // material index overide (-1) - no material (new)
	float *avg_radius; // [OBJECT4DV2_MAX_FRAMES];   // average radius of object used for collision detection
	float *max_radius; // [OBJECT4DV2_MAX_FRAMES];   // maximum radius of object

	POINT4D world_pos; // position of object in world

	VECTOR4D dir; // rotation angles of object in local
				  // cords or unit direction vector user defined???

	VECTOR4D ux, uy, uz; // local axes to track full orientation
						 // this is updated automatically during
						 // rotation calls

	int num_vertices;   // number of vertices per frame of this object
	int num_frames;		// number of frames
	int total_vertices; // total vertices, redudant, but it saves a multiply in a lot of places
	int curr_frame;		// current animation frame (0) if single frame

	VERTEX4DTV1_PTR vlist_local; // [OBJECT4DV1_MAX_VERTICES]; // array of local vertices
	VERTEX4DTV1_PTR vlist_trans; // [OBJECT4DV1_MAX_VERTICES]; // array of transformed vertices

	// these are needed to track the "head" of the vertex list for mult-frame objects
	VERTEX4DTV1_PTR head_vlist_local;
	VERTEX4DTV1_PTR head_vlist_trans;

	// texture coordinates list (new)
	POINT2D_PTR tlist; // 3*num polys at max

	BITMAP_IMAGE_PTR texture; // pointer to the texture information for simple texture mapping (new)

	int num_polys;		// number of polygons in object mesh
	POLY4DV2_PTR plist; // ptr to polygons (new)

	int ivar1, ivar2;   // auxiliary vars
	float fvar1, fvar2; // auxiliary vars

	// METHODS //////////////////////////////////////////////////

	// setting the frame is so important that it should be a member function
	// calling functions without doing this can wreak havok!
	int Set_Frame(int frame);

} OBJECT4DV2, *OBJECT4DV2_PTR;

// first light structure
typedef struct LIGHTV1_TYP
{
	int state; // state of light
	int id;	// id of light
	int attr;  // type of light, and extra qualifiers

	RGBAV1 c_ambient;  // ambient light intensity
	RGBAV1 c_diffuse;  // diffuse light intensity
	RGBAV1 c_specular; // specular light intensity

	POINT4D pos;	  // position of light
	VECTOR4D dir;	 // direction of light
	float kc, kl, kq; // attenuation factors
	float spot_inner; // inner angle for spot light
	float spot_outer; // outer angle for spot light
	float pf;		  // power factor/falloff for spot lights

	int iaux1, iaux2; // auxiliary vars for future expansion
	float faux1, faux2;
	void *ptr;

} LIGHTV1, *LIGHTV1_PTR;

typedef struct POLYF4DV2_TYP
{
	int state;				  // state information
	int attr;				  // physical attributes of polygon
	int color;				  // color of polygon
	int lit_color[3];		  // holds colors after lighting, 0 for flat shading
							  // 0,1,2 for vertex colors after vertex lighting
	BITMAP_IMAGE_PTR texture; // pointer to the texture information for simple texture mapping

	int mati; // material index (-1) for no material  (new)

	float nlength;   // length of the polygon normal if not normalized (new)
	VECTOR4D normal; // the general polygon normal (new)

	float avg_z; // average z of vertices, used for simple sorting (new)

	VERTEX4DTV1 vlist[3];  // the vertices of this triangle
	VERTEX4DTV1 tvlist[3]; // the vertices after transformation if needed

	POLYF4DV2_TYP *next; // pointer to next polygon in list??
	POLYF4DV2_TYP *prev; // pointer to previous polygon in list??

} POLYF4DV2, *POLYF4DV2_PTR;

typedef struct RENDERLIST4DV2_TYP
{
	int state; // state of renderlist ???
	int attr;  // attributes of renderlist ???

	// the render list is an array of pointers each pointing to
	// a self contained "renderable" polygon face POLYF4DV2
	POLYF4DV2_PTR poly_ptrs[RENDERLIST4DV2_MAX_POLYS];

	// additionally to cut down on allocatation, de-allocation
	// of polygons each frame, here's where the actual polygon
	// faces will be stored
	POLYF4DV2 poly_data[RENDERLIST4DV2_MAX_POLYS];

	int num_polys; // number of polys in render list

} RENDERLIST4DV2, *RENDERLIST4DV2_PTR;


// a 2D vertex
typedef struct VERTEX2DF_TYP
{
	float x, y; // the vertex
} VERTEX2DF, *VERTEX2DF_PTR;