#pragma once

#include <math.h>

typedef unsigned short USHORT;
typedef unsigned char UCHAR;
typedef unsigned int IUINT32;


#define FLT_MIN 1.175494351e-38F 
#define FLT_MAX 3.402823466e+38F 

#define INT_MAX 2147483647
#define INT_MIN (-INT_MAX - 1)


#define _RGB565FROM16BIT(RGB, r,g,b) { *r = ( ((RGB) >> 11) & 0x1f); *g = (((RGB) >> 5) & 0x3f); *b = ((RGB) & 0x1f); 
#define _RGB16BIT565(r, g, b) ((b & 31) + ((g & 63) << 5) + ((r & 31) << 11))

#define EPSILON_E3 (float)(1E-3)
#define EPSILON_E4 (float)(1E-4)
#define EPSILON_E5 (float)(1E-5)
#define EPSILON_E6 (float)(1E-6)
#define FCMP(a,b) ( (fabs(a-b) < EPSILON_E3) ? 1 : 0)
#define SWAP(a,b,t) {t=a; a=b; b=t;}

#define WINDOW_WIDTH 1000 // size of window
#define WINDOW_HEIGHT 1000

#define min_clip_x  0 // clipping rectangle
#define max_clip_x  (WINDOW_WIDTH - 1)
#define min_clip_y  0
#define max_clip_y  (WINDOW_HEIGHT - 1)

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

// #define WINDOW_WIDTH 400 // size of window
// #define WINDOW_HEIGHT 400

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

#define MAT_COPY_4X4(src_mat, dest_mat)                                   \
{                                                                     \
memcpy((void *)(dest_mat), (void *)(src_mat), sizeof(MATRIX4X4)); \
}


#define VERTEX_FLAGS_INVERT_TEXTURE_U 0x0080 // invert u texture coordinate
#define VERTEX_FLAGS_INVERT_TEXTURE_V 0x0100 // invert v texture coordinate
#define VERTEX_FLAGS_INVERT_SWAP_UV 0x0800   // swap u and v texture coordinate
