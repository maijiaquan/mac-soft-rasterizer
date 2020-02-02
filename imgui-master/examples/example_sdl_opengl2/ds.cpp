#include "ds.h"

//extern 变量
int num_lights;             // current number of lights
LIGHTV1 lights[MAX_LIGHTS]; // lights in system

// storage for our lookup tables
float cos_look[361]; // 1 extra element so we can store 0-360 inclusive
float sin_look[361];

char texture_path[80] = "./"; // root path to ALL textures, make current directory for now
// BITMAP_FILE bitmap16bit; // a 16 bit bitmap file

MATV1 materials[MAX_MATERIALS]; // materials in system
int num_materials;              // current number of materials

USHORT RGB16Bit565(int r, int g, int b)
{
    // this function simply builds a 5.6.5 format 16 bit pixel
    // assumes input is RGB 0-255 each channel
    r >>= 3;
    g >>= 2;
    b >>= 3;
    return (_RGB16BIT565((r), (g), (b)));

} // end RGB16Bit565


void Color2RGB(const int &c, int &r, int &g, int &b)
{
    r = (0xff << 16 & c) >> 16;
    g = (0xff << 8 & c) >> 8;
    b = 0xff & c;
}

void RGB2Color(int &c, const int &r, const int &g, const int &b)
{
    c = (r << 16) | (g << 8) | b;
}

void Reset_RENDERLIST4DV1(RENDERLIST4DV1_PTR renderList)
{
    renderList->num_polys = 0;
}

void Init_CAM4DV1(CAM4DV1_PTR cam,        // the camera object
                  int cam_attr,           // attributes
                  POINT4D_PTR cam_pos,    // initial camera position
                  VECTOR4D_PTR cam_dir,   // initial camera angles
                  POINT4D_PTR cam_target, // UVN target
                  float near_clip_z,      // near and far clipping planes
                  float far_clip_z,
                  float fov,            // field of view in degrees
                  float viewport_width, // size of final screen viewport
                  float viewport_height)
{
    // this function initializes the camera object cam, the function
    // doesn't do a lot of error checking or sanity checking since
    // I want to allow you to create projections as you wish, also
    // I tried to minimize the number of parameters the functions needs

    // first set up parms that are no brainers
    cam->attr = cam_attr; // camera attributes

    VECTOR4D_COPY(&cam->pos, cam_pos); // positions
    VECTOR4D_COPY(&cam->dir, cam_dir); // direction vector or angles for
    // euler camera
    // for UVN camera
    VECTOR4D_INITXYZ(&cam->u, 1, 0, 0); // set to +x
    VECTOR4D_INITXYZ(&cam->v, 0, 1, 0); // set to +y
    VECTOR4D_INITXYZ(&cam->n, 0, 0, 1); // set to +z

    if (cam_target != nullptr)
        VECTOR4D_COPY(&cam->target, cam_target); // UVN target
    else
        VECTOR4D_ZERO(&cam->target);

    cam->near_clip_z = near_clip_z; // near z=constant clipping plane
    cam->far_clip_z = far_clip_z;   // far z=constant clipping plane

    cam->viewport_width = viewport_width; // dimensions of viewport
    cam->viewport_height = viewport_height;

    cam->viewport_center_x = (viewport_width - 1) / 2; // center of viewport
    cam->viewport_center_y = (viewport_height - 1) / 2;

    cam->aspect_ratio = (float)viewport_width / (float)viewport_height;

    // set all camera matrices to identity matrix
    MAT_IDENTITY_4X4(&cam->mcam);
    MAT_IDENTITY_4X4(&cam->mper);
    MAT_IDENTITY_4X4(&cam->mscr);

    // set independent vars
    cam->fov = fov;

    // set the viewplane dimensions up, they will be 2 x (2/ar)
    cam->viewplane_width = 2.0;
    cam->viewplane_height = 2.0 / cam->aspect_ratio;

    // now we know fov and we know the viewplane dimensions plug into formula and
    // solve for view distance parameters
    float tan_fov_div2 = tan(DEG_TO_RAD(fov / 2));

    cam->view_dist = (0.5) * (cam->viewplane_width) * tan_fov_div2;

    // test for 90 fov first since it's easy :)
    if (fov == 90.0)
    {
        // set up the clipping planes -- easy for 90 degrees!
        POINT3D pt_origin; // point on the plane
        VECTOR3D_INITXYZ(&pt_origin, 0, 0, 0);

        VECTOR3D vn; // normal to plane

        // right clipping plane
        VECTOR3D_INITXYZ(&vn, 1, 0, -1); // x=z plane
        PLANE3D_Init(&cam->rt_clip_plane, &pt_origin, &vn, 1);

        // left clipping plane
        VECTOR3D_INITXYZ(&vn, -1, 0, -1); // -x=z plane
        PLANE3D_Init(&cam->lt_clip_plane, &pt_origin, &vn, 1);

        // top clipping plane
        VECTOR3D_INITXYZ(&vn, 0, 1, -1); // y=z plane
        PLANE3D_Init(&cam->tp_clip_plane, &pt_origin, &vn, 1);

        // bottom clipping plane
        VECTOR3D_INITXYZ(&vn, 0, -1, -1); // -y=z plane
        PLANE3D_Init(&cam->bt_clip_plane, &pt_origin, &vn, 1);
    } // end if d=1
    else
    {
        // now compute clipping planes yuck!
        POINT3D pt_origin; // point on the plane
        VECTOR3D_INITXYZ(&pt_origin, 0, 0, 0);

        VECTOR3D vn; // normal to plane

        // since we don't have a 90 fov, computing the normals
        // are a bit tricky, there are a number of geometric constructions
        // that solve the problem, but I'm going to solve for the
        // vectors that represent the 2D projections of the frustrum planes
        // on the x-z and y-z planes and then find perpendiculars to them

        // right clipping plane, check the math on graph paper
        VECTOR3D_INITXYZ(&vn, cam->view_dist, 0, -cam->viewplane_width / 2.0);
        PLANE3D_Init(&cam->rt_clip_plane, &pt_origin, &vn, 1);

        // left clipping plane, we can simply reflect the right normal about
        // the z axis since the planes are symetric about the z axis
        // thus invert x only
        VECTOR3D_INITXYZ(&vn, -cam->view_dist, 0, -cam->viewplane_width / 2.0);
        PLANE3D_Init(&cam->lt_clip_plane, &pt_origin, &vn, 1);

        // top clipping plane, same construction
        VECTOR3D_INITXYZ(&vn, 0, cam->view_dist, -cam->viewplane_width / 2.0);
        PLANE3D_Init(&cam->tp_clip_plane, &pt_origin, &vn, 1);

        // bottom clipping plane, same inversion
        VECTOR3D_INITXYZ(&vn, 0, -cam->view_dist, -cam->viewplane_width / 2.0);
        PLANE3D_Init(&cam->bt_clip_plane, &pt_origin, &vn, 1);
    } // end else

} // end Init_CAM4DV1

void MAT_IDENTITY_4X4(MATRIX4X4_PTR m)
{
    memcpy((void *)(m), (void *)&IMAT_4X4, sizeof(MATRIX4X4));
}

void PLANE3D_Init(PLANE3D_PTR plane, POINT3D_PTR p0,
                  VECTOR3D_PTR normal, int normalize = 0)
{
    // this function initializes a 3d plane

    // copy the point
    POINT3D_COPY(&plane->p0, p0);

    // if normalize is 1 then the normal is made into a unit vector
    if (!normalize)
        VECTOR3D_COPY(&plane->n, normal);
    else
    {
        // make normal into unit vector
        VECTOR3D_Normalize(normal, &plane->n);
    } // end else
}

void VECTOR3D_Normalize(VECTOR3D_PTR va)
{
    // normalizes the sent vector in placew

    // compute length
    float length = sqrtf(va->x * va->x + va->y * va->y + va->z * va->z);

    // test for zero length vector
    // if found return zero vector
    if (length < EPSILON_E5)
        return;

    float length_inv = 1 / length;

    // compute normalized version of vector
    va->x *= length_inv;
    va->y *= length_inv;
    va->z *= length_inv;

} // end VECTOR3D_Normalize

void VECTOR3D_Normalize(VECTOR3D_PTR va, VECTOR3D_PTR vn)
{
    // normalizes the sent vector and returns the result in vn

    VECTOR3D_ZERO(vn);

    // compute length
    float length = VECTOR3D_Length(va);

    // test for zero length vector
    // if found return zero vector
    if (length < EPSILON_E5)
        return;

    float length_inv = 1.0 / length;

    // compute normalized version of vector
    vn->x = va->x * length_inv;
    vn->y = va->y * length_inv;
    vn->z = va->z * length_inv;

} // end VECTOR3D_Normalize

float VECTOR3D_Length(VECTOR3D_PTR va)
{
    // computes the magnitude of a vector, slow

    return ((float)sqrtf(va->x * va->x + va->y * va->y + va->z * va->z));

} // end VECTOR3D_Length

int Insert_POLYF4DV1_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,
                                    POLYF4DV1_PTR poly)
{
    // inserts the sent polyface POLYF4DV1 into the render list

    // step 0: are we full?
    if (rend_list->num_polys >= RENDERLIST4DV1_MAX_POLYS)
        return (0);

    // step 1: copy polygon into next opening in polygon render list

    // point pointer to polygon structure
    rend_list->poly_ptrs[rend_list->num_polys] = &rend_list->poly_data[rend_list->num_polys];

    // copy face right into array, thats it
    memcpy((void *)&rend_list->poly_data[rend_list->num_polys], (void *)poly, sizeof(POLYF4DV1));

    // now the polygon is loaded into the next free array position, but
    // we need to fix up the links
    // test if this is the first entry
    if (rend_list->num_polys == 0)
    {
        // set pointers to null, could loop them around though to self
        rend_list->poly_data[0].next = NULL;
        rend_list->poly_data[0].prev = NULL;
    } // end if
    else
    {
        // first set this node to point to previous node and next node (null)
        rend_list->poly_data[rend_list->num_polys].next = NULL;
        rend_list->poly_data[rend_list->num_polys].prev =
            &rend_list->poly_data[rend_list->num_polys - 1];

        // now set previous node to point to this node
        rend_list->poly_data[rend_list->num_polys - 1].next =
            &rend_list->poly_data[rend_list->num_polys];
    } // end else

    // increment number of polys in list
    rend_list->num_polys++;

    // return successful insertion
    return (1);

} // end Insert_POLYF4DV1_RENDERLIST4DV1

void Build_XYZ_Rotation_MATRIX4X4(float theta_x, // euler angles
                                  float theta_y,
                                  float theta_z,
                                  MATRIX4X4_PTR mrot) // output
{
    // this helper function takes a set if euler angles and computes
    // a rotation matrix from them, usefull for object and camera
    // work, also  we will do a little testing in the function to determine
    // the rotations that need to be performed, since there's no
    // reason to perform extra matrix multiplies if the angles are
    // zero!

    MATRIX4X4 mx, my, mz, mtmp;         // working matrices
    float sin_theta = 0, cos_theta = 0; // used to initialize matrices
    int rot_seq = 0;                    // 1 for x, 2 for y, 4 for z

    // step 0: fill in with identity matrix
    MAT_IDENTITY_4X4(mrot);

    // step 1: based on zero and non-zero rotation angles, determine
    // rotation sequence
    if (fabs(theta_x) > EPSILON_E5) // x
        rot_seq = rot_seq | 1;

    if (fabs(theta_y) > EPSILON_E5) // y
        rot_seq = rot_seq | 2;

    if (fabs(theta_z) > EPSILON_E5) // z
        rot_seq = rot_seq | 4;

    // now case on sequence
    switch (rot_seq)
    {
    case 0: // no rotation
    {
        // what a waste!
        return;
    }
    break;

    case 1: // x rotation
    {
        // compute the sine and cosine of the angle
        cos_theta = Fast_Cos(theta_x);
        sin_theta = Fast_Sin(theta_x);

        // set the matrix up
        Mat_Init_4X4(&mx, 1, 0, 0, 0,
                     0, cos_theta, sin_theta, 0,
                     0, -sin_theta, cos_theta, 0,
                     0, 0, 0, 1);

        // that's it, copy to output matrix
        MAT_COPY_4X4(&mx, mrot);
        return;
    }
    break;

    case 2: // y rotation
    {
        // compute the sine and cosine of the angle
        cos_theta = Fast_Cos(theta_y);
        sin_theta = Fast_Sin(theta_y);

        // set the matrix up
        Mat_Init_4X4(&my, cos_theta, 0, -sin_theta, 0,
                     0, 1, 0, 0,
                     sin_theta, 0, cos_theta, 0,
                     0, 0, 0, 1);

        // that's it, copy to output matrix
        MAT_COPY_4X4(&my, mrot);
        return;
    }
    break;

    case 3: // xy rotation
    {
        // compute the sine and cosine of the angle for x
        cos_theta = Fast_Cos(theta_x);
        sin_theta = Fast_Sin(theta_x);

        // set the matrix up
        Mat_Init_4X4(&mx, 1, 0, 0, 0,
                     0, cos_theta, sin_theta, 0,
                     0, -sin_theta, cos_theta, 0,
                     0, 0, 0, 1);

        // compute the sine and cosine of the angle for y
        cos_theta = Fast_Cos(theta_y);
        sin_theta = Fast_Sin(theta_y);

        // set the matrix up
        Mat_Init_4X4(&my, cos_theta, 0, -sin_theta, 0,
                     0, 1, 0, 0,
                     sin_theta, 0, cos_theta, 0,
                     0, 0, 0, 1);

        // concatenate matrices
        Mat_Mul_4X4(&mx, &my, mrot);
        return;
    }
    break;

    case 4: // z rotation
    {
        // compute the sine and cosine of the angle
        cos_theta = Fast_Cos(theta_z);
        sin_theta = Fast_Sin(theta_z);

        // set the matrix up
        Mat_Init_4X4(&mz, cos_theta, sin_theta, 0, 0,
                     -sin_theta, cos_theta, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1);

        // that's it, copy to output matrix
        MAT_COPY_4X4(&mz, mrot);
        return;
    }
    break;

    case 5: // xz rotation
    {
        // compute the sine and cosine of the angle x
        cos_theta = Fast_Cos(theta_x);
        sin_theta = Fast_Sin(theta_x);

        // set the matrix up
        Mat_Init_4X4(&mx, 1, 0, 0, 0,
                     0, cos_theta, sin_theta, 0,
                     0, -sin_theta, cos_theta, 0,
                     0, 0, 0, 1);

        // compute the sine and cosine of the angle z
        cos_theta = Fast_Cos(theta_z);
        sin_theta = Fast_Sin(theta_z);

        // set the matrix up
        Mat_Init_4X4(&mz, cos_theta, sin_theta, 0, 0,
                     -sin_theta, cos_theta, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1);

        // concatenate matrices
        Mat_Mul_4X4(&mx, &mz, mrot);
        return;
    }
    break;

    case 6: // yz rotation
    {
        // compute the sine and cosine of the angle y
        cos_theta = Fast_Cos(theta_y);
        sin_theta = Fast_Sin(theta_y);

        // set the matrix up
        Mat_Init_4X4(&my, cos_theta, 0, -sin_theta, 0,
                     0, 1, 0, 0,
                     sin_theta, 0, cos_theta, 0,
                     0, 0, 0, 1);

        // compute the sine and cosine of the angle z
        cos_theta = Fast_Cos(theta_z);
        sin_theta = Fast_Sin(theta_z);

        // set the matrix up
        Mat_Init_4X4(&mz, cos_theta, sin_theta, 0, 0,
                     -sin_theta, cos_theta, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1);

        // concatenate matrices
        Mat_Mul_4X4(&my, &mz, mrot);
        return;
    }
    break;

    case 7: // xyz rotation
    {
        // compute the sine and cosine of the angle x
        cos_theta = Fast_Cos(theta_x);
        sin_theta = Fast_Sin(theta_x);

        // set the matrix up
        Mat_Init_4X4(&mx, 1, 0, 0, 0,
                     0, cos_theta, sin_theta, 0,
                     0, -sin_theta, cos_theta, 0,
                     0, 0, 0, 1);

        // compute the sine and cosine of the angle y
        cos_theta = Fast_Cos(theta_y);
        sin_theta = Fast_Sin(theta_y);

        // set the matrix up
        Mat_Init_4X4(&my, cos_theta, 0, -sin_theta, 0,
                     0, 1, 0, 0,
                     sin_theta, 0, cos_theta, 0,
                     0, 0, 0, 1);

        // compute the sine and cosine of the angle z
        cos_theta = Fast_Cos(theta_z);
        sin_theta = Fast_Sin(theta_z);

        // set the matrix up
        Mat_Init_4X4(&mz, cos_theta, sin_theta, 0, 0,
                     -sin_theta, cos_theta, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1);

        // concatenate matrices, watch order!
        Mat_Mul_4X4(&mx, &my, &mtmp);
        Mat_Mul_4X4(&mtmp, &mz, mrot);
    }
    break;

    default:
        break;

    } // end switch

} // end Build_XYZ_Rotation_MATRIX4X4

float Fast_Sin(float theta)
{
    // this function uses the sin_look[] lookup table, but
    // has logic to handle negative angles as well as fractional
    // angles via interpolation, use this for a more robust
    // sin computation that the blind lookup, but with with
    // a slight hit in speed

    // convert angle to 0-359
    theta = fmodf(theta, 360);

    // make angle positive
    if (theta < 0)
        theta += 360.0;

    // compute floor of theta and fractional part to interpolate
    int theta_int = (int)theta;
    float theta_frac = theta - theta_int;

    // now compute the value of sin(angle) using the lookup tables
    // and interpolating the fractional part, note that if theta_int
    // is equal to 359 then theta_int+1=360, but this is fine since the
    // table was made with the entries 0-360 inclusive
    return (sin_look[theta_int] +
            theta_frac * (sin_look[theta_int + 1] - sin_look[theta_int]));

} // end Fast_Sin

///////////////////////////////////////////////////////////////

float Fast_Cos(float theta)
{
    // this function uses the cos_look[] lookup table, but
    // has logic to handle negative angles as well as fractional
    // angles via interpolation, use this for a more robust
    // cos computation that the blind lookup, but with with
    // a slight hit in speed

    // convert angle to 0-359
    theta = fmodf(theta, 360);

    // make angle positive
    if (theta < 0)
        theta += 360.0;

    // compute floor of theta and fractional part to interpolate
    int theta_int = (int)theta;
    float theta_frac = theta - theta_int;

    // now compute the value of sin(angle) using the lookup tables
    // and interpolating the fractional part, note that if theta_int
    // is equal to 359 then theta_int+1=360, but this is fine since the
    // table was made with the entries 0-360 inclusive
    return (cos_look[theta_int] +
            theta_frac * (cos_look[theta_int + 1] - cos_look[theta_int]));

} // end Fast_Cos

void Mat_Init_4X4(MATRIX4X4_PTR ma,
                  float m00, float m01, float m02, float m03,
                  float m10, float m11, float m12, float m13,
                  float m20, float m21, float m22, float m23,
                  float m30, float m31, float m32, float m33)

{
    // this function fills a 4x4 matrix with the sent data in
    // row major form
    ma->M00 = m00;
    ma->M01 = m01;
    ma->M02 = m02;
    ma->M03 = m03;
    ma->M10 = m10;
    ma->M11 = m11;
    ma->M12 = m12;
    ma->M13 = m13;
    ma->M20 = m20;
    ma->M21 = m21;
    ma->M22 = m22;
    ma->M23 = m23;
    ma->M30 = m30;
    ma->M31 = m31;
    ma->M32 = m32;
    ma->M33 = m33;

} // end Mat_Init_4X4

void Mat_Mul_4X4(MATRIX4X4_PTR ma,
                 MATRIX4X4_PTR mb,
                 MATRIX4X4_PTR mprod)
{
    // this function multiplies two 4x4 matrices together and
    // and stores the result in mprod
    // note later we will take advantage of the fact that we know
    // that w=1 always, and that the last column of a 4x4 is
    // always 0

    for (int row = 0; row < 4; row++)
    {
        for (int col = 0; col < 4; col++)
        {
            // compute dot product from row of ma
            // and column of mb

            float sum = 0; // used to hold result

            for (int index = 0; index < 4; index++)
            {
                // add in next product pair
                sum += (ma->M[row][index] * mb->M[index][col]);
            } // end for index

            // insert resulting row,col element
            mprod->M[row][col] = sum;

        } // end for col

    } // end for row

} // end Mat_Mul_4X4

void Transform_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, // render list to transform
                              MATRIX4X4_PTR mt,             // transformation matrix
                              int coord_select)             // selects coords to transform
{
    // this function simply transforms all of the polygons vertices in the local or trans
    // array of the render list by the sent matrix

    // what coordinates should be transformed?
    switch (coord_select)
    {
    case TRANSFORM_LOCAL_ONLY:
    {
        for (int poly = 0; poly < rend_list->num_polys; poly++)
        {
            // acquire current polygon
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

            // is this polygon valid?
            // transform this polygon if and only if it's not clipped, not culled,
            // active, and visible, note however the concept of "backface" is
            // irrelevant in a wire frame engine though
            if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
                (curr_poly->state & POLY4DV1_STATE_BACKFACE))
                continue; // move onto next poly

            // all good, let's transform
            for (int vertex = 0; vertex < 3; vertex++)
            {
                // transform the vertex by mt
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&curr_poly->vlist[vertex], mt, &presult);

                // store result back
                VECTOR4D_COPY(&curr_poly->vlist[vertex], &presult);
            } // end for vertex

        } // end for poly
    }
    break;

    case TRANSFORM_TRANS_ONLY:
    {
        // transform each "transformed" vertex of the render list
        // remember, the idea of the tvlist[] array is to accumulate
        // transformations
        for (int poly = 0; poly < rend_list->num_polys; poly++)
        {
            // acquire current polygon
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

            // is this polygon valid?
            // transform this polygon if and only if it's not clipped, not culled,
            // active, and visible, note however the concept of "backface" is
            // irrelevant in a wire frame engine though
            if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
                (curr_poly->state & POLY4DV1_STATE_BACKFACE))
                continue; // move onto next poly

            // all good, let's transform
            for (int vertex = 0; vertex < 3; vertex++)
            {
                // transform the vertex by mt
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&curr_poly->tvlist[vertex], mt, &presult);

                // store result back
                VECTOR4D_COPY(&curr_poly->tvlist[vertex], &presult);
            } // end for vertex

        } // end for poly
    }
    break;

    case TRANSFORM_LOCAL_TO_TRANS:
    {
        // transform each local/model vertex of the render list and store result
        // in "transformed" vertex list
        for (int poly = 0; poly < rend_list->num_polys; poly++)
        {
            // acquire current polygon
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

            // is this polygon valid?
            // transform this polygon if and only if it's not clipped, not culled,
            // active, and visible, note however the concept of "backface" is
            // irrelevant in a wire frame engine though
            if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
                (curr_poly->state & POLY4DV1_STATE_BACKFACE))
                continue; // move onto next poly

            // all good, let's transform
            for (int vertex = 0; vertex < 3; vertex++)
            {
                // transform the vertex by mt
                Mat_Mul_VECTOR4D_4X4(&curr_poly->vlist[vertex], mt, &curr_poly->tvlist[vertex]);
            } // end for vertex

        } // end for poly
    }
    break;

    default:
        break;

    } // end switch

} // end Transform_RENDERLIST4DV1

void Model_To_World_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,
                                   POINT4D_PTR world_pos,
                                   int coord_select)
{
    // NOTE: Not matrix based
    // this function converts the local model coordinates of the
    // sent render list into world coordinates, the results are stored
    // in the transformed vertex list (tvlist) within the renderlist

    // interate thru vertex list and transform all the model/local
    // coords to world coords by translating the vertex list by
    // the amount world_pos and storing the results in tvlist[]
    // is this polygon valid?

    if (coord_select == TRANSFORM_LOCAL_TO_TRANS)
    {
        for (int poly = 0; poly < rend_list->num_polys; poly++)
        {
            // acquire current polygon
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

            // transform this polygon if and only if it's not clipped, not culled,
            // active, and visible, note however the concept of "backface" is
            // irrelevant in a wire frame engine though
            if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
                (curr_poly->state & POLY4DV1_STATE_BACKFACE))
                continue; // move onto next poly

            // all good, let's transform
            for (int vertex = 0; vertex < 3; vertex++)
            {
                // translate vertex
                VECTOR4D_Add(&curr_poly->vlist[vertex], world_pos, &curr_poly->tvlist[vertex]);
            } // end for vertex

        } // end for poly
    }     // end if local
    else  // TRANSFORM_TRANS_ONLY
    {
        for (int poly = 0; poly < rend_list->num_polys; poly++)
        {
            // acquire current polygon
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

            // transform this polygon if and only if it's not clipped, not culled,
            // active, and visible, note however the concept of "backface" is
            // irrelevant in a wire frame engine though
            if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
                (curr_poly->state & POLY4DV1_STATE_BACKFACE))
                continue; // move onto next poly

            for (int vertex = 0; vertex < 3; vertex++)
            {
                // translate vertex
                VECTOR4D_Add(&curr_poly->tvlist[vertex], world_pos, &curr_poly->tvlist[vertex]);
            } // end for vertex

        } // end for poly

    } // end else

} // end Model_To_World_RENDERLIST4DV1

void Build_CAM4DV1_Matrix_Euler(CAM4DV1_PTR cam, int cam_rot_seq)
{
    // this creates a camera matrix based on Euler angles
    // and stores it in the sent camera object
    // if you recall from chapter 6 to create the camera matrix
    // we need to create a transformation matrix that looks like:

    // Mcam = mt(-1) * my(-1) * mx(-1) * mz(-1)
    // that is the inverse of the camera translation matrix mutilplied
    // by the inverses of yxz, in that order, however, the order of
    // the rotation matrices is really up to you, so we aren't going
    // to force any order, thus its programmable based on the value
    // of cam_rot_seq which can be any value CAM_ROT_SEQ_XYZ where
    // XYZ can be in any order, YXZ, ZXY, etc.

    MATRIX4X4 mt_inv, // inverse camera translation matrix
        mx_inv,       // inverse camera x axis rotation matrix
        my_inv,       // inverse camera y axis rotation matrix
        mz_inv,       // inverse camera z axis rotation matrix
        mrot,         // concatenated inverse rotation matrices
        mtmp;         // temporary working matrix

    // step 1: create the inverse translation matrix for the camera
    // position
    Mat_Init_4X4(&mt_inv, 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1, 0,
                 -cam->pos.x, -cam->pos.y, -cam->pos.z, 1);

    // step 2: create the inverse rotation sequence for the camera
    // rember either the transpose of the normal rotation matrix or
    // plugging negative values into each of the rotations will result
    // in an inverse matrix

    // first compute all 3 rotation matrices

    // extract out euler angles
    float theta_x = cam->dir.x;
    float theta_y = cam->dir.y;
    float theta_z = cam->dir.z;

    // compute the sine and cosine of the angle x
    float cos_theta = Fast_Cos(theta_x);  // no change since cos(-x) = cos(x)
    float sin_theta = -Fast_Sin(theta_x); // sin(-x) = -sin(x)

    // set the matrix up
    Mat_Init_4X4(&mx_inv, 1, 0, 0, 0,
                 0, cos_theta, sin_theta, 0,
                 0, -sin_theta, cos_theta, 0,
                 0, 0, 0, 1);

    // compute the sine and cosine of the angle y
    cos_theta = Fast_Cos(theta_y);  // no change since cos(-x) = cos(x)
    sin_theta = -Fast_Sin(theta_y); // sin(-x) = -sin(x)

    // set the matrix up
    Mat_Init_4X4(&my_inv, cos_theta, 0, -sin_theta, 0,
                 0, 1, 0, 0,
                 sin_theta, 0, cos_theta, 0,
                 0, 0, 0, 1);

    // compute the sine and cosine of the angle z
    cos_theta = Fast_Cos(theta_z);  // no change since cos(-x) = cos(x)
    sin_theta = -Fast_Sin(theta_z); // sin(-x) = -sin(x)

    // set the matrix up
    Mat_Init_4X4(&mz_inv, cos_theta, sin_theta, 0, 0,
                 -sin_theta, cos_theta, 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1);

    // now compute inverse camera rotation sequence
    switch (cam_rot_seq)
    {
    case CAM_ROT_SEQ_XYZ:
    {
        Mat_Mul_4X4(&mx_inv, &my_inv, &mtmp);
        Mat_Mul_4X4(&mtmp, &mz_inv, &mrot);
    }
    break;

    case CAM_ROT_SEQ_YXZ:
    {
        Mat_Mul_4X4(&my_inv, &mx_inv, &mtmp);
        Mat_Mul_4X4(&mtmp, &mz_inv, &mrot);
    }
    break;

    case CAM_ROT_SEQ_XZY:
    {
        Mat_Mul_4X4(&mx_inv, &mz_inv, &mtmp);
        Mat_Mul_4X4(&mtmp, &my_inv, &mrot);
    }
    break;

    case CAM_ROT_SEQ_YZX:
    {
        Mat_Mul_4X4(&my_inv, &mz_inv, &mtmp);
        Mat_Mul_4X4(&mtmp, &mx_inv, &mrot);
    }
    break;

    case CAM_ROT_SEQ_ZYX:
    {
        Mat_Mul_4X4(&mz_inv, &my_inv, &mtmp);
        Mat_Mul_4X4(&mtmp, &mx_inv, &mrot);
    }
    break;

    case CAM_ROT_SEQ_ZXY:
    {
        Mat_Mul_4X4(&mz_inv, &mx_inv, &mtmp);
        Mat_Mul_4X4(&mtmp, &my_inv, &mrot);
    }
    break;

    default:
        break;
    } // end switch

    // now mrot holds the concatenated product of inverse rotation matrices
    // multiply the inverse translation matrix against it and store in the
    // camera objects' camera transform matrix we are done!
    Mat_Mul_4X4(&mt_inv, &mrot, &cam->mcam);

} // end Build_CAM4DV1_Matrix_Euler

void World_To_Camera_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,
                                    CAM4DV1_PTR cam)
{
    // NOTE: this is a matrix based function
    // this function transforms each polygon in the global render list
    // to camera coordinates based on the sent camera transform matrix
    // you would use this function instead of the object based function
    // if you decided earlier in the pipeline to turn each object into
    // a list of polygons and then add them to the global render list
    // the conversion of an object into polygons probably would have
    // happened after object culling, local transforms, local to world
    // and backface culling, so the minimum number of polygons from
    // each object are in the list, note that the function assumes
    // that at LEAST the local to world transform has been called
    // and the polygon data is in the transformed list tvlist of
    // the POLYF4DV1 object

    // transform each polygon in the render list into camera coordinates
    // assumes the render list has already been transformed to world
    // coordinates and the result is in tvlist[] of each polygon object

    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire current polygon
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // transform this polygon if and only if it's not clipped, not culled,
        // active, and visible, note however the concept of "backface" is
        // irrelevant in a wire frame engine though
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV1_STATE_BACKFACE))
            continue; // move onto next poly

        // all good, let's transform
        for (int vertex = 0; vertex < 3; vertex++)
        {
            // transform the vertex by the mcam matrix within the camera
            // it better be valid!
            POINT4D presult; // hold result of each transformation

            // transform point
            Mat_Mul_VECTOR4D_4X4(&curr_poly->tvlist[vertex], &cam->mcam, &presult);

            // store result back
            VECTOR4D_COPY(&curr_poly->tvlist[vertex], &presult);
        } // end for vertex

    } // end for poly

} // end World_To_Camera_RENDERLIST4DV1

void Camera_To_Perspective_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,
                                          CAM4DV1_PTR cam)
{
    // NOTE: this is not a matrix based function
    // this function transforms each polygon in the global render list
    // into perspective coordinates, based on the
    // sent camera object,
    // you would use this function instead of the object based function
    // if you decided earlier in the pipeline to turn each object into
    // a list of polygons and then add them to the global render list

    // transform each polygon in the render list into camera coordinates
    // assumes the render list has already been transformed to world
    // coordinates and the result is in tvlist[] of each polygon object

    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire current polygon
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // transform this polygon if and only if it's not clipped, not culled,
        // active, and visible, note however the concept of "backface" is
        // irrelevant in a wire frame engine though
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV1_STATE_BACKFACE))
            continue; // move onto next poly

        // all good, let's transform
        for (int vertex = 0; vertex < 3; vertex++)
        {
            float z = curr_poly->tvlist[vertex].z;

            // transform the vertex by the view parameters in the camera
            curr_poly->tvlist[vertex].x = cam->view_dist * curr_poly->tvlist[vertex].x / z;
            curr_poly->tvlist[vertex].y = cam->view_dist * curr_poly->tvlist[vertex].y * cam->aspect_ratio / z;
            // z = z, so no change

            // not that we are NOT dividing by the homogenous w coordinate since
            // we are not using a matrix operation for this version of the function

        } // end for vertex

    } // end for poly

} // end Camera_To_Perspective_RENDERLIST4DV1

void Perspective_To_Screen_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, CAM4DV1_PTR cam)
{
    // NOTE: this is not a matrix based function
    // this function transforms the perspective coordinates of the render
    // list into screen coordinates, based on the sent viewport in the camera
    // assuming that the viewplane coordinates were normalized
    // you would use this function instead of the object based function
    // if you decided earlier in the pipeline to turn each object into
    // a list of polygons and then add them to the global render list
    // you would only call this function if you previously performed
    // a normalized perspective transform

    // transform each polygon in the render list from perspective to screen
    // coordinates assumes the render list has already been transformed
    // to normalized perspective coordinates and the result is in tvlist[]
    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire current polygon
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // transform this polygon if and only if it's not clipped, not culled,
        // active, and visible, note however the concept of "backface" is
        // irrelevant in a wire frame engine though
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV1_STATE_BACKFACE))
            continue; // move onto next poly

        float alpha = (0.5 * cam->viewport_width - 0.5);
        float beta = (0.5 * cam->viewport_height - 0.5);

        // all good, let's transform
        for (int vertex = 0; vertex < 3; vertex++)
        {
            // the vertex is in perspective normalized coords from -1 to 1
            // on each axis, simple scale them and invert y axis and project
            // to screen

            // transform the vertex by the view parameters in the camera
            curr_poly->tvlist[vertex].x = alpha + alpha * curr_poly->tvlist[vertex].x;
            curr_poly->tvlist[vertex].y = beta - beta * curr_poly->tvlist[vertex].y;
        } // end for vertex

    } // end for poly

} // end Perspective_To_Screen_RENDERLIST4DV1

void Mat_Mul_VECTOR4D_4X4(VECTOR4D_PTR va,
                          MATRIX4X4_PTR mb,
                          VECTOR4D_PTR vprod)
{
    // this function multiplies a VECTOR4D against a
    // 4x4 matrix - ma*mb and stores the result in mprod
    // the function makes no assumptions

    for (int col = 0; col < 4; col++)
    {
        // compute dot product from row of ma
        // and column of mb
        float sum = 0; // used to hold result

        for (int row = 0; row < 4; row++)
        {
            // add in next product pair
            sum += (va->M[row] * mb->M[row][col]);
        } // end for index

        // insert resulting col element
        vprod->M[col] = sum;

    } // end for col

} // end Mat_Mul_VECTOR4D_4X4

void VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vsum)
{
    // this function adds va+vb and return it in vsum
    vsum->x = va->x + vb->x;
    vsum->y = va->y + vb->y;
    vsum->z = va->z + vb->z;
    vsum->w = 1;

} // end VECTOR4D_Add

////////////////////////////////////////////////////////////

VECTOR4D VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb)
{
    // this function adds va+vb and returns the result on
    // the stack
    VECTOR4D vsum;

    vsum.x = va->x + vb->x;
    vsum.y = va->y + vb->y;
    vsum.z = va->z + vb->z;
    vsum.w = 1;

    // return result
    return (vsum);

} // end VECTOR4D_Add

//transform
void Reset_OBJECT4DV1(OBJECT4DV1_PTR obj)
{
    RESET_BIT(obj->state, OBJECT4DV1_STATE_CULLED);

    for (int poly = 0; poly < obj->num_polys; poly++)
    {
        POLY4DV1_PTR curr_poly = &obj->plist[poly];

        if (!(curr_poly->state & POLY4DV1_STATE_ACTIVE))
            continue;
        RESET_BIT(curr_poly->state, POLY4DV1_STATE_CLIPPED);
        RESET_BIT(curr_poly->state, POLY4DV1_STATE_BACKFACE);
    }
}

//通用变换函数
void Transform_OBJECT4DV1(OBJECT4DV1_PTR obj, MATRIX4X4_PTR mt, int coord_select, int transform_basis)
{
    switch (coord_select)
    {
    case TRANSFORM_LOCAL_ONLY: //local到local
    {
        for (int vertex = 0; vertex < obj->num_vertices; vertex++)
        {
            POINT4D presult;
            Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex], mt, &presult);
            VECTOR4D_COPY(&obj->vlist_local[vertex], &presult);
        }
    }
    break;

    case TRANSFORM_TRANS_ONLY: //trans到trans
    {
        for (int vertex = 0; vertex < obj->num_vertices; vertex++)
        {
            POINT4D presult;
            Mat_Mul_VECTOR4D_4X4(&obj->vlist_trans[vertex], mt, &presult);
            VECTOR4D_COPY(&obj->vlist_trans[vertex], &presult);
        }
    }
    break;

    case TRANSFORM_LOCAL_TO_TRANS: //local到trans
    {
        for (int vertex = 0; vertex < obj->num_vertices; vertex++)
        {
            POINT4D presult;
            Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex], mt, &obj->vlist_trans[vertex]);
        }
    }
    break;

    default:
        break;
    }

    if (transform_basis)
    {
        VECTOR4D vresult;

        // rotate ux of basis
        Mat_Mul_VECTOR4D_4X4(&obj->ux, mt, &vresult);
        VECTOR4D_COPY(&obj->ux, &vresult);

        // rotate uy of basis
        Mat_Mul_VECTOR4D_4X4(&obj->uy, mt, &vresult);
        VECTOR4D_COPY(&obj->uy, &vresult);

        // rotate uz of basis
        Mat_Mul_VECTOR4D_4X4(&obj->uz, mt, &vresult);
        VECTOR4D_COPY(&obj->uz, &vresult);
    }
}

//世界坐标的变换，结果保存在obj->vlist_trans
void Model_To_World_OBJECT4DV1(OBJECT4DV1_PTR obj, int coord_select)
{
    if (coord_select == TRANSFORM_LOCAL_TO_TRANS)
    {
        for (int vertex = 0; vertex < obj->num_vertices; vertex++)
        {
            VECTOR4D_Add(&obj->vlist_local[vertex], &obj->world_pos, &obj->vlist_trans[vertex]);
        }
    }
    else
    {
        for (int vertex = 0; vertex < obj->num_vertices; vertex++)
        {
            VECTOR4D_Add(&obj->vlist_trans[vertex], &obj->world_pos, &obj->vlist_trans[vertex]);
        }
    }
}

//对所有obj->vlist_trans，执行相机变换，结果保存在obj->vlist_trans
void World_To_Camera_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam)
{
    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        POINT4D presult;
        Mat_Mul_VECTOR4D_4X4(&obj->vlist_trans[vertex], &cam->mcam, &presult);
        VECTOR4D_COPY(&obj->vlist_trans[vertex], &presult);
    }
}

//对所有obj->vlist_trans，执行透视变换，结果保存在obj->vlist_trans
void Camera_To_Perspective_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam)
{
    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        float z = obj->vlist_trans[vertex].z;
        obj->vlist_trans[vertex].x = cam->view_dist * obj->vlist_trans[vertex].x / z;
        obj->vlist_trans[vertex].y = cam->view_dist * obj->vlist_trans[vertex].y * cam->aspect_ratio / z;
    }
}

//对所有obj->vlist_trans，执行屏幕变换，结果保存在obj->vlist_trans
void Perspective_To_Screen_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam)
{

    float alpha = (0.5 * cam->viewport_width - 0.5);
    float beta = (0.5 * cam->viewport_height - 0.5);

    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        obj->vlist_trans[vertex].x = alpha + alpha * obj->vlist_trans[vertex].x;
        obj->vlist_trans[vertex].y = beta - beta * obj->vlist_trans[vertex].y;
    }
}

void Remove_Backfaces_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam)
{
    if (obj->state & OBJECT4DV1_STATE_CULLED)
        return;

    for (int poly = 0; poly < obj->num_polys; poly++)
    {
        POLY4DV1_PTR curr_poly = &obj->plist[poly];

        if (!(curr_poly->state & POLY4DV1_STATE_ACTIVE) || (curr_poly->state & POLY4DV1_STATE_CLIPPED) || (curr_poly->attr & POLY4DV1_ATTR_2SIDED) || (curr_poly->state & POLY4DV1_STATE_BACKFACE))
            continue;

        int vindex_0 = curr_poly->vert[0];
        int vindex_1 = curr_poly->vert[1];
        int vindex_2 = curr_poly->vert[2];

        VECTOR4D u, v, n;

        // 计算向量u和v
        VECTOR4D_Build(&obj->vlist_trans[vindex_0], &obj->vlist_trans[vindex_1], &u);
        VECTOR4D_Build(&obj->vlist_trans[vindex_0], &obj->vlist_trans[vindex_2], &v);

        // 计算法向量
        VECTOR4D_Cross(&u, &v, &n);

        //计算观察向量
        VECTOR4D view;
        VECTOR4D_Build(&obj->vlist_trans[vindex_0], &cam->pos, &view);

        float dp = VECTOR4D_Dot(&n, &view);

        if (dp <= 0.0)
            SET_BIT(curr_poly->state, POLY4DV1_STATE_BACKFACE);
    }
}

void Remove_Backfaces_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, CAM4DV1_PTR cam)
{

    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {

        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
            (curr_poly->attr & POLY4DV1_ATTR_2SIDED) ||
            (curr_poly->state & POLY4DV1_STATE_BACKFACE))
            continue;

        VECTOR4D u, v, n;

        VECTOR4D_Build(&curr_poly->tvlist[0], &curr_poly->tvlist[1], &u);
        VECTOR4D_Build(&curr_poly->tvlist[0], &curr_poly->tvlist[2], &v);

        VECTOR4D_Cross(&u, &v, &n);

        VECTOR4D view;
        VECTOR4D_Build(&curr_poly->tvlist[0], &cam->pos, &view);

        float dp = VECTOR4D_Dot(&n, &view);

        if (dp <= 0.0)
            SET_BIT(curr_poly->state, POLY4DV1_STATE_BACKFACE);
    }
}

int Cull_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam, int cull_flags)
{

    POINT4D sphere_pos;
    Mat_Mul_VECTOR4D_4X4(&obj->world_pos, &cam->mcam, &sphere_pos); //将包围球转换到相机空间

    if (cull_flags & CULL_OBJECT_Z_PLANE)
    {

        if (((sphere_pos.z - obj->max_radius) > cam->far_clip_z) || ((sphere_pos.z + obj->max_radius) < cam->near_clip_z))
        {
            SET_BIT(obj->state, OBJECT4DV1_STATE_CULLED);
            return (1);
        }
    }

    if (cull_flags & CULL_OBJECT_X_PLANE)
    {

        float z_test = (0.5) * cam->viewplane_width * sphere_pos.z / cam->view_dist;

        if (((sphere_pos.x - obj->max_radius) > z_test) || ((sphere_pos.x + obj->max_radius) < -z_test))
        {
            SET_BIT(obj->state, OBJECT4DV1_STATE_CULLED);
            return (1);
        }
    }

    if (cull_flags & CULL_OBJECT_Y_PLANE)
    {

        float z_test = (0.5) * cam->viewplane_height * sphere_pos.z / cam->view_dist;

        if (((sphere_pos.y - obj->max_radius) > z_test) || ((sphere_pos.y + obj->max_radius) < -z_test))
        {
            SET_BIT(obj->state, OBJECT4DV1_STATE_CULLED);
            return (1);
        }
    }

    return (0);
}

int Insert_POLY4DV1_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, POLY4DV1_PTR poly)
{

    if (rend_list->num_polys >= RENDERLIST4DV1_MAX_POLYS)
        return (0);

    rend_list->poly_ptrs[rend_list->num_polys] = &rend_list->poly_data[rend_list->num_polys];

    rend_list->poly_data[rend_list->num_polys].state = poly->state;
    rend_list->poly_data[rend_list->num_polys].attr = poly->attr;
    rend_list->poly_data[rend_list->num_polys].color = poly->color;

    VECTOR4D_COPY(&rend_list->poly_data[rend_list->num_polys].tvlist[0],
                  &poly->vlist[poly->vert[0]]);

    VECTOR4D_COPY(&rend_list->poly_data[rend_list->num_polys].tvlist[1],
                  &poly->vlist[poly->vert[1]]);

    VECTOR4D_COPY(&rend_list->poly_data[rend_list->num_polys].tvlist[2],
                  &poly->vlist[poly->vert[2]]);

    VECTOR4D_COPY(&rend_list->poly_data[rend_list->num_polys].vlist[0],
                  &poly->vlist[poly->vert[0]]);

    VECTOR4D_COPY(&rend_list->poly_data[rend_list->num_polys].vlist[1],
                  &poly->vlist[poly->vert[1]]);

    VECTOR4D_COPY(&rend_list->poly_data[rend_list->num_polys].vlist[2],
                  &poly->vlist[poly->vert[2]]);

    if (rend_list->num_polys == 0)
    {

        rend_list->poly_data[0].next = NULL;
        rend_list->poly_data[0].prev = NULL;
    }
    else
    {

        rend_list->poly_data[rend_list->num_polys].next = NULL;
        rend_list->poly_data[rend_list->num_polys].prev =
            &rend_list->poly_data[rend_list->num_polys - 1];

        rend_list->poly_data[rend_list->num_polys - 1].next =
            &rend_list->poly_data[rend_list->num_polys];
    }

    rend_list->num_polys++;

    return (1);
}

int Insert_OBJECT4DV1_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, OBJECT4DV1_PTR obj, int insert_local)
{

    if (!(obj->state & OBJECT4DV1_STATE_ACTIVE) ||
        (obj->state & OBJECT4DV1_STATE_CULLED) ||
        !(obj->state & OBJECT4DV1_STATE_VISIBLE))
        return (0);

    for (int poly = 0; poly < obj->num_polys; poly++)
    {

        POLY4DV1_PTR curr_poly = &obj->plist[poly];

        if (!(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV1_STATE_BACKFACE))
            continue;

        POINT4D_PTR vlist_old = curr_poly->vlist;

        if (insert_local)
            curr_poly->vlist = obj->vlist_local;
        else
            curr_poly->vlist = obj->vlist_trans;

        if (!Insert_POLY4DV1_RENDERLIST4DV1(rend_list, curr_poly))
        {

            curr_poly->vlist = vlist_old;

            return (0);
        }

        curr_poly->vlist = vlist_old;
    }

    return (1);
}

//参数1指向参数2的向量，保存到参数3
void VECTOR4D_Build(VECTOR4D_PTR init, VECTOR4D_PTR term, VECTOR4D_PTR result)
{
    // build a 4d vector
    result->x = term->x - init->x;
    result->y = term->y - init->y;
    result->z = term->z - init->z;
    result->w = 1;

} // end VECTOR4D_Buil

//向量叉乘
void VECTOR4D_Cross(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vn)
{
    vn->x = ((va->y * vb->z) - (va->z * vb->y));
    vn->y = -((va->x * vb->z) - (va->z * vb->x));
    vn->z = ((va->x * vb->y) - (va->y * vb->x));
    vn->w = 1;
}

//向量叉乘
VECTOR4D VECTOR4D_Cross(VECTOR4D_PTR va, VECTOR4D_PTR vb)
{
    VECTOR4D vn;
    vn.x = ((va->y * vb->z) - (va->z * vb->y));
    vn.y = -((va->x * vb->z) - (va->z * vb->x));
    vn.z = ((va->x * vb->y) - (va->y * vb->x));
    vn.w = 1;
    return (vn);
}

//向量点乘
float VECTOR4D_Dot(VECTOR4D_PTR va, VECTOR4D_PTR vb)
{
    return ((va->x * vb->x) + (va->y * vb->y) + (va->z * vb->z));
}


int Init_Light_LIGHTV1(int index,          // index of light to create (0..MAX_LIGHTS-1)
                       int _state,         // state of light
                       int _attr,          // type of light, and extra qualifiers
                       RGBAV1 _c_ambient,  // ambient light intensity
                       RGBAV1 _c_diffuse,  // diffuse light intensity
                       RGBAV1 _c_specular, // specular light intensity
                       POINT4D_PTR _pos,   // position of light
                       VECTOR4D_PTR _dir,  // direction of light
                       float _kc,          // attenuation factors
                       float _kl,
                       float _kq,
                       float _spot_inner, // inner angle for spot light
                       float _spot_outer, // outer angle for spot light
                       float _pf)         // power factor/falloff for spot lights
{
    // this function initializes a light based on the flags sent in _attr, values that
    // aren't needed are set to 0 by caller

    // make sure light is in range
    if (index < 0 || index >= MAX_LIGHTS)
        return (0);

    // all good, initialize the light (many fields may be dead)
    lights[index].state = _state; // state of light
    lights[index].id = index;     // id of light
    lights[index].attr = _attr;   // type of light, and extra qualifiers

    lights[index].c_ambient = _c_ambient;   // ambient light intensity
    lights[index].c_diffuse = _c_diffuse;   // diffuse light intensity
    lights[index].c_specular = _c_specular; // specular light intensity

    lights[index].kc = _kc; // constant, linear, and quadratic attenuation factors
    lights[index].kl = _kl;
    lights[index].kq = _kq;

    if (_pos)
        VECTOR4D_COPY(&lights[index].pos, _pos); // position of light

    if (_dir)
    {
        VECTOR4D_COPY(&lights[index].dir, _dir); // direction of light
        // normalize it
        VECTOR4D_Normalize(&lights[index].dir);

    } // end if

    lights[index].spot_inner = _spot_inner; // inner angle for spot light
    lights[index].spot_outer = _spot_outer; // outer angle for spot light
    lights[index].pf = _pf;                 // power factor/falloff for spot lights

    // return light index as success
    return (index);

} // end Create_Light_LIGHTV1

int Reset_Lights_LIGHTV1(void)
{
    // this function simply resets all lights in the system
    static int first_time = 1;

    memset(lights, 0, MAX_LIGHTS * sizeof(LIGHTV1));

    // reset number of lights
    num_lights = 0;

    // reset first time
    first_time = 0;

    // return success
    return (1);

} // end Reset_Lights_LIGHTV1

void Reset_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list)
{
    rend_list->num_polys = 0; // that was hard!
} // end Reset_RENDERLIST4DV2

void Reset_OBJECT4DV2(OBJECT4DV2_PTR obj)
{
    // this function resets the sent object and redies it for
    // transformations, basically just resets the culled, clipped and
    // backface flags, but here's where you would add stuff
    // to ready any object for the pipeline
    // the object is valid, let's rip it apart polygon by polygon
    // note: works on the entire object, all frames

    // reset object's culled flag
    RESET_BIT(obj->state, OBJECT4DV2_STATE_CULLED);

    // now the clipped and backface flags for the polygons
    for (int poly = 0; poly < obj->num_polys; poly++)
    {
        // acquire polygon
        POLY4DV2_PTR curr_poly = &obj->plist[poly];

        // first is this polygon even visible?
        if (!(curr_poly->state & POLY4DV2_STATE_ACTIVE))
            continue; // move onto next poly

        // reset clipped and backface flags
        RESET_BIT(curr_poly->state, POLY4DV2_STATE_CLIPPED);
        RESET_BIT(curr_poly->state, POLY4DV2_STATE_BACKFACE);
        RESET_BIT(curr_poly->state, POLY4DV2_STATE_LIT);

    } // end for poly

} // end Reset_OBJECT4DV2

void Transform_OBJECT4DV2(OBJECT4DV2_PTR obj,  // object to transform
                          MATRIX4X4_PTR mt,    // transformation matrix
                          int coord_select,    // selects coords to transform
                          int transform_basis, // flags if vector orientation
                                               // should be transformed too
                          int all_frames)      // should all frames be transformed

{
    // this function simply transforms all of the vertices in the local or trans
    // array by the sent matrix, since the object may have multiple frames, it
    // takes that into consideration
    // also vertex normals are rotated, however, if there is a translation factor
    // in the sent matrix that will corrupt the normals, later we might want to
    // null out the last row of the matrix before transforming the normals?
    // future optimization: set flag in object attributes, and objects without
    // vertex normals can be rotated without the test in line

    // single frame or all frames?
    if (!all_frames)
    {
        // what coordinates should be transformed?
        switch (coord_select)
        {
        case TRANSFORM_LOCAL_ONLY:
        {
            // transform each local/model vertex of the object mesh in place
            for (int vertex = 0; vertex < obj->num_vertices; vertex++)
            {
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex].v, mt, &presult);

                // store result back
                VECTOR4D_COPY(&obj->vlist_local[vertex].v, &presult);

                // transform vertex normal if needed
                if (obj->vlist_local[vertex].attr & VERTEX4DTV1_ATTR_NORMAL)
                {
                    // transform normal
                    Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex].n, mt, &presult);

                    // store result back
                    VECTOR4D_COPY(&obj->vlist_local[vertex].n, &presult);
                } // end if

            } // end for index
        }
        break;

        case TRANSFORM_TRANS_ONLY:
        {
            // transform each "transformed" vertex of the object mesh in place
            // remember, the idea of the vlist_trans[] array is to accumulate
            // transformations
            for (int vertex = 0; vertex < obj->num_vertices; vertex++)
            {
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&obj->vlist_trans[vertex].v, mt, &presult);

                // store result back
                VECTOR4D_COPY(&obj->vlist_trans[vertex].v, &presult);

                // transform vertex normal if needed
                if (obj->vlist_trans[vertex].attr & VERTEX4DTV1_ATTR_NORMAL)
                {
                    // transform normal
                    Mat_Mul_VECTOR4D_4X4(&obj->vlist_trans[vertex].n, mt, &presult);

                    // store result back
                    VECTOR4D_COPY(&obj->vlist_trans[vertex].n, &presult);
                } // end if

            } // end for index
        }
        break;

        case TRANSFORM_LOCAL_TO_TRANS:
        {
            // transform each local/model vertex of the object mesh and store result
            // in "transformed" vertex list
            for (int vertex = 0; vertex < obj->num_vertices; vertex++)
            {
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex].v, mt, &obj->vlist_trans[vertex].v);

                // transform vertex normal if needed
                if (obj->vlist_local[vertex].attr & VERTEX4DTV1_ATTR_NORMAL)
                {
                    // transform point
                    Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex].n, mt, &obj->vlist_trans[vertex].n);
                } // end if

            } // end for index
        }
        break;

        default:
            break;

        } // end switch

    }    // end if single frame
    else // transform all frames
    {
        // what coordinates should be transformed?
        switch (coord_select)
        {
        case TRANSFORM_LOCAL_ONLY:
        {
            // transform each local/model vertex of the object mesh in place
            for (int vertex = 0; vertex < obj->total_vertices; vertex++)
            {
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&obj->head_vlist_local[vertex].v, mt, &presult);

                // store result back
                VECTOR4D_COPY(&obj->head_vlist_local[vertex].v, &presult);

                // transform vertex normal if needed
                if (obj->head_vlist_local[vertex].attr & VERTEX4DTV1_ATTR_NORMAL)
                {
                    // transform normal
                    Mat_Mul_VECTOR4D_4X4(&obj->head_vlist_local[vertex].n, mt, &presult);

                    // store result back
                    VECTOR4D_COPY(&obj->head_vlist_local[vertex].n, &presult);
                } // end if

            } // end for index
        }
        break;

        case TRANSFORM_TRANS_ONLY:
        {
            // transform each "transformed" vertex of the object mesh in place
            // remember, the idea of the vlist_trans[] array is to accumulate
            // transformations
            for (int vertex = 0; vertex < obj->total_vertices; vertex++)
            {
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&obj->head_vlist_trans[vertex].v, mt, &presult);

                // store result back
                VECTOR4D_COPY(&obj->head_vlist_trans[vertex].v, &presult);

                // transform vertex normal if needed
                if (obj->head_vlist_trans[vertex].attr & VERTEX4DTV1_ATTR_NORMAL)
                {
                    // transform normal
                    Mat_Mul_VECTOR4D_4X4(&obj->head_vlist_trans[vertex].n, mt, &presult);

                    // store result back
                    VECTOR4D_COPY(&obj->head_vlist_trans[vertex].n, &presult);
                } // end if

            } // end for index
        }
        break;

        case TRANSFORM_LOCAL_TO_TRANS:
        {
            // transform each local/model vertex of the object mesh and store result
            // in "transformed" vertex list
            for (int vertex = 0; vertex < obj->total_vertices; vertex++)
            {
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&obj->head_vlist_local[vertex].v, mt, &obj->head_vlist_trans[vertex].v);

                // transform vertex normal if needed
                if (obj->head_vlist_local[vertex].attr & VERTEX4DTV1_ATTR_NORMAL)
                {
                    // transform point
                    Mat_Mul_VECTOR4D_4X4(&obj->head_vlist_local[vertex].n, mt, &obj->head_vlist_trans[vertex].n);
                } // end if

            } // end for index
        }
        break;

        default:
            break;

        } // end switch

    } // end else multiple frames

    // finally, test if transform should be applied to orientation basis
    // hopefully this is a rotation, otherwise the basis will get corrupted
    if (transform_basis)
    {
        // now rotate orientation basis for object
        VECTOR4D vresult; // use to rotate each orientation vector axis

        // rotate ux of basis
        Mat_Mul_VECTOR4D_4X4(&obj->ux, mt, &vresult);
        VECTOR4D_COPY(&obj->ux, &vresult);

        // rotate uy of basis
        Mat_Mul_VECTOR4D_4X4(&obj->uy, mt, &vresult);
        VECTOR4D_COPY(&obj->uy, &vresult);

        // rotate uz of basis
        Mat_Mul_VECTOR4D_4X4(&obj->uz, mt, &vresult);
        VECTOR4D_COPY(&obj->uz, &vresult);
    } // end if

} // end Transform_OBJECT4DV2

void Model_To_World_OBJECT4DV2(OBJECT4DV2_PTR obj,
                               int coord_select,
                               int all_frames)
{
    // NOTE: Not matrix based
    // this function converts the local model coordinates of the
    // sent object into world coordinates, the results are stored
    // in the transformed vertex list (vlist_trans) within the object

    // interate thru vertex list and transform all the model/local
    // coords to world coords by translating the vertex list by
    // the amount world_pos and storing the results in vlist_trans[]
    // no need to transform vertex normals, they are invariant of position

    if (!all_frames)
    {
        if (coord_select == TRANSFORM_LOCAL_TO_TRANS)
        {
            for (int vertex = 0; vertex < obj->num_vertices; vertex++)
            {
                // translate vertex
                VECTOR4D_Add(&obj->vlist_local[vertex].v, &obj->world_pos, &obj->vlist_trans[vertex].v);
                // copy normal
                VECTOR4D_COPY(&obj->vlist_trans[vertex].n, &obj->vlist_local[vertex].n);

            } // end for vertex
        }     // end if local
        else
        { // TRANSFORM_TRANS_ONLY
            for (int vertex = 0; vertex < obj->num_vertices; vertex++)
            {
                // std::cout<<"TRANSFORM_TRANS_ONLY"<<std::endl;
                // translate vertex
                VECTOR4D_Add(&obj->vlist_trans[vertex].v, &obj->world_pos, &obj->vlist_trans[vertex].v);
            } // end for vertex
        }     // end else trans

    }    // end if single frame
    else // all frames
    {
        if (coord_select == TRANSFORM_LOCAL_TO_TRANS)
        {
            for (int vertex = 0; vertex < obj->total_vertices; vertex++)
            {
                // translate vertex
                VECTOR4D_Add(&obj->head_vlist_local[vertex].v, &obj->world_pos, &obj->head_vlist_trans[vertex].v);
                // copy normal
                VECTOR4D_COPY(&obj->head_vlist_trans[vertex].n, &obj->head_vlist_local[vertex].n);
            } // end for vertex
        }     // end if local
        else
        { // TRANSFORM_TRANS_ONLY
            for (int vertex = 0; vertex < obj->total_vertices; vertex++)
            {
                // translate vertex
                VECTOR4D_Add(&obj->head_vlist_trans[vertex].v, &obj->world_pos, &obj->head_vlist_trans[vertex].v);
            } // end for vertex
        }     // end else trans

    } // end if all frames

} // end Model_To_World_OBJECT4DV2

int Insert_OBJECT4DV2_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
                                     OBJECT4DV2_PTR obj,
                                     int insert_local = 0)

{
    // { andre work in progress, rewrite with materials...}

    // converts the entire object into a face list and then inserts
    // the visible, active, non-clipped, non-culled polygons into
    // the render list, also note the flag insert_local control
    // whether or not the vlist_local or vlist_trans vertex list
    // is used, thus you can insert an object "raw" totally untranformed
    // if you set insert_local to 1, default is 0, that is you would
    // only insert an object after at least the local to world transform
    // the last parameter is used to control if their has been
    // a lighting step that has generated a light value stored
    // in the upper 16-bits of color, if lighting_on = 1 then
    // this value is used to overwrite the base color of the
    // polygon when its sent to the rendering list

    unsigned int base_color; // save base color of polygon

    // is this objective inactive or culled or invisible?
    if (!(obj->state & OBJECT4DV2_STATE_ACTIVE) ||
        (obj->state & OBJECT4DV2_STATE_CULLED) ||
        !(obj->state & OBJECT4DV2_STATE_VISIBLE))
        return (0);

    // the object is valid, let's rip it apart polygon by polygon
    for (int poly = 0; poly < obj->num_polys; poly++)
    {
        // acquire polygon
        POLY4DV2_PTR curr_poly = &obj->plist[poly];
        // std::cout << curr_poly->color << std::endl;

        // first is this polygon even visible?
        if (!(curr_poly->state & POLY4DV2_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV2_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV2_STATE_BACKFACE))
            continue; // move onto next poly

        // override vertex list polygon refers to
        // the case that you want the local coords used
        // first save old pointer

        VERTEX4DTV1_PTR vlist_old = curr_poly->vlist;

        if (insert_local)
            curr_poly->vlist = obj->vlist_local;
        else
            curr_poly->vlist = obj->vlist_trans;

        // 丢到了 obj->vlist_trans[vertex].v
        // now insert this polygon
        if (!Insert_POLY4DV2_RENDERLIST4DV2(rend_list, curr_poly))
        {
            // fix vertex list pointer
            curr_poly->vlist = vlist_old;

            // the whole object didn't fit!
            return (0);
        } // end if

        // fix vertex list pointer
        curr_poly->vlist = vlist_old;

    } // end for

    // return success
    return (1);

} // end Insert_OBJECT4DV2_RENDERLIST4DV2

void Remove_Backfaces_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list, CAM4DV1_PTR cam)
{
    // NOTE: this is not a matrix based function
    // this function removes the backfaces from polygon list
    // the function does this based on the polygon list data
    // tvlist along with the camera position (only)
    // note that only the backface state is set in each polygon

    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire current polygon
        POLYF4DV2_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // test this polygon if and only if it's not clipped, not culled,
        // active, and visible and not 2 sided. Note we test for backface in the event that
        // a previous call might have already determined this, so why work
        // harder!
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV2_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV2_STATE_CLIPPED) ||
            (curr_poly->attr & POLY4DV2_ATTR_2SIDED) ||
            (curr_poly->state & POLY4DV2_STATE_BACKFACE))
            continue; // move onto next poly

        // we need to compute the normal of this polygon face, and recall
        // that the vertices are in cw order, u = p0->p1, v=p0->p2, n=uxv
        VECTOR4D u, v, n;

        // build u, v
        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[1].v, &u);
        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[2].v, &v);

        // compute cross product
        VECTOR4D_Cross(&u, &v, &n);

        // now create eye vector to viewpoint
        VECTOR4D view;
        VECTOR4D_Build(&curr_poly->tvlist[0].v, &cam->pos, &view);

        // and finally, compute the dot product
        float dp = VECTOR4D_Dot(&n, &view);

        // if the sign is > 0 then visible, 0 = scathing, < 0 invisible
        if (dp <= 0.0)
            SET_BIT(curr_poly->state, POLY4DV2_STATE_BACKFACE);

    } // end for poly

} // end Remove_Backfaces_RENDERLIST4DV2

void World_To_Camera_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
                                    CAM4DV1_PTR cam)
{
    // NOTE: this is a matrix based function
    // this function transforms each polygon in the global render list
    // to camera coordinates based on the sent camera transform matrix
    // you would use this function instead of the object based function
    // if you decided earlier in the pipeline to turn each object into
    // a list of polygons and then add them to the global render list
    // the conversion of an object into polygons probably would have
    // happened after object culling, local transforms, local to world
    // and backface culling, so the minimum number of polygons from
    // each object are in the list, note that the function assumes
    // that at LEAST the local to world transform has been called
    // and the polygon data is in the transformed list tvlist of
    // the POLYF4DV1 object

    // transform each polygon in the render list into camera coordinates
    // assumes the render list has already been transformed to world
    // coordinates and the result is in tvlist[] of each polygon object

    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire current polygon
        POLYF4DV2_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // transform this polygon if and only if it's not clipped, not culled,
        // active, and visible, note however the concept of "backface" is
        // irrelevant in a wire frame engine though
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV2_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV2_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV2_STATE_BACKFACE))
            continue; // move onto next poly

        // all good, let's transform
        for (int vertex = 0; vertex < 3; vertex++)
        {
            // transform the vertex by the mcam matrix within the camera
            // it better be valid!
            POINT4D presult; // hold result of each transformation

            // transform point
            Mat_Mul_VECTOR4D_4X4(&curr_poly->tvlist[vertex].v, &cam->mcam, &presult);

            // store result back
            VECTOR4D_COPY(&curr_poly->tvlist[vertex].v, &presult);
        } // end for vertex

    } // end for poly

} // end World_To_Camera_RENDERLIST4DV2

void Sort_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list, int sort_method = SORT_POLYLIST_AVGZ)
{
    // this function sorts the rendering list based on the polygon z-values
    // the specific sorting method is controlled by sending in control flags
    // #define SORT_POLYLIST_AVGZ  0 - sorts on average of all vertices
    // #define SORT_POLYLIST_NEARZ 1 - sorts on closest z vertex of each poly
    // #define SORT_POLYLIST_FARZ  2 - sorts on farthest z vertex of each poly

    switch (sort_method)
    {
    case SORT_POLYLIST_AVGZ: //  - sorts on average of all vertices
    {
        qsort((void *)rend_list->poly_ptrs, rend_list->num_polys, sizeof(POLYF4DV2_PTR), Compare_AvgZ_POLYF4DV2);
    }
    break;

    case SORT_POLYLIST_NEARZ: // - sorts on closest z vertex of each poly
    {
        qsort((void *)rend_list->poly_ptrs, rend_list->num_polys, sizeof(POLYF4DV2_PTR), Compare_NearZ_POLYF4DV2);
    }
    break;

    case SORT_POLYLIST_FARZ: //  - sorts on farthest z vertex of each poly
    {
        qsort((void *)rend_list->poly_ptrs, rend_list->num_polys, sizeof(POLYF4DV2_PTR), Compare_FarZ_POLYF4DV2);
    }
    break;

    default:
        break;
    } // end switch

} // end Sort_RENDERLIST4DV2

void Camera_To_Perspective_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
                                          CAM4DV1_PTR cam)
{
    // NOTE: this is not a matrix based function
    // this function transforms each polygon in the global render list
    // into perspective coordinates, based on the
    // sent camera object,
    // you would use this function instead of the object based function
    // if you decided earlier in the pipeline to turn each object into
    // a list of polygons and then add them to the global render list

    // transform each polygon in the render list into camera coordinates
    // assumes the render list has already been transformed to world
    // coordinates and the result is in tvlist[] of each polygon object

    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire current polygon
        POLYF4DV2_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // transform this polygon if and only if it's not clipped, not culled,
        // active, and visible, note however the concept of "backface" is
        // irrelevant in a wire frame engine though
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV2_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV2_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV2_STATE_BACKFACE))
            continue; // move onto next poly

        // all good, let's transform
        for (int vertex = 0; vertex < 3; vertex++)
        {
            float z = curr_poly->tvlist[vertex].z;

            // transform the vertex by the view parameters in the camera
            curr_poly->tvlist[vertex].x = cam->view_dist * curr_poly->tvlist[vertex].x / z;
            curr_poly->tvlist[vertex].y = cam->view_dist * curr_poly->tvlist[vertex].y * cam->aspect_ratio / z;
            // z = z, so no change

            // not that we are NOT dividing by the homogenous w coordinate since
            // we are not using a matrix operation for this version of the function

        } // end for vertex

    } // end for poly

} // end Camera_To_Perspective_RENDERLIST4DV2

void Perspective_To_Screen_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
                                          CAM4DV1_PTR cam)
{
    // NOTE: this is not a matrix based function
    // this function transforms the perspective coordinates of the render
    // list into screen coordinates, based on the sent viewport in the camera
    // assuming that the viewplane coordinates were normalized
    // you would use this function instead of the object based function
    // if you decided earlier in the pipeline to turn each object into
    // a list of polygons and then add them to the global render list
    // you would only call this function if you previously performed
    // a normalized perspective transform

    // transform each polygon in the render list from perspective to screen
    // coordinates assumes the render list has already been transformed
    // to normalized perspective coordinates and the result is in tvlist[]
    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire current polygon
        POLYF4DV2_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // transform this polygon if and only if it's not clipped, not culled,
        // active, and visible, note however the concept of "backface" is
        // irrelevant in a wire frame engine though
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV2_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV2_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV2_STATE_BACKFACE))
            continue; // move onto next poly

        float alpha = (0.5 * cam->viewport_width - 0.5);
        float beta = (0.5 * cam->viewport_height - 0.5);

        // all good, let's transform
        for (int vertex = 0; vertex < 3; vertex++)
        {
            // the vertex is in perspective normalized coords from -1 to 1
            // on each axis, simple scale them and invert y axis and project
            // to screen

            // transform the vertex by the view parameters in the camera
            curr_poly->tvlist[vertex].x = alpha + alpha * curr_poly->tvlist[vertex].x;
            curr_poly->tvlist[vertex].y = beta - beta * curr_poly->tvlist[vertex].y;
        } // end for vertex

    } // end for poly

} // end Perspective_To_Screen_RENDERLIST4DV2

int Init_OBJECT4DV2(OBJECT4DV2_PTR obj, // object to allocate
                    int _num_vertices,
                    int _num_polys,
                    int _num_frames,
                    int destroy)
{
    // this function does nothing more than allocate the memory for an OBJECT4DV2
    // based on the sent data, later we may want to create more robust initializers
    // but the problem is that we don't want to tie the initializer to anthing yet
    // in 99% of cases this all will be done by the call to load object
    // we just might need this function if we manually want to build an object???

    // first destroy the object if it exists
    if (destroy)
        Destroy_OBJECT4DV2(obj);

    // allocate memory for vertex lists
    if (!(obj->vlist_local = (VERTEX4DTV1_PTR)malloc(sizeof(VERTEX4DTV1) * _num_vertices * _num_frames)))
        return (0);

    // clear data
    memset((void *)obj->vlist_local, 0, sizeof(VERTEX4DTV1) * _num_vertices * _num_frames);

    if (!(obj->vlist_trans = (VERTEX4DTV1_PTR)malloc(sizeof(VERTEX4DTV1) * _num_vertices * _num_frames)))
        return (0);

    // clear data
    memset((void *)obj->vlist_trans, 0, sizeof(VERTEX4DTV1) * _num_vertices * _num_frames);

    // number of texture coordinates always 3*number of polys
    if (!(obj->tlist = (POINT2D_PTR)malloc(sizeof(POINT2D) * _num_polys * 3)))
        return (0);

    // clear data
    memset((void *)obj->tlist, 0, sizeof(POINT2D) * _num_polys * 3);

    // allocate memory for radii arrays
    if (!(obj->avg_radius = (float *)malloc(sizeof(float) * _num_frames)))
        return (0);

    // clear data
    memset((void *)obj->avg_radius, 0, sizeof(float) * _num_frames);

    if (!(obj->max_radius = (float *)malloc(sizeof(float) * _num_frames)))
        return (0);

    // clear data
    memset((void *)obj->max_radius, 0, sizeof(float) * _num_frames);

    // allocate memory for polygon list
    if (!(obj->plist = (POLY4DV2_PTR)malloc(sizeof(POLY4DV2) * _num_polys)))
        return (0);

    // clear data
    memset((void *)obj->plist, 0, sizeof(POLY4DV2) * _num_polys);

    // alias head pointers
    obj->head_vlist_local = obj->vlist_local;
    obj->head_vlist_trans = obj->vlist_trans;

    // set some internal variables
    obj->num_frames = _num_frames;
    obj->num_polys = _num_polys;
    obj->num_vertices = _num_vertices;
    obj->total_vertices = _num_vertices * _num_frames;

    // return success
    return (1);

} // end Init_OBJECT4DV2


float Compute_OBJECT4DV2_Radius(OBJECT4DV2_PTR obj)
{
    // this function computes the average and maximum radius for
    // sent object and opdates the object data for the "current frame"
    // it's up to the caller to make sure Set_Frame() for this object
    // has been called to set the object up properly

    // reset incase there's any residue
    obj->avg_radius[obj->curr_frame] = 0;
    obj->max_radius[obj->curr_frame] = 0;

    // loop thru and compute radius
    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        // update the average and maximum radius
        float dist_to_vertex =
            sqrt(obj->vlist_local[vertex].x * obj->vlist_local[vertex].x +
                 obj->vlist_local[vertex].y * obj->vlist_local[vertex].y +
                 obj->vlist_local[vertex].z * obj->vlist_local[vertex].z);

        // accumulate total radius
        obj->avg_radius[obj->curr_frame] += dist_to_vertex;

        // update maximum radius
        if (dist_to_vertex > obj->max_radius[obj->curr_frame])
            obj->max_radius[obj->curr_frame] = dist_to_vertex;

    } // end for vertex

    // finallize average radius computation
    obj->avg_radius[obj->curr_frame] /= obj->num_vertices;

    // return max radius of frame 0
    return (obj->max_radius[0]);

} // end Compute_OBJECT4DV2_Radius

int Compute_OBJECT4DV2_Poly_Normals(OBJECT4DV2_PTR obj)
{
    // the normal of a polygon is commonly needed in a number
    // of functions, however, to store a normal turns out to
    // be counterproductive in most cases since the transformation
    // to rotate the normal ends up taking as long as computing the
    // normal -- HOWEVER, if the normal must have unit length, then
    // pre-computing the length of the normal, and then in real-time
    // dividing by this save a length computation, so we get the
    // best of both worlds... thus, this function computes the length
    // of a polygon's normal, but care must be taken, so that we compute
    // the length based on the EXACT same two vectors that all other
    // functions will use when computing the normal
    // in most cases the functions of interest are the lighting functions
    // if we can pre-compute the normal length
    // for all these functions then that will save at least:
    // num_polys_per_frame * (time to compute length of vector)

    // the way we have written the engine, in all cases the normals
    // during lighting are computed as u = v0->v1, and v = v1->v2
    // so as long as we follow that convention we are fine.
    // also, since the new OBJECT4DV2 format supports multiple frames
    // we must perform these calculations for EACH frame of the animation
    // since although the poly indices don't change, the vertice positions
    // do and thus, so do the normals!!!

    // is this object valid
    if (!obj)
        return (0);

    // iterate thru the poly list of the object and compute normals
    // each polygon
    for (int poly = 0; poly < obj->num_polys; poly++)
    {

        // extract vertex indices into master list, rember the polygons are
        // NOT self contained, but based on the vertex list stored in the object
        // itself
        int vindex_0 = obj->plist[poly].vert[0];
        int vindex_1 = obj->plist[poly].vert[1];
        int vindex_2 = obj->plist[poly].vert[2];

        // we need to compute the normal of this polygon face, and recall
        // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv
        VECTOR4D u, v, n;

        // build u, v
        VECTOR4D_Build(&obj->vlist_local[vindex_0].v, &obj->vlist_local[vindex_1].v, &u);
        VECTOR4D_Build(&obj->vlist_local[vindex_0].v, &obj->vlist_local[vindex_2].v, &v);

        // compute cross product
        VECTOR4D_Cross(&u, &v, &n);

        // compute length of normal accurately and store in poly nlength
        // +- epsilon later to fix over/underflows
        obj->plist[poly].nlength = VECTOR4D_Length(&n);
    } // end for poly

    // return success
    return (1);

} // end Compute_OBJECT4DV2_Poly_Normals

int Compute_OBJECT4DV2_Vertex_Normals(OBJECT4DV2_PTR obj)
{
    // the vertex normals of each polygon are commonly needed in a number
    // functions, most importantly lighting calculations for gouraud shading
    // however, we only need to compute the vertex normals for polygons that are
    // gouraud shader, so for every vertex we must determine the polygons that
    // share the vertex then compute the average normal, to determine if a polygon
    // contributes we look at the shading flags for the polygon

    // is this object valid
    if (!obj)
        return (0);

    // algorithm: we are going to scan the polygon list and for every polygon
    // that needs normals we are going to "accumulate" the surface normal into all
    // vertices that the polygon touches, and increment a counter to track how many
    // polys contribute to vertex, then when the scan is done the counts will be used
    // to average the accumulated values, so instead of an O(n^2) algorithm, we get a O(c*n)

    // this tracks the polygon indices that touch a particular vertex
    // the array is used to count the number of contributors to the vertex
    // so at the end of the process we can divide each "accumulated" normal
    // and average
    int polys_touch_vertex[OBJECT4DV2_MAX_VERTICES];
    memset((void *)polys_touch_vertex, 0, sizeof(int) * OBJECT4DV2_MAX_VERTICES);

    // iterate thru the poly list of the object, compute its normal, then add
    // each vertice that composes it to the "touching" vertex array
    // while accumulating the normal in the vertex normal array

    for (int poly = 0; poly < obj->num_polys; poly++)
    {
        //Write_Error("\nprocessing poly %d", poly);

        // test if this polygon needs vertex normals
        if (obj->plist[poly].attr & POLY4DV2_ATTR_SHADE_MODE_GOURAUD)
        {
            // extract vertex indices into master list, rember the polygons are
            // NOT self contained, but based on the vertex list stored in the object
            // itself
            int vindex_0 = obj->plist[poly].vert[0];
            int vindex_1 = obj->plist[poly].vert[1];
            int vindex_2 = obj->plist[poly].vert[2];

            //Write_Error("\nTouches vertices: %d, %d, %d", vindex_0, vindex_1, vindex_2);

            // we need to compute the normal of this polygon face, and recall
            // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv
            VECTOR4D u, v, n;

            // build u, v
            VECTOR4D_Build(&obj->vlist_local[vindex_0].v, &obj->vlist_local[vindex_1].v, &u);
            VECTOR4D_Build(&obj->vlist_local[vindex_0].v, &obj->vlist_local[vindex_2].v, &v);

            // compute cross product
            VECTOR4D_Cross(&u, &v, &n);

            // update vertex array to flag this polygon as a contributor
            polys_touch_vertex[vindex_0]++;
            polys_touch_vertex[vindex_1]++;
            polys_touch_vertex[vindex_2]++;

            //Write_Error("\nPoly touch array v[%d] = %d,  v[%d] = %d,  v[%d] = %d", vindex_0, polys_touch_vertex[vindex_0],                                                                           vindex_1, polys_touch_vertex[vindex_1],                                                                           vindex_2, polys_touch_vertex[vindex_2]);

            // now accumulate the normal into the vertex normal itself
            // note, we do NOT normalize at this point since we want the length of the normal
            // to weight on the average, and since the length is in fact the area of the parallelogram
            // constructed by uxv, so we are taking the "influence" of the area into consideration
            VECTOR4D_Add(&obj->vlist_local[vindex_0].n, &n, &obj->vlist_local[vindex_0].n);
            VECTOR4D_Add(&obj->vlist_local[vindex_1].n, &n, &obj->vlist_local[vindex_1].n);
            VECTOR4D_Add(&obj->vlist_local[vindex_2].n, &n, &obj->vlist_local[vindex_2].n);
        } // end for poly

    } // end if needs vertex normals

    // now we are almost done, we have accumulated all the vertex normals, but need to average them
    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        // if this vertex has any contributors then it must need averaging, OR we could check
        // the shading hints flags, they should be one to one
        //Write_Error("\nProcessing vertex: %d, attr: %d, contributors: %d", vertex,                                                                        obj->vlist_local[vertex].attr,                                                                        polys_touch_vertex[vertex]);

        // test if this vertex has a normal and needs averaging
        if (polys_touch_vertex[vertex] >= 1)
        {
            obj->vlist_local[vertex].nx /= polys_touch_vertex[vertex];
            obj->vlist_local[vertex].ny /= polys_touch_vertex[vertex];
            obj->vlist_local[vertex].nz /= polys_touch_vertex[vertex];

            // now normalize the normal
            VECTOR4D_Normalize(&obj->vlist_local[vertex].n);

            //Write_Error("\nAvg Vertex normal: [%f, %f, %f]", obj->vlist_local[vertex].nx,                                                        obj->vlist_local[vertex].ny,                                                        obj->vlist_local[vertex].nz);

        } // end if

    } // end for

    // return success
    return (1);

} // end Compute_OBJECT4DV2_Vertex_Normals

void VECTOR4D_Normalize(VECTOR4D_PTR va)
{
    // normalizes the sent vector and returns the result

    // compute length
    float length = sqrtf(va->x * va->x + va->y * va->y + va->z * va->z);

    // test for zero length vector
    // if found return zero vector
    if (length < EPSILON_E5)
        return;

    float length_inv = 1.0 / length;

    // compute normalized version of vector
    va->x *= length_inv;
    va->y *= length_inv;
    va->z *= length_inv;
    va->w = 1;

} // end VECTOR4D_Normalize

float VECTOR4D_Length_Fast(VECTOR4D_PTR va)
{
    // computes the magnitude of a vector using an approximation
    // very fast
    return (Fast_Distance_3D(va->x, va->y, va->z));

} // end VECTOR4D_Length_Fast

float Fast_Distance_3D(float fx, float fy, float fz)
{
    // this function computes the distance from the origin to x,y,z

    int temp;    // used for swaping
    int x, y, z; // used for algorithm

    // make sure values are all positive
    x = fabs(fx) * 1024;
    y = fabs(fy) * 1024;
    z = fabs(fz) * 1024;

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

} // end Fast_Distance_3D

int Insert_POLY4DV2_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
                                   POLY4DV2_PTR poly)
{
    // converts the sent POLY4DV2 into a POLYF4DV2 and inserts it
    // into the render list, this function needs optmizing

    // step 0: are we full?
    if (rend_list->num_polys >= RENDERLIST4DV2_MAX_POLYS)
        return (0);

    // step 1: copy polygon into next opening in polygon render list

    // point pointer to polygon structure
    rend_list->poly_ptrs[rend_list->num_polys] = &rend_list->poly_data[rend_list->num_polys];

    // copy fields { ??????????? make sure ALL fields are copied, normals, textures, etc!!!  }
    rend_list->poly_data[rend_list->num_polys].state = poly->state;
    rend_list->poly_data[rend_list->num_polys].attr = poly->attr;
    rend_list->poly_data[rend_list->num_polys].color = poly->color;
    rend_list->poly_data[rend_list->num_polys].nlength = poly->nlength;
    rend_list->poly_data[rend_list->num_polys].texture = poly->texture;

    // poly could be lit, so copy these too...
    rend_list->poly_data[rend_list->num_polys].lit_color[0] = poly->lit_color[0];
    rend_list->poly_data[rend_list->num_polys].lit_color[1] = poly->lit_color[1];
    rend_list->poly_data[rend_list->num_polys].lit_color[2] = poly->lit_color[2];

    // now copy vertices, be careful! later put a loop, but for now
    // know there are 3 vertices always!
    VERTEX4DTV1_COPY(&rend_list->poly_data[rend_list->num_polys].tvlist[0],
                     &poly->vlist[poly->vert[0]]);

    VERTEX4DTV1_COPY(&rend_list->poly_data[rend_list->num_polys].tvlist[1],
                     &poly->vlist[poly->vert[1]]);

    VERTEX4DTV1_COPY(&rend_list->poly_data[rend_list->num_polys].tvlist[2],
                     &poly->vlist[poly->vert[2]]);

    // and copy into local vertices too
    VERTEX4DTV1_COPY(&rend_list->poly_data[rend_list->num_polys].vlist[0],
                     &poly->vlist[poly->vert[0]]);

    VERTEX4DTV1_COPY(&rend_list->poly_data[rend_list->num_polys].vlist[1],
                     &poly->vlist[poly->vert[1]]);

    VERTEX4DTV1_COPY(&rend_list->poly_data[rend_list->num_polys].vlist[2],
                     &poly->vlist[poly->vert[2]]);

    // finally the texture coordinates, this has to be performed manually
    // since at this point in the pipeline the vertices do NOT have texture
    // coordinate, the polygons DO, however, now, there are 3 vertices for
    // EVERY polygon, rather than vertex sharing, so we can copy the texture
    // coordinates out of the indexed arrays into the VERTEX4DTV1 structures
    rend_list->poly_data[rend_list->num_polys].tvlist[0].t = poly->tlist[poly->text[0]];
    rend_list->poly_data[rend_list->num_polys].tvlist[1].t = poly->tlist[poly->text[1]];
    rend_list->poly_data[rend_list->num_polys].tvlist[2].t = poly->tlist[poly->text[2]];

    rend_list->poly_data[rend_list->num_polys].vlist[0].t = poly->tlist[poly->text[0]];
    rend_list->poly_data[rend_list->num_polys].vlist[1].t = poly->tlist[poly->text[1]];
    rend_list->poly_data[rend_list->num_polys].vlist[2].t = poly->tlist[poly->text[2]];

    // now the polygon is loaded into the next free array position, but
    // we need to fix up the links

    // test if this is the first entry
    if (rend_list->num_polys == 0)
    {
        // set pointers to null, could loop them around though to self
        rend_list->poly_data[0].next = NULL;
        rend_list->poly_data[0].prev = NULL;
    } // end if
    else
    {
        // first set this node to point to previous node and next node (null)
        rend_list->poly_data[rend_list->num_polys].next = NULL;
        rend_list->poly_data[rend_list->num_polys].prev =
            &rend_list->poly_data[rend_list->num_polys - 1];

        // now set previous node to point to this node
        rend_list->poly_data[rend_list->num_polys - 1].next =
            &rend_list->poly_data[rend_list->num_polys];
    } // end else

    // increment number of polys in list
    rend_list->num_polys++;

    // return successful insertion
    return (1);

} // end Insert_POLY4DV2_RENDERLIST4DV2

int Compare_AvgZ_POLYF4DV2(const void *arg1, const void *arg2)
{
    // this function comapares the average z's of two polygons and is used by the
    // depth sort surface ordering algorithm

    float z1, z2;

    POLYF4DV2_PTR poly_1, poly_2;

    // dereference the poly pointers
    poly_1 = *((POLYF4DV2_PTR *)(arg1));
    poly_2 = *((POLYF4DV2_PTR *)(arg2));

    // compute z average of each polygon
    z1 = (float)0.33333 * (poly_1->tvlist[0].z + poly_1->tvlist[1].z + poly_1->tvlist[2].z);

    // now polygon 2
    z2 = (float)0.33333 * (poly_2->tvlist[0].z + poly_2->tvlist[1].z + poly_2->tvlist[2].z);

    // compare z1 and z2, such that polys' will be sorted in descending Z order
    if (z1 > z2)
        return (-1);
    else if (z1 < z2)
        return (1);
    else
        return (0);

} // end Compare_AvgZ_POLYF4DV2

////////////////////////////////////////////////////////////////////////////////

int Compare_NearZ_POLYF4DV2(const void *arg1, const void *arg2)
{
    // this function comapares the closest z's of two polygons and is used by the
    // depth sort surface ordering algorithm

    float z1, z2;

    POLYF4DV2_PTR poly_1, poly_2;

    // dereference the poly pointers
    poly_1 = *((POLYF4DV2_PTR *)(arg1));
    poly_2 = *((POLYF4DV2_PTR *)(arg2));

    // compute the near z of each polygon
    z1 = MIN(poly_1->tvlist[0].z, poly_1->tvlist[1].z);
    z1 = MIN(z1, poly_1->tvlist[2].z);

    z2 = MIN(poly_2->tvlist[0].z, poly_2->tvlist[1].z);
    z2 = MIN(z2, poly_2->tvlist[2].z);

    // compare z1 and z2, such that polys' will be sorted in descending Z order
    if (z1 > z2)
        return (-1);
    else if (z1 < z2)
        return (1);
    else
        return (0);

} // end Compare_NearZ_POLYF4DV2

////////////////////////////////////////////////////////////////////////////////
int Destroy_OBJECT4DV2(OBJECT4DV2_PTR obj) // object to destroy
{
    // this function destroys the sent object, basically frees the memory
    // if any that has been allocated

    // local vertex list
    if (obj->head_vlist_local)
        free(obj->head_vlist_local);

    // transformed vertex list
    if (obj->head_vlist_trans)
        free(obj->head_vlist_trans);

    // texture coordinate list
    if (obj->tlist)
        free(obj->tlist);

    // polygon list
    if (obj->plist)
        free(obj->plist);

    // object radii arrays
    if (obj->avg_radius)
        free(obj->avg_radius);

    if (obj->max_radius)
        free(obj->max_radius);

    // now clear out object completely
    memset((void *)obj, 0, sizeof(OBJECT4DV2));

    // return success
    return (1);

} // end Destroy_OBJECT4DV

int Compare_FarZ_POLYF4DV2(const void *arg1, const void *arg2)
{
    // this function comapares the farthest z's of two polygons and is used by the
    // depth sort surface ordering algorithm

    float z1, z2;

    POLYF4DV2_PTR poly_1, poly_2;

    // dereference the poly pointers
    poly_1 = *((POLYF4DV2_PTR *)(arg1));
    poly_2 = *((POLYF4DV2_PTR *)(arg2));

    // compute the near z of each polygon
    z1 = MAX(poly_1->tvlist[0].z, poly_1->tvlist[1].z);
    z1 = MAX(z1, poly_1->tvlist[2].z);

    z2 = MAX(poly_2->tvlist[0].z, poly_2->tvlist[1].z);
    z2 = MAX(z2, poly_2->tvlist[2].z);

    // compare z1 and z2, such that polys' will be sorted in descending Z order
    if (z1 > z2)
        return (-1);
    else if (z1 < z2)
        return (1);
    else
        return (0);

} // end Compare_FarZ_POLYF4DV2

float VECTOR4D_Length(VECTOR4D_PTR va)
{
    // computes the magnitude of a vector, slow

    return (sqrtf(va->x * va->x + va->y * va->y + va->z * va->z));

} // end VECTOR4D_Length

//计算出当前poly的color
int Light_RENDERLIST4DV2_World16(RENDERLIST4DV2_PTR rend_list, // list to process
                                 CAM4DV1_PTR cam,              // camera position
                                 LIGHTV1_PTR lights,           // light list (might have more than one)
                                 int max_lights)               // maximum lights in list
{

    // 16-bit version of function
    // function lights the entire rendering list based on the sent lights and camera. the function supports
    // constant/pure shading (emmisive), flat shading with ambient, infinite, point lights, and spot lights
    // note that this lighting function is rather brute force and simply follows the math, however
    // there are some clever integer operations that are used in scale 256 rather than going to floating
    // point, but why? floating point and ints are the same speed, HOWEVER, the conversion to and from floating
    // point can be cycle intensive, so if you can keep your calcs in ints then you can gain some speed
    // also note, type 1 spot lights are simply point lights with direction, the "cone" is more of a function
    // of the falloff due to attenuation, but they still look like spot lights
    // type 2 spot lights are implemented with the intensity having a dot product relationship with the
    // angle from the surface point to the light direction just like in the optimized model, but the pf term
    // that is used for a concentration control must be 1,2,3,.... integral and non-fractional
    // this function now performs emissive, flat, and gouraud lighting, results are stored in the
    // lit_color[] array of each polygon

    unsigned int r_base, g_base, b_base, // base color being lit
        r_sum, g_sum, b_sum,             // sum of lighting process over all lights
        r_sum0, g_sum0, b_sum0,
        r_sum1, g_sum1, b_sum1,
        r_sum2, g_sum2, b_sum2,
        ri, gi, bi,
        shaded_color; // final color

    float dp, // dot product
        dist, // distance from light to surface
        dists,
        i,     // general intensities
        nl,    // length of normal
        atten; // attenuation computations

    VECTOR4D u, v, n, l, d, s; // used for cross product and light vector calculations

    //Write_Error("\nEntering lighting function");

    // for each valid poly, light it...
    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire polygon
        POLYF4DV2_PTR curr_poly = rend_list->poly_ptrs[poly];

        // light this polygon if and only if it's not clipped, not culled,
        // active, and visible
        if (!(curr_poly->state & POLY4DV2_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV2_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV2_STATE_BACKFACE) ||
            (curr_poly->state & POLY4DV2_STATE_LIT))
            continue; // move onto next poly

        //Write_Error("\npoly %d",poly);

        // set state of polygon to lit
        SET_BIT(curr_poly->state, POLY4DV2_STATE_LIT);

        // we will use the transformed polygon vertex list since the backface removal
        // only makes sense at the world coord stage further of the pipeline

        // test the lighting mode of the polygon (use flat for flat, gouraud))
        if (curr_poly->attr & POLY4DV2_ATTR_SHADE_MODE_FLAT)
        {
            //Write_Error("\nEntering Flat Shader");

            // step 1: extract the base color out in RGB mode
            // assume 565 format
            _RGB565FROM16BIT(curr_poly->color, &r_base, &g_base, &b_base);
            // scale to 8 bit
            r_base <<= 3;
            g_base <<= 2;
            b_base <<= 3;

            // std::cout<<"color = " << r_base<<" "<< g_base <<" "<< b_base<<std::endl;
            //Write_Error("\nBase color=%d,%d,%d", r_base, g_base, b_base);

            // initialize color sum
            r_sum = 0;
            g_sum = 0;
            b_sum = 0;

            //Write_Error("\nsum color=%d,%d,%d", r_sum, g_sum, b_sum);

            // new optimization:
            // when there are multiple lights in the system we will end up performing numerous
            // redundant calculations to minimize this my strategy is to set key variables to
            // to MAX values on each loop, then during the lighting calcs to test the vars for
            // the max value, if they are the max value then the first light that needs the math
            // will do it, and then save the information into the variable (causing it to change state
            // from an invalid number) then any other lights that need the math can use the previously
            // computed value

            // set surface normal.z to FLT_MAX to flag it as non-computed
            n.z = FLT_MAX;

            // loop thru lights
            for (int curr_light = 0; curr_light < max_lights; curr_light++)
            {
                // is this light active
                if (lights[curr_light].state == LIGHTV1_STATE_OFF)
                    continue;

                //Write_Error("\nprocessing light %d",curr_light);

                // what kind of light are we dealing with
                if (lights[curr_light].attr & LIGHTV1_ATTR_AMBIENT)
                {
                    //Write_Error("\nEntering ambient light...");

                    // simply multiply each channel against the color of the
                    // polygon then divide by 256 to scale back to 0..255
                    // use a shift in real life!!! >> 8
                    r_sum += ((lights[curr_light].c_ambient.r * r_base) / 256);
                    g_sum += ((lights[curr_light].c_ambient.g * g_base) / 256);
                    b_sum += ((lights[curr_light].c_ambient.b * b_base) / 256);

                    //Write_Error("\nambient sum=%d,%d,%d", r_sum, g_sum, b_sum);

                    // there better only be one ambient light!

                }                                                         // end if
                else if (lights[curr_light].attr & LIGHTV1_ATTR_INFINITE) ///////////////////////////////////////////
                {
                    //Write_Error("\nEntering infinite light...");

                    // infinite lighting, we need the surface normal, and the direction
                    // of the light source

                    // test if we already computed poly normal in previous calculation
                    if (n.z == FLT_MAX)
                    {
                        // we need to compute the normal of this polygon face, and recall
                        // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv

                        // build u, v
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[1].v, &u);
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[2].v, &v);

                        // compute cross product
                        VECTOR4D_Cross(&u, &v, &n);
                    } // end if

                    // at this point, we are almost ready, but we have to normalize the normal vector!
                    // this is a key optimization we can make later, we can pre-compute the length of all polygon
                    // normals, so this step can be optimized
                    // compute length of normal
                    //nl = VECTOR4D_Length_Fast2(&n);
                    nl = curr_poly->nlength;

                    // ok, recalling the lighting model for infinite lights
                    // I(d)dir = I0dir * Cldir
                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles

                    dp = VECTOR4D_Dot(&n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        i = 128 * dp / nl;
                        r_sum += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                    //Write_Error("\ninfinite sum=%d,%d,%d", r_sum, g_sum, b_sum);

                }                                                      // end if infinite light
                else if (lights[curr_light].attr & LIGHTV1_ATTR_POINT) ///////////////////////////////////////
                {
                    //Write_Error("\nEntering point light...");

                    // perform point light computations
                    // light model for point light is once again:
                    //              I0point * Clpoint
                    //  I(d)point = ___________________
                    //              kc +  kl*d + kq*d2
                    //
                    //  Where d = |p - s|
                    // thus it's almost identical to the infinite light, but attenuates as a function
                    // of distance from the point source to the surface point being lit

                    // test if we already computed poly normal in previous calculation
                    if (n.z == FLT_MAX)
                    {
                        // we need to compute the normal of this polygon face, and recall
                        // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv

                        // build u, v
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[1].v, &u);
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[2].v, &v);

                        // compute cross product
                        VECTOR4D_Cross(&u, &v, &n);
                    } // end if

                    // at this point, we are almost ready, but we have to normalize the normal vector!
                    // this is a key optimization we can make later, we can pre-compute the length of all polygon
                    // normals, so this step can be optimized
                    // compute length of normal
                    //nl = VECTOR4D_Length_Fast2(&n);
                    nl = curr_poly->nlength;

                    // compute vector from surface to light
                    VECTOR4D_Build(&curr_poly->tvlist[0].v, &lights[curr_light].pos, &l);

                    // compute distance and attenuation
                    dist = VECTOR4D_Length_Fast2(&l);

                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles
                    dp = VECTOR4D_Dot(&n, &l);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (nl * dist * atten);

                        r_sum += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                    //Write_Error("\npoint sum=%d,%d,%d",r_sum,g_sum,b_sum);

                }                                                           // end if point
                else if (lights[curr_light].attr & LIGHTV1_ATTR_SPOTLIGHT1) ////////////////////////////////////
                {
                    //Write_Error("\nentering spot light1...");

                    // perform spotlight/point computations simplified model that uses
                    // point light WITH a direction to simulate a spotlight
                    // light model for point light is once again:
                    //              I0point * Clpoint
                    //  I(d)point = ___________________
                    //              kc +  kl*d + kq*d2
                    //
                    //  Where d = |p - s|
                    // thus it's almost identical to the infinite light, but attenuates as a function
                    // of distance from the point source to the surface point being lit

                    // test if we already computed poly normal in previous calculation
                    if (n.z == FLT_MAX)
                    {
                        // we need to compute the normal of this polygon face, and recall
                        // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv

                        // build u, v
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[1].v, &u);
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[2].v, &v);

                        // compute cross product
                        VECTOR4D_Cross(&u, &v, &n);
                    } // end if

                    // at this point, we are almost ready, but we have to normalize the normal vector!
                    // this is a key optimization we can make later, we can pre-compute the length of all polygon
                    // normals, so this step can be optimized
                    // compute length of normal
                    //nl = VECTOR4D_Length_Fast2(&n);
                    nl = curr_poly->nlength;

                    // compute vector from surface to light
                    VECTOR4D_Build(&curr_poly->tvlist[0].v, &lights[curr_light].pos, &l);

                    // compute distance and attenuation
                    dist = VECTOR4D_Length_Fast2(&l);

                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles

                    // note that I use the direction of the light here rather than a the vector to the light
                    // thus we are taking orientation into account which is similar to the spotlight model
                    dp = VECTOR4D_Dot(&n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (nl * atten);

                        r_sum += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                    //Write_Error("\nspotlight sum=%d,%d,%d",r_sum, g_sum, b_sum);

                }                                                           // end if spotlight1
                else if (lights[curr_light].attr & LIGHTV1_ATTR_SPOTLIGHT2) // simple version ////////////////////
                {
                    //Write_Error("\nEntering spotlight2 ...");

                    // perform spot light computations
                    // light model for spot light simple version is once again:
                    //         	     I0spotlight * Clspotlight * MAX( (l . s), 0)^pf
                    // I(d)spotlight = __________________________________________
                    //               		 kc + kl*d + kq*d2
                    // Where d = |p - s|, and pf = power factor

                    // thus it's almost identical to the point, but has the extra term in the numerator
                    // relating the angle between the light source and the point on the surface

                    // test if we already computed poly normal in previous calculation
                    if (n.z == FLT_MAX)
                    {
                        // we need to compute the normal of this polygon face, and recall
                        // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv

                        // build u, v
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[1].v, &u);
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[2].v, &v);

                        // compute cross product
                        VECTOR4D_Cross(&u, &v, &n);
                    } // end if

                    // at this point, we are almost ready, but we have to normalize the normal vector!
                    // this is a key optimization we can make later, we can pre-compute the length of all polygon
                    // normals, so this step can be optimized
                    // compute length of normal
                    //nl = VECTOR4D_Length_Fast2(&n);
                    nl = curr_poly->nlength;

                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles
                    dp = VECTOR4D_Dot(&n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        // compute vector from light to surface (different from l which IS the light dir)
                        VECTOR4D_Build(&lights[curr_light].pos, &curr_poly->tvlist[0].v, &s);

                        // compute length of s (distance to light source) to normalize s for lighting calc
                        dists = VECTOR4D_Length_Fast2(&s);

                        // compute spot light term (s . l)
                        float dpsl = VECTOR4D_Dot(&s, &lights[curr_light].dir) / dists;

                        // proceed only if term is positive
                        if (dpsl > 0)
                        {
                            // compute attenuation
                            atten = (lights[curr_light].kc + lights[curr_light].kl * dists + lights[curr_light].kq * dists * dists);

                            // for speed reasons, pf exponents that are less that 1.0 are out of the question, and exponents
                            // must be integral
                            float dpsl_exp = dpsl;

                            // exponentiate for positive integral powers
                            for (int e_index = 1; e_index < (int)lights[curr_light].pf; e_index++)
                                dpsl_exp *= dpsl;

                            // now dpsl_exp holds (dpsl)^pf power which is of course (s . l)^pf

                            i = 128 * dp * dpsl_exp / (nl * atten);

                            r_sum += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                            g_sum += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                            b_sum += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);

                        } // end if

                    } // end if

                    //Write_Error("\nSpotlight sum=%d,%d,%d",r_sum, g_sum, b_sum);

                } // end if spot light

            } // end for light

            // make sure colors aren't out of range
            if (r_sum > 255)
                r_sum = 255;
            if (g_sum > 255)
                g_sum = 255;
            if (b_sum > 255)
                b_sum = 255;

            //Write_Error("\nWriting final values to polygon %d = %d,%d,%d", poly, r_sum, g_sum, b_sum);

            // write the color over current color
            curr_poly->lit_color[0] = RGB16Bit565(r_sum, g_sum, b_sum);

        }                                                       // end if
        if (curr_poly->attr & POLY4DV2_ATTR_SHADE_MODE_GOURAUD) /////////////////////////////////
        {
            // gouraud shade, unfortunetly at this point in the pipeline, we have lost the original
            // mesh, and only have triangles, thus, many triangles will share the same vertices and
            // they will get lit 2x since we don't have any way to tell this, alas, performing lighting
            // at the object level is a better idea when gouraud shading is performed since the
            // commonality of vertices is still intact, in any case, lighting here is similar to polygon
            // flat shaded, but we do it 3 times, once for each vertex, additionally there are lots
            // of opportunities for optimization, but I am going to lay off them for now, so the code
            // is intelligible, later we will optimize

            //Write_Error("\nEntering gouraud shader...");

            // step 1: extract the base color out in RGB mode
            // assume 565 format
            _RGB565FROM16BIT(curr_poly->color, &r_base, &g_base, &b_base);

            // scale to 8 bit
            r_base <<= 3;
            g_base <<= 2;
            b_base <<= 3;

            //Write_Error("\nBase color=%d, %d, %d", r_base, g_base, b_base);

            // initialize color sum(s) for vertices
            r_sum0 = 0;
            g_sum0 = 0;
            b_sum0 = 0;

            r_sum1 = 0;
            g_sum1 = 0;
            b_sum1 = 0;

            r_sum2 = 0;
            g_sum2 = 0;
            b_sum2 = 0;

            //Write_Error("\nColor sums rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,   r_sum1, g_sum1, b_sum1, r_sum2, g_sum2, b_sum2);

            // new optimization:
            // when there are multiple lights in the system we will end up performing numerous
            // redundant calculations to minimize this my strategy is to set key variables to
            // to MAX values on each loop, then during the lighting calcs to test the vars for
            // the max value, if they are the max value then the first light that needs the math
            // will do it, and then save the information into the variable (causing it to change state
            // from an invalid number) then any other lights that need the math can use the previously
            // computed value, however, since we already have the normals, not much here to cache on
            // a large scale, but small scale stuff is there, however, we will optimize those later

            // loop thru lights
            for (int curr_light = 0; curr_light < max_lights; curr_light++)
            {
                // is this light active
                if (lights[curr_light].state == LIGHTV1_STATE_OFF)
                    continue;

                //Write_Error("\nprocessing light %d", curr_light);

                // what kind of light are we dealing with
                if (lights[curr_light].attr & LIGHTV1_ATTR_AMBIENT) ///////////////////////////////
                {
                    //Write_Error("\nEntering ambient light....");

                    // simply multiply each channel against the color of the
                    // polygon then divide by 256 to scale back to 0..255
                    // use a shift in real life!!! >> 8
                    ri = ((lights[curr_light].c_ambient.r * r_base) / 256);
                    gi = ((lights[curr_light].c_ambient.g * g_base) / 256);
                    bi = ((lights[curr_light].c_ambient.b * b_base) / 256);

                    // ambient light has the same affect on each vertex
                    r_sum0 += ri;
                    g_sum0 += gi;
                    b_sum0 += bi;

                    r_sum1 += ri;
                    g_sum1 += gi;
                    b_sum1 += bi;

                    r_sum2 += ri;
                    g_sum2 += gi;
                    b_sum2 += bi;

                    // there better only be one ambient light!
                    //Write_Error("\nexiting ambient ,sums rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,  r_sum1, g_sum1, b_sum1,  r_sum2, g_sum2, b_sum2);

                }                                                         // end if
                else if (lights[curr_light].attr & LIGHTV1_ATTR_INFINITE) /////////////////////////////////
                {
                    //Write_Error("\nentering infinite light...");

                    // infinite lighting, we need the surface normal, and the direction
                    // of the light source

                    // no longer need to compute normal or length, we already have the vertex normal
                    // and it's length is 1.0
                    // ....

                    // ok, recalling the lighting model for infinite lights
                    // I(d)dir = I0dir * Cldir
                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles

                    // need to perform lighting for each vertex (lots of redundant math, optimize later!)

                    //Write_Error("\nv0=[%f, %f, %f]=%f, v1=[%f, %f, %f]=%f, v2=[%f, %f, %f]=%f",
                    // curr_poly->tvlist[0].n.x, curr_poly->tvlist[0].n.y,curr_poly->tvlist[0].n.z, VECTOR4D_Length(&curr_poly->tvlist[0].n),
                    // curr_poly->tvlist[1].n.x, curr_poly->tvlist[1].n.y,curr_poly->tvlist[1].n.z, VECTOR4D_Length(&curr_poly->tvlist[1].n),
                    // curr_poly->tvlist[2].n.x, curr_poly->tvlist[2].n.y,curr_poly->tvlist[2].n.z, VECTOR4D_Length(&curr_poly->tvlist[2].n) );

                    // vertex 0
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[0].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        i = 128 * dp;
                        r_sum0 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum0 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum0 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                    // vertex 1
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[1].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        i = 128 * dp;
                        r_sum1 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum1 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum1 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                    // vertex 2
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[2].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        i = 128 * dp;
                        r_sum2 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum2 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum2 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                    //Write_Error("\nexiting infinite, color sums rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,  r_sum1, g_sum1, b_sum1,  r_sum2, g_sum2, b_sum2);

                }                                                      // end if infinite light
                else if (lights[curr_light].attr & LIGHTV1_ATTR_POINT) //////////////////////////////////////
                {
                    // perform point light computations
                    // light model for point light is once again:
                    //              I0point * Clpoint
                    //  I(d)point = ___________________
                    //              kc +  kl*d + kq*d2
                    //
                    //  Where d = |p - s|
                    // thus it's almost identical to the infinite light, but attenuates as a function
                    // of distance from the point source to the surface point being lit

                    // .. normal already in vertex

                    //Write_Error("\nEntering point light....");

                    // compute vector from surface to light
                    VECTOR4D_Build(&curr_poly->tvlist[0].v, &lights[curr_light].pos, &l);

                    // compute distance and attenuation
                    dist = VECTOR4D_Length_Fast2(&l);

                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles

                    // perform the calculation for all 3 vertices

                    // vertex 0
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[0].n, &l);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (dist * atten);

                        r_sum0 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum0 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum0 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                        // r_sum0 += (lights[curr_light].c_diffuse.r * dp * r_base)/256;
                        // g_sum0 += (lights[curr_light].c_diffuse.g * dp * g_base)/256;
                        // b_sum0 += (lights[curr_light].c_diffuse.b * dp * b_base)/256;
                    } // end if

                    // vertex 1
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[1].n, &l);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (dist * atten);

                        r_sum1 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum1 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum1 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                        // r_sum1 += (lights[curr_light].c_diffuse.r * dp * r_base)/256;
                        // g_sum1 += (lights[curr_light].c_diffuse.g * dp * g_base)/256;
                        // b_sum1 += (lights[curr_light].c_diffuse.b * dp * b_base)/256;
                    } // end if

                    // vertex 2
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[2].n, &l);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (dist * atten);
                        r_sum2 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum2 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum2 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                        // r_sum2 += (lights[curr_light].c_diffuse.r * dp * r_base)/256;
                        // g_sum2 += (lights[curr_light].c_diffuse.g * dp * g_base)/256;
                        // b_sum2 += (lights[curr_light].c_diffuse.b * dp * b_base)/256;
                    } // end if

                    //Write_Error("\nexiting point light, rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,  r_sum1, g_sum1, b_sum1, r_sum2, g_sum2, b_sum2);

                }                                                           // end if point
                else if (lights[curr_light].attr & LIGHTV1_ATTR_SPOTLIGHT1) ///////////////////////////////////////
                {
                    // perform spotlight/point computations simplified model that uses
                    // point light WITH a direction to simulate a spotlight
                    // light model for point light is once again:
                    //              I0point * Clpoint
                    //  I(d)point = ___________________
                    //              kc +  kl*d + kq*d2
                    //
                    //  Where d = |p - s|
                    // thus it's almost identical to the infinite light, but attenuates as a function
                    // of distance from the point source to the surface point being lit

                    //Write_Error("\nentering spotlight1....");

                    // .. normal is already computed

                    // compute vector from surface to light
                    VECTOR4D_Build(&curr_poly->tvlist[0].v, &lights[curr_light].pos, &l);

                    // compute distance and attenuation
                    dist = VECTOR4D_Length_Fast2(&l);

                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles

                    // note that I use the direction of the light here rather than a the vector to the light
                    // thus we are taking orientation into account which is similar to the spotlight model

                    // vertex 0
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[0].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (atten);

                        r_sum0 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum0 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum0 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                    // vertex 1
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[1].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (atten);

                        r_sum1 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum1 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum1 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end i

                    // vertex 2
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[2].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (atten);

                        r_sum2 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum2 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum2 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end i

                    //Write_Error("\nexiting spotlight1, sums rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,  r_sum1, g_sum1, b_sum1,  r_sum2, g_sum2, b_sum2);

                }                                                           // end if spotlight1
                else if (lights[curr_light].attr & LIGHTV1_ATTR_SPOTLIGHT2) // simple version //////////////////////////
                {
                    // perform spot light computations
                    // light model for spot light simple version is once again:
                    //         	     I0spotlight * Clspotlight * MAX( (l . s), 0)^pf
                    // I(d)spotlight = __________________________________________
                    //               		 kc + kl*d + kq*d2
                    // Where d = |p - s|, and pf = power factor

                    // thus it's almost identical to the point, but has the extra term in the numerator
                    // relating the angle between the light source and the point on the surface

                    // .. already have normals and length are 1.0

                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles

                    //Write_Error("\nEntering spotlight2...");

                    // tons of redundant math here! lots to optimize later!

                    // vertex 0
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[0].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        // compute vector from light to surface (different from l which IS the light dir)
                        VECTOR4D_Build(&lights[curr_light].pos, &curr_poly->tvlist[0].v, &s);

                        // compute length of s (distance to light source) to normalize s for lighting calc
                        dists = VECTOR4D_Length_Fast2(&s);

                        // compute spot light term (s . l)
                        float dpsl = VECTOR4D_Dot(&s, &lights[curr_light].dir) / dists;

                        // proceed only if term is positive
                        if (dpsl > 0)
                        {
                            // compute attenuation
                            atten = (lights[curr_light].kc + lights[curr_light].kl * dists + lights[curr_light].kq * dists * dists);

                            // for speed reasons, pf exponents that are less that 1.0 are out of the question, and exponents
                            // must be integral
                            float dpsl_exp = dpsl;

                            // exponentiate for positive integral powers
                            for (int e_index = 1; e_index < (int)lights[curr_light].pf; e_index++)
                                dpsl_exp *= dpsl;

                            // now dpsl_exp holds (dpsl)^pf power which is of course (s . l)^pf

                            i = 128 * dp * dpsl_exp / (atten);

                            r_sum0 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                            g_sum0 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                            b_sum0 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);

                        } // end if

                    } // end if

                    // vertex 1
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[1].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        // compute vector from light to surface (different from l which IS the light dir)
                        VECTOR4D_Build(&lights[curr_light].pos, &curr_poly->tvlist[1].v, &s);

                        // compute length of s (distance to light source) to normalize s for lighting calc
                        dists = VECTOR4D_Length_Fast2(&s);

                        // compute spot light term (s . l)
                        float dpsl = VECTOR4D_Dot(&s, &lights[curr_light].dir) / dists;

                        // proceed only if term is positive
                        if (dpsl > 0)
                        {
                            // compute attenuation
                            atten = (lights[curr_light].kc + lights[curr_light].kl * dists + lights[curr_light].kq * dists * dists);

                            // for speed reasons, pf exponents that are less that 1.0 are out of the question, and exponents
                            // must be integral
                            float dpsl_exp = dpsl;

                            // exponentiate for positive integral powers
                            for (int e_index = 1; e_index < (int)lights[curr_light].pf; e_index++)
                                dpsl_exp *= dpsl;

                            // now dpsl_exp holds (dpsl)^pf power which is of course (s . l)^pf

                            i = 128 * dp * dpsl_exp / (atten);

                            r_sum1 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                            g_sum1 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                            b_sum1 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);

                        } // end if

                    } // end if

                    // vertex 2
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[2].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        // compute vector from light to surface (different from l which IS the light dir)
                        VECTOR4D_Build(&lights[curr_light].pos, &curr_poly->tvlist[2].v, &s);

                        // compute length of s (distance to light source) to normalize s for lighting calc
                        dists = VECTOR4D_Length_Fast2(&s);

                        // compute spot light term (s . l)
                        float dpsl = VECTOR4D_Dot(&s, &lights[curr_light].dir) / dists;

                        // proceed only if term is positive
                        if (dpsl > 0)
                        {
                            // compute attenuation
                            atten = (lights[curr_light].kc + lights[curr_light].kl * dists + lights[curr_light].kq * dists * dists);

                            // for speed reasons, pf exponents that are less that 1.0 are out of the question, and exponents
                            // must be integral
                            float dpsl_exp = dpsl;

                            // exponentiate for positive integral powers
                            for (int e_index = 1; e_index < (int)lights[curr_light].pf; e_index++)
                                dpsl_exp *= dpsl;

                            // now dpsl_exp holds (dpsl)^pf power which is of course (s . l)^pf

                            i = 128 * dp * dpsl_exp / (atten);

                            r_sum2 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                            g_sum2 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                            b_sum2 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                        } // end if

                    } // end if

                    //Write_Error("\nexiting spotlight2, sums rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,   r_sum1, g_sum1, b_sum1,  r_sum2, g_sum2, b_sum2);

                } // end if spot light

            } // end for light

            // make sure colors aren't out of range
            if (r_sum0 > 255)
                r_sum0 = 255;
            if (g_sum0 > 255)
                g_sum0 = 255;
            if (b_sum0 > 255)
                b_sum0 = 255;

            if (r_sum1 > 255)
                r_sum1 = 255;
            if (g_sum1 > 255)
                g_sum1 = 255;
            if (b_sum1 > 255)
                b_sum1 = 255;

            if (r_sum2 > 255)
                r_sum2 = 255;
            if (g_sum2 > 255)
                g_sum2 = 255;
            if (b_sum2 > 255)
                b_sum2 = 255;

            //Write_Error("\nwriting color for poly %d", poly);

            //Write_Error("\n******** final sums rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,  r_sum1, g_sum1, b_sum1, r_sum2, g_sum2, b_sum2);

            // write the colors
            curr_poly->lit_color[0] = RGB16Bit565(r_sum0, g_sum0, b_sum0);
            curr_poly->lit_color[1] = RGB16Bit565(r_sum1, g_sum1, b_sum1);
            curr_poly->lit_color[2] = RGB16Bit565(r_sum2, g_sum2, b_sum2);

        } // end if
        // else // assume POLY4DV2_ATTR_SHADE_MODE_CONSTANT
        // {
        //     // emmisive shading only, do nothing
        //     // ...
        //     curr_poly->lit_color[0] = curr_poly->color;

        //     //Write_Error("\nentering constant shader, and exiting...");

        // } // end if

    } // end for poly
}
}

// return success
return (1);

} // end Light_RENDERLIST4DV2_World16
