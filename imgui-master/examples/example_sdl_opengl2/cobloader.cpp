#include "cobloader.h"

char *Extract_Filename_From_Path(char *filepath, char *filename)
{
    // this function extracts the filename from a complete path and file
    // "../folder/.../filname.ext"
    // the function operates by scanning backward and looking for the first
    // occurance of "\" or "/" then copies the filename from there to the end
    // test of filepath is valid
    if (!filepath || strlen(filepath) == 0)
        return (NULL);

    int index_end = strlen(filepath) - 1;

    // find filename
    while ((filepath[index_end] != '\\') &&
           (filepath[index_end] != '/') &&
           (filepath[index_end] > 0))
        index_end--;

    // copy file name out into filename var
    memcpy(filename, &filepath[index_end + 1], strlen(filepath) - index_end);

    // return result
    return (filename);

} // end Extract_Filename_From_Path // end ReplaceChars

int Load_OBJECT4DV2_COB(OBJECT4DV2_PTR obj, // pointer to object
                        char *filename,     // filename of Caligari COB file
                        VECTOR4D_PTR scale, // initial scaling factors
                        VECTOR4D_PTR pos,   // initial position
                        VECTOR4D_PTR rot,   // initial rotations
                        int vertex_flags)   // flags to re-order vertices
                                            // and perform transforms
{
    // this function loads a Caligari TrueSpace .COB file object in off disk, additionally
    // it allows the caller to scale, position, and rotate the object
    // to save extra calls later for non-dynamic objects, note that this function
    // works with a OBJECT4DV2 which has support for textures, but not materials, etc,
    // however we will still parse out the material stuff and get them ready for the
    // next incarnation objects, so we can re-use this code to support those features
    // also, since this version IS going to read in the texture map and texture coordinates
    // we have a couple issues to think about, first COB format like absolute texture paths
    // we can't have that, so we will simple extract out ONLY the texture map bitmap name
    // and use the global texture path variable to build a real file path, also texture
    // coordinates are in 0..1 0..1 form, I still haven't decided if I want to use absolute
    // coordinates or 0..1 0..1, but right now the affine texture mapper uses

    // create a parser object
    CPARSERV1 parser;

    char seps[16];          // seperators for token scanning
    char token_buffer[256]; // used as working buffer for token
    char *token;            // pointer to next token

    int r, g, b; // working colors

    // cache for texture vertices
    VERTEX2DF texture_vertices[OBJECT4DV2_MAX_VERTICES];

    int num_texture_vertices = 0;

    MATRIX4X4 mat_local, // storage for local transform if user requests it in cob format
        mat_world;       // "   " for local to world " "

    // initialize matrices
    MAT_IDENTITY_4X4(&mat_local);
    MAT_IDENTITY_4X4(&mat_world);

    // Step 1: clear out the object and initialize it a bit
    memset(obj, 0, sizeof(OBJECT4DV2));

    // set state of object to active and visible
    obj->state = OBJECT4DV2_STATE_ACTIVE | OBJECT4DV2_STATE_VISIBLE;

    // set number of frames
    obj->num_frames = 1;
    obj->curr_frame = 0;
    obj->attr = OBJECT4DV2_ATTR_SINGLE_FRAME;

    // set position of object is caller requested position
    if (pos)
    {
        // set position of object
        obj->world_pos.x = pos->x;
        obj->world_pos.y = pos->y;
        obj->world_pos.z = pos->z;
        obj->world_pos.w = pos->w;
    } // end
    else
    {
        // set it to (0,0,0,1)
        obj->world_pos.x = 0;
        obj->world_pos.y = 0;
        obj->world_pos.z = 0;
        obj->world_pos.w = 1;
    } // end else

    // Step 2: open the file for reading using the parser
    if (!parser.Open(filename))
    {
        //Write_Error("Couldn't open .COB file %s.", filename);
        return (0);
    } // end if

    // Step 3:

    // lets find the name of the object first
    while (1)
    {
        // get the next line, we are looking for "Name"
        if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
        {
            //Write_Error("Image 'name' not found in .COB file %s.", filename);
            return (0);
        } // end if

        // check for pattern?
        if (parser.Pattern_Match(parser.buffer, "['Name'] [s>0]"))
        {
            // name should be in second string variable, index 1
            strcpy(obj->name, parser.pstrings[1]);
            //Write_Error("\nCOB Reader Object Name: %s", obj->name);

            break;
        } // end if

    } // end while

    // step 4: get local and world transforms and store them

    // center 0 0 0
    // x axis 1 0 0
    // y axis 0 1 0
    // z axis 0 0 1

    while (1)
    {
        // get the next line, we are looking for "center"
        if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
        {
            //Write_Error("Center not found in .COB file %s.", filename);
            return (0);
        } // end if

        // check for pattern?
        if (parser.Pattern_Match(parser.buffer, "['center'] [f] [f] [f]"))
        {
            // the "center" holds the translation factors, so place in
            // last row of homogeneous matrix, note that these are row vectors
            // that we need to drop in each column of matrix
            mat_local.M[3][0] = -parser.pfloats[0]; // center x
            mat_local.M[3][1] = -parser.pfloats[1]; // center y
            mat_local.M[3][2] = -parser.pfloats[2]; // center z

            // ok now, the next 3 lines should be the x,y,z transform vectors
            // so build up

            // "x axis"
            parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS);
            parser.Pattern_Match(parser.buffer, "['x'] ['axis'] [f] [f] [f]");

            // place row in x column of transform matrix
            mat_local.M[0][0] = parser.pfloats[0]; // rxx
            mat_local.M[1][0] = parser.pfloats[1]; // rxy
            mat_local.M[2][0] = parser.pfloats[2]; // rxz

            // "y axis"
            parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS);
            parser.Pattern_Match(parser.buffer, "['y'] ['axis'] [f] [f] [f]");

            // place row in y column of transform matrix
            mat_local.M[0][1] = parser.pfloats[0]; // ryx
            mat_local.M[1][1] = parser.pfloats[1]; // ryy
            mat_local.M[2][1] = parser.pfloats[2]; // ryz

            // "z axis"
            parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS);
            parser.Pattern_Match(parser.buffer, "['z'] ['axis'] [f] [f] [f]");

            // place row in z column of transform matrix
            mat_local.M[0][2] = parser.pfloats[0]; // rzx
            mat_local.M[1][2] = parser.pfloats[1]; // rzy
            mat_local.M[2][2] = parser.pfloats[2]; // rzz

            // Print_Mat_4X4(&mat_local, "Local COB Matrix:");

            break;
        } // end if

    } // end while

    // now "Transform"
    while (1)
    {
        // get the next line, we are looking for "Transform"
        if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
        {
            //Write_Error("Transform not found in .COB file %s.", filename);
            return (0);
        } // end if

        // check for pattern?
        if (parser.Pattern_Match(parser.buffer, "['Transform']"))
        {

            // "x axis"
            parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS);
            parser.Pattern_Match(parser.buffer, "[f] [f] [f]");

            // place row in x column of transform matrix
            mat_world.M[0][0] = parser.pfloats[0]; // rxx
            mat_world.M[1][0] = parser.pfloats[1]; // rxy
            mat_world.M[2][0] = parser.pfloats[2]; // rxz

            // "y axis"
            parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS);
            parser.Pattern_Match(parser.buffer, "[f] [f] [f]");

            // place row in y column of transform matrix
            mat_world.M[0][1] = parser.pfloats[0]; // ryx
            mat_world.M[1][1] = parser.pfloats[1]; // ryy
            mat_world.M[2][1] = parser.pfloats[2]; // ryz

            // "z axis"
            parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS);
            parser.Pattern_Match(parser.buffer, "[f] [f] [f]");

            // place row in z column of transform matrix
            mat_world.M[0][2] = parser.pfloats[0]; // rzx
            mat_world.M[1][2] = parser.pfloats[1]; // rzy
            mat_world.M[2][2] = parser.pfloats[2]; // rzz

            // Print_Mat_4X4(&mat_world, "World COB Matrix:");

            // no need to read in last row, since it's always 0,0,0,1 and we don't use it anyway
            break;

        } // end if

    } // end while

    // step 6: get number of vertices and polys in object
    while (1)
    {
        // get the next line, we are looking for "World Vertices"
        if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
        {
            //Write_Error("'World Vertices' line not found in .COB file %s.", filename);
            return (0);
        } // end if

        // check for pattern?
        if (parser.Pattern_Match(parser.buffer, "['World'] ['Vertices'] [i]"))
        {
            // simply extract the number of vertices from the pattern matching
            // output arrays
            obj->num_vertices = parser.pints[0];

            //Write_Error("\nCOB Reader Num Vertices: %d", obj->num_vertices);
            break;

        } // end if

    } // end while

    // allocate the memory for the vertices and number of polys (unknown, so use 3*num_vertices)
    // the call parameters are redundant in this case, but who cares
    if (!Init_OBJECT4DV2(obj, // object to allocate
                         obj->num_vertices,
                         obj->num_vertices * 3,
                         obj->num_frames))
    {
        //Write_Error("\nASC file error with file %s (can't allocate memory).",filename);
    } // end if

    // Step 7: load the vertex list
    // now read in vertex list, format:
    // "d.d d.d d.d"
    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        // hunt for vertex
        while (1)
        {
            // get the next vertex
            if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
            {
                //Write_Error("\nVertex list ended abruptly! in .COB file %s.", filename);
                return (0);
            } // end if

            // check for pattern?
            if (parser.Pattern_Match(parser.buffer, "[f] [f] [f]"))
            {
                // at this point we have the x,y,z in the the pfloats array locations 0,1,2
                obj->vlist_local[vertex].x = parser.pfloats[0];
                obj->vlist_local[vertex].y = parser.pfloats[1];
                obj->vlist_local[vertex].z = parser.pfloats[2];
                obj->vlist_local[vertex].w = 1;

                // do vertex swapping right here, allow muliple swaps, why not!
                // defines for vertex re-ordering flags

                //#define VERTEX_FLAGS_INVERT_X   1    // inverts the Z-coordinates
                //#define VERTEX_FLAGS_INVERT_Y   2    // inverts the Z-coordinates
                //#define VERTEX_FLAGS_INVERT_Z   4    // inverts the Z-coordinates
                //#define VERTEX_FLAGS_SWAP_YZ    8    // transforms a RHS model to a LHS model
                //#define VERTEX_FLAGS_SWAP_XZ    16   // ???
                //#define VERTEX_FLAGS_SWAP_XY    32
                //#define VERTEX_FLAGS_INVERT_WINDING_ORDER 64  // invert winding order from cw to ccw or ccw to cc
                //#define VERTEX_FLAGS_TRANSFORM_LOCAL         512   // if file format has local transform then do it!
                //#define VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD  1024  // if file format has local to world then do it!

                VECTOR4D temp_vector; // temp for calculations

                // now apply local and world transformations encoded in COB format
                if (vertex_flags & VERTEX_FLAGS_TRANSFORM_LOCAL)
                {
                    Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex].v, &mat_local, &temp_vector);
                    VECTOR4D_COPY(&obj->vlist_local[vertex].v, &temp_vector);
                } // end if

                if (vertex_flags & VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD)
                {
                    Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex].v, &mat_world, &temp_vector);
                    VECTOR4D_COPY(&obj->vlist_local[vertex].v, &temp_vector);
                } // end if

                float temp_f; // used for swapping

                // invert signs?
                if (vertex_flags & VERTEX_FLAGS_INVERT_X)
                    obj->vlist_local[vertex].x = -obj->vlist_local[vertex].x;

                if (vertex_flags & VERTEX_FLAGS_INVERT_Y)
                    obj->vlist_local[vertex].y = -obj->vlist_local[vertex].y;

                if (vertex_flags & VERTEX_FLAGS_INVERT_Z)
                    obj->vlist_local[vertex].z = -obj->vlist_local[vertex].z;

                // swap any axes?
                if (vertex_flags & VERTEX_FLAGS_SWAP_YZ)
                    SWAP(obj->vlist_local[vertex].y, obj->vlist_local[vertex].z, temp_f);

                if (vertex_flags & VERTEX_FLAGS_SWAP_XZ)
                    SWAP(obj->vlist_local[vertex].x, obj->vlist_local[vertex].z, temp_f);

                if (vertex_flags & VERTEX_FLAGS_SWAP_XY)
                    SWAP(obj->vlist_local[vertex].x, obj->vlist_local[vertex].y, temp_f);

                // scale vertices
                if (scale)
                {
                    obj->vlist_local[vertex].x *= scale->x;
                    obj->vlist_local[vertex].y *= scale->y;
                    obj->vlist_local[vertex].z *= scale->z;
                } // end if

                //Write_Error("\nVertex %d = %f, %f, %f, %f", vertex,obj->vlist_local[vertex].x, obj->vlist_local[vertex].y, obj->vlist_local[vertex].z,obj->vlist_local[vertex].w);

                // set point field in this vertex, we need that at least
                SET_BIT(obj->vlist_local[vertex].attr, VERTEX4DTV1_ATTR_POINT);

                // found vertex, break out of while for next pass
                break;

            } // end if

        } // end while

    } // end for vertex

    // compute average and max radius
    Compute_OBJECT4DV2_Radius(obj);

    //Write_Error("\nObject average radius = %f, max radius = %f",

    // step 8: get number of texture vertices
    while (1)
    {
        // get the next line, we are looking for "Texture Vertices ddd"
        if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
        {
            //Write_Error("'Texture Vertices' line not found in .COB file %s.", filename);
            return (0);
        } // end if

        // check for pattern?
        if (parser.Pattern_Match(parser.buffer, "['Texture'] ['Vertices'] [i]"))
        {
            // simply extract the number of texture vertices from the pattern matching
            // output arrays
            num_texture_vertices = parser.pints[0];

            //Write_Error("\nCOB Reader Texture Vertices: %d", num_texture_vertices);
            break;

        } // end if

    } // end while

    // Step 9: load the texture vertex list in format "U V"
    // "d.d d.d"
    for (int tvertex = 0; tvertex < num_texture_vertices; tvertex++)
    {
        // hunt for texture
        while (1)
        {
            // get the next vertex
            if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
            {
                //Write_Error("\nTexture Vertex list ended abruptly! in .COB file %s.", filename);
                return (0);
            } // end if

            // check for pattern?
            if (parser.Pattern_Match(parser.buffer, "[f] [f]"))
            {
                // at this point we have the U V in the the pfloats array locations 0,1 for this
                // texture vertex, store in texture coordinate list
                // note texture coords are in 0..1 format, and must be scaled to texture size
                // after we load the texture
                obj->tlist[tvertex].x = parser.pfloats[0];
                obj->tlist[tvertex].y = parser.pfloats[1];

                //Write_Error("\nTexture Vertex %d: U=%f, V=%f", tvertex, obj->tlist[tvertex].x, obj->tlist[tvertex].y );

                // found vertex, break out of while for next pass
                break;

            } // end if

        } // end while

    } // end for

    // when we load in the polygons then we will copy the texture vertices into the polygon
    // vertices assuming that each vertex has a SINGLE texture coordinate, this means that
    // you must NOT use multiple textures on an object! in other words think "skin" this is
    // inline with Quake II md2 format, in 99% of the cases a single object can be textured
    // with a single skin and the texture coordinates can be unique for each vertex and 1:1

    int poly_material[OBJECT4DV2_MAX_POLYS]; // this holds the material index for each polygon
                                             // we need these indices since when reading the file
                                             // we read the polygons BEFORE the materials, so we need
                                             // this data, so we can go back later and extract the material
                                             // that each poly WAS assigned and get the colors out, since
                                             // objects and polygons do not currently support materials

    int material_index_referenced[MAX_MATERIALS]; // used to track if an index has been used yet as a material
                                                  // reference. since we don't know how many materials, we need
                                                  // a way to count them up, but if we have seen a material reference
                                                  // more than once then we don't increment the total number of materials
                                                  // this array is for this

    // clear out reference array
    memset(material_index_referenced, 0, sizeof(material_index_referenced));

    // step 10: load in the polygons
    // poly list starts off with:
    // "Faces ddd:"
    while (1)
    {
        // get next line
        if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
        {
            //Write_Error("\n'Faces' line not found in .COB file %s.", filename);
            return (0);
        } // end if

        // check for pattern?
        if (parser.Pattern_Match(parser.buffer, "['Faces'] [i]"))
        {
            //Write_Error("\nCOB Reader found face list in .COB file %s.", filename);

            // finally set number of polys
            obj->num_polys = parser.pints[0];

            break;
        } // end if
    }     // end while

    // now read each face in format:
    // Face verts nn flags ff mat mm
    // the nn is the number of vertices, always 3
    // the ff is the flags, unused for now, has to do with holes
    // the mm is the material index number

    int poly_surface_desc = 0;    // ASC surface descriptor/material in this case
    int poly_num_verts = 0;       // number of vertices for current poly (always 3)
    int num_materials_object = 0; // number of materials for this object

    for (int poly = 0; poly < obj->num_polys; poly++)
    {
        //Write_Error("\nPolygon %d:", poly);
        // hunt until next face is found
        while (1)
        {
            // get the next polygon face
            if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
            {
                //Write_Error("\nface list ended abruptly! in .COB file %s.", filename);
                return (0);
            } // end if

            // check for pattern?
            if (parser.Pattern_Match(parser.buffer, "['Face'] ['verts'] [i] ['flags'] [i] ['mat'] [i]"))
            {
                // at this point we have the number of vertices for the polygon, the flags, and it's material index
                // in the integer output array locations 0,1,2

                // store the material index for this polygon for retrieval later, but make sure adjust the
                // the index to take into consideration that the data in parser.pints[2] is 0 based, and we need
                // an index relative to the entire library, so we simply need to add num_materials to offset the
                // index properly, but we will leave this reference zero based for now... and fix up later
                poly_material[poly] = parser.pints[2];

                // update the reference array
                if (material_index_referenced[poly_material[poly]] == 0)
                {
                    // mark as referenced
                    material_index_referenced[poly_material[poly]] = 1;

                    // increment total number of materials for this object
                    num_materials_object++;
                } // end if

                // test if number of vertices is 3
                if (parser.pints[0] != 3)
                {
                    //Write_Error("\nface not a triangle! in .COB file %s.", filename);
                    return (0);
                } // end if

                // now read out the vertex indices and texture indices format:
                // <vindex0, tindex0>  <vindex1, tindex1> <vindex1, tindex1>
                if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                {
                    //Write_Error("\nface list ended abruptly! in .COB file %s.", filename);
                    return (0);
                } // end if

                // lets replace ",<>" with ' ' to make extraction easy
                ReplaceChars(parser.buffer, parser.buffer, ",<>", ' ');
                parser.Pattern_Match(parser.buffer, "[i] [i] [i] [i] [i] [i]");

                // 0,2,4 holds vertex indices
                // 1,3,5 holds texture indices

                // insert polygon, check for winding order invert
                if (vertex_flags & VERTEX_FLAGS_INVERT_WINDING_ORDER)
                {
                    poly_num_verts = 3;
                    obj->plist[poly].vert[0] = parser.pints[4];
                    obj->plist[poly].vert[1] = parser.pints[2];
                    obj->plist[poly].vert[2] = parser.pints[0];

                    // now copy the texture coordinates into the vertices, this
                    // may not be needed if the polygon doesn't have texture mapping
                    // enabled, etc.,

                    // so here's the deal the texture coordinates that
                    // map to vertex 0,1,2 have indices stored in the odd
                    // numbered pints[] locations, so we simply need to copy
                    // the right texture coordinate into the right vertex
                    obj->plist[poly].text[0] = parser.pints[5];
                    obj->plist[poly].text[1] = parser.pints[3];
                    obj->plist[poly].text[2] = parser.pints[1];

                } // end if
                else
                { // leave winding order alone
                    poly_num_verts = 3;
                    obj->plist[poly].vert[0] = parser.pints[0];
                    obj->plist[poly].vert[1] = parser.pints[2];
                    obj->plist[poly].vert[2] = parser.pints[4];

                    // now copy the texture coordinates into the vertices, this
                    // may not be needed if the polygon doesn't have texture mapping
                    // enabled, etc.,

                    // so here's the deal the texture coordinates that
                    // map to vertex 0,1,2 have indices stored in the odd
                    // numbered pints[] locations, so we simply need to copy
                    // the right texture coordinate into the right vertex
                    obj->plist[poly].text[0] = parser.pints[1];
                    obj->plist[poly].text[1] = parser.pints[3];
                    obj->plist[poly].text[2] = parser.pints[5];

                } // end else

                // point polygon vertex list to object's vertex list
                // note that this is redundant since the polylist is contained
                // within the object in this case and its up to the user to select
                // whether the local or transformed vertex list is used when building up
                // polygon geometry, might be a better idea to set to NULL in the context
                // of polygons that are part of an object
                obj->plist[poly].vlist = obj->vlist_local;

                // set texture coordinate list, this is needed
                obj->plist[poly].tlist = obj->tlist;

                // set polygon to active
                obj->plist[poly].state = POLY4DV2_STATE_ACTIVE;

                // found the face, break out of while for another pass
                break;

            } // end if

        } // end while

    } // end for poly

    // now find materials!!! and we are out of here!
    for (int curr_material = 0; curr_material < num_materials_object; curr_material++)
    {
        // hunt for the material header "mat# ddd"
        while (1)
        {
            // get the next polygon material
            if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
            {
                //Write_Error("\nmaterial list ended abruptly! in .COB file %s.", filename);
                return (0);
            } // end if

            // check for pattern?
            if (parser.Pattern_Match(parser.buffer, "['mat#'] [i]"))
            {
                // extract the material that is being defined
                int material_index = parser.pints[0];

                // get color of polygon, although it might be irrelevant for a textured surface
                while (1)
                {
                    // get the next line
                    if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                    {
                        //Write_Error("\nRGB color ended abruptly! in .COB file %s.", filename);
                        return (0);
                    } // end if

                    // replace the , comma's if there are any with spaces
                    ReplaceChars(parser.buffer, parser.buffer, ",", ' ', 1);

                    // look for "rgb float,float,float"
                    if (parser.Pattern_Match(parser.buffer, "['rgb'] [f] [f] [f]"))
                    {
                        // extract data and store color in material libary
                        // pfloats[] 0,1,2,3, has data
                        materials[material_index + num_materials].color.r = (int)(parser.pfloats[0] * 255 + 0.5);
                        materials[material_index + num_materials].color.g = (int)(parser.pfloats[1] * 255 + 0.5);
                        materials[material_index + num_materials].color.b = (int)(parser.pfloats[2] * 255 + 0.5);

                        break; // while looking for rgb
                    }          // end if

                } // end while

                // extract out lighting constants for the heck of it, they are on a line like this:
                // "alpha float ka float ks float exp float ior float"
                // alpha is transparency           0 - 1
                // ka is ambient coefficient       0 - 1
                // ks is specular coefficient      0 - 1
                // exp is highlight power exponent 0 - 1
                // ior is index of refraction (unused)

                // although our engine will have minimal support for these, we might as well get them
                while (1)
                {
                    // get the next line
                    if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                    {
                        //Write_Error("\nmaterial properties ended abruptly! in .COB file %s.", filename);
                        return (0);
                    } // end if

                    // look for "alpha float ka float ks float exp float ior float"
                    if (parser.Pattern_Match(parser.buffer, "['alpha'] [f] ['ka'] [f] ['ks'] [f] ['exp'] [f]"))
                    {
                        // extract data and store in material libary
                        // pfloats[] 0,1,2,3, has data
                        materials[material_index + num_materials].color.a = (UCHAR)(parser.pfloats[0] * 255 + 0.5);
                        materials[material_index + num_materials].ka = parser.pfloats[1];
                        materials[material_index + num_materials].kd = 1; // hard code for now
                        materials[material_index + num_materials].ks = parser.pfloats[2];
                        materials[material_index + num_materials].power = parser.pfloats[3];

                        // compute material reflectivities in pre-multiplied format to help engine
                        for (int rgb_index = 0; rgb_index < 3; rgb_index++)
                        {
                            // ambient reflectivity
                            materials[material_index + num_materials].ra.rgba_M[rgb_index] =
                                ((UCHAR)(materials[material_index + num_materials].ka *
                                             (float)materials[material_index + num_materials].color.rgba_M[rgb_index] +
                                         0.5));

                            // diffuse reflectivity
                            materials[material_index + num_materials].rd.rgba_M[rgb_index] =
                                ((UCHAR)(materials[material_index + num_materials].kd *
                                             (float)materials[material_index + num_materials].color.rgba_M[rgb_index] +
                                         0.5));

                            // specular reflectivity
                            materials[material_index + num_materials].rs.rgba_M[rgb_index] =
                                ((UCHAR)(materials[material_index + num_materials].ks *
                                             (float)materials[material_index + num_materials].color.rgba_M[rgb_index] +
                                         0.5));

                        } // end for rgb_index

                        break;
                    } // end if

                } // end while

                // now we need to know the shading model, it's a bit tricky, we need to look for the lines
                // "Shader class: color" first, then after this line is:
                // "Shader name: "xxxxxx" (xxxxxx) "
                // where the xxxxx part will be "plain color" and "plain" for colored polys
                // or "texture map" and "caligari texture"  for textures
                // THEN based on that we hunt for "Shader class: reflectance" which is where the type
                // of shading is encoded, we look for the "Shader name: "xxxxxx" (xxxxxx) " again,
                // and based on it's value we map it to our shading system as follows:
                // "constant" -> MATV1_ATTR_SHADE_MODE_CONSTANT
                // "matte"    -> MATV1_ATTR_SHADE_MODE_FLAT
                // "plastic"  -> MATV1_ATTR_SHADE_MODE_GOURAUD
                // "phong"    -> MATV1_ATTR_SHADE_MODE_FASTPHONG
                // and in the case that in the "color" class, we found a "texture map" then the "shading mode" is
                // "texture map" -> MATV1_ATTR_SHADE_MODE_TEXTURE
                // which must be logically or'ed with the other previous modes

                //  look for the "shader class: color"
                while (1)
                {
                    // get the next line
                    if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                    {
                        //Write_Error("\nshader class ended abruptly! in .COB file %s.", filename);
                        return (0);
                    } // end if

                    if (parser.Pattern_Match(parser.buffer, "['Shader'] ['class:'] ['color']"))
                    {
                        break;
                    } // end if

                } // end while

                // now look for the shader name for this class
                // Shader name: "plain color" or Shader name: "texture map"
                while (1)
                {
                    // get the next line
                    if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                    {
                        //Write_Error("\nshader name ended abruptly! in .COB file %s.", filename);
                        return (0);
                    } // end if

                    // replace the " with spaces
                    ReplaceChars(parser.buffer, parser.buffer, "\"", ' ', 1);

                    // is this a "plain color" poly?
                    if (parser.Pattern_Match(parser.buffer, "['Shader'] ['name:'] ['plain'] ['color']"))
                    {
                        // not much to do this is default, we need to wait for the reflectance type
                        // to tell us the shading mode

                        break;
                    } // end if

                    // is this a "texture map" poly?
                    if (parser.Pattern_Match(parser.buffer, "['Shader'] ['name:'] ['texture'] ['map']"))
                    {
                        // set the texture mapping flag in material
                        SET_BIT(materials[material_index + num_materials].attr, MATV1_ATTR_SHADE_MODE_TEXTURE);

                        // almost done, we need the file name of the darn texture map, its in this format:
                        // file name: string "D:\Source\..\models\textures\wall01.bmp"

                        // of course the filename in the quotes will change
                        // so lets hunt until we find it...
                        while (1)
                        {
                            // get the next line
                            if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                            {
                                //Write_Error("\ncouldnt find texture name! in .COB file %s.", filename);
                                return (0);
                            } // end if

                            // replace the " with spaces
                            ReplaceChars(parser.buffer, parser.buffer, "\"", ' ', 1);

                            // is this the file name?
                            if (parser.Pattern_Match(parser.buffer, "['file'] ['name:'] ['string']"))
                            {
                                // and save the FULL filename (useless though since its the path from the
                                // machine that created it, but later we might want some of the info).
                                // filename and path starts at char position 19, 0 indexed
                                memcpy(materials[material_index + num_materials].texture_file, &parser.buffer[18], strlen(parser.buffer) - 18 + 2);

                                // the OBJECT4DV2 is only allowed a single texture, although we are loading in all
                                // the materials, if this is the first texture map, load it, and set a flag disallowing
                                // any more texture loads for the object
                                if (!obj->texture)
                                {
                                    // step 1: allocate memory for bitmap
                                    obj->texture = (BITMAP_IMAGE_PTR)malloc(sizeof(BITMAP_IMAGE));

                                    // load the texture, just use the final file name and the absolute global
                                    // texture path
                                    char filename[80];
                                    char path_filename[80];
                                    // get the filename
                                    Extract_Filename_From_Path(materials[material_index + num_materials].texture_file, filename);

                                    // build the filename with root path
                                    strcpy(path_filename, texture_path);
                                    strcat(path_filename, filename);

                                    // buffer now holds final texture path and file name
                                    // load the bitmap(8/16 bit)
                                    // Load_Bitmap_File(&bitmap16bit, path_filename);

                                    // create a proper size and bitdepth bitmap
                                    // Create_Bitmap(obj->texture, 0, 0,
                                    //   bitmap16bit.bitmapinfoheader.biWidth,
                                    //   bitmap16bit.bitmapinfoheader.biHeight,
                                    //   bitmap16bit.bitmapinfoheader.biBitCount);

                                    // load the bitmap image (later make this 8/16 bit)
                                    if (obj->texture->bpp == 16)
                                        // Load_Image_Bitmap16(obj->texture, &bitmap16bit, 0, 0, BITMAP_EXTRACT_MODE_ABS);
                                        // else
                                        // {
                                        //     // Load_Image_Bitmap(obj->texture, &bitmap16bit, 0, 0, BITMAP_EXTRACT_MODE_ABS);
                                        // } // end else 8 bit

                                        // done, so unload the bitmap
                                        // Unload_Bitmap_File(&bitmap16bit);

                                        // flag object as having textures
                                        SET_BIT(obj->attr, OBJECT4DV2_ATTR_TEXTURES);

                                } // end if

                                break;
                            } // end if

                        } // end while

                        break;
                    } // end if

                } // end while

                // alright, finally! Now we need to know what the actual shader type, now in the COB format
                // I have decided that in the "reflectance" class that's where we will look at what kind
                // of shader is supposed to be used on the polygon

                //  look for the "Shader class: reflectance"
                while (1)
                {
                    // get the next line
                    if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                    {
                        //Write_Error("\nshader reflectance class not found in .COB file %s.", filename);
                        return (0);
                    } // end if

                    // look for "Shader class: reflectance"
                    if (parser.Pattern_Match(parser.buffer, "['Shader'] ['class:'] ['reflectance']"))
                    {
                        // now we know the next "shader name" is what we are looking for so, break

                        break;
                    } // end if

                } // end while

                // looking for "Shader name: "xxxxxx" (xxxxxx) " again,
                // and based on it's value we map it to our shading system as follows:
                // "constant" -> MATV1_ATTR_SHADE_MODE_CONSTANT
                // "matte"    -> MATV1_ATTR_SHADE_MODE_FLAT
                // "plastic"  -> MATV1_ATTR_SHADE_MODE_GOURAUD
                // "phong"    -> MATV1_ATTR_SHADE_MODE_FASTPHONG
                // and in the case that in the "color" class, we found a "texture map" then the "shading mode" is
                // "texture map" -> MATV1_ATTR_SHADE_MODE_TEXTURE
                // which must be logically or'ed with the other previous modes
                while (1)
                {
                    // get the next line
                    if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                    {
                        //Write_Error("\nshader name ended abruptly! in .COB file %s.", filename);
                        return (0);
                    } // end if

                    // get rid of those quotes
                    ReplaceChars(parser.buffer, parser.buffer, "\"", ' ', 1);

                    // did we find the name?
                    if (parser.Pattern_Match(parser.buffer, "['Shader'] ['name:'] [s>0]"))
                    {
                        // figure out which shader to use
                        if (strcmp(parser.pstrings[2], "constant") == 0)
                        {
                            // set the shading mode flag in material
                            SET_BIT(materials[material_index + num_materials].attr, MATV1_ATTR_SHADE_MODE_CONSTANT);
                        } // end if
                        else if (strcmp(parser.pstrings[2], "matte") == 0)
                        {
                            // set the shading mode flag in material
                            SET_BIT(materials[material_index + num_materials].attr, MATV1_ATTR_SHADE_MODE_FLAT);
                        } // end if
                        else if (strcmp(parser.pstrings[2], "plastic") == 0)
                        {
                            // set the shading mode flag in material
                            SET_BIT(materials[curr_material + num_materials].attr, MATV1_ATTR_SHADE_MODE_GOURAUD);
                        } // end if
                        else if (strcmp(parser.pstrings[2], "phong") == 0)
                        {
                            // set the shading mode flag in material
                            SET_BIT(materials[material_index + num_materials].attr, MATV1_ATTR_SHADE_MODE_FASTPHONG);
                        } // end if
                        else
                        {
                            // set the shading mode flag in material
                            SET_BIT(materials[material_index + num_materials].attr, MATV1_ATTR_SHADE_MODE_FLAT);
                        } // end else

                        break;
                    } // end if

                } // end while

                // found the material, break out of while for another pass
                break;

            } // end if found material

        } // end while looking for mat#1

    } // end for curr_material

    // at this point poly_material[] holds all the indices for the polygon materials (zero based, so they need fix up)
    // and we must access the materials array to fill in each polygon with the polygon color, etc.
    // now that we finally have the material libary loaded
    for (int curr_poly = 0; curr_poly < obj->num_polys; curr_poly++)
    {

        // fix up offset
        poly_material[curr_poly] = poly_material[curr_poly] + num_materials;

        // we need to know what color depth we are dealing with, so check
        // the bits per pixel, this assumes that the system has already
        // made the call to DDraw_Init() or set the bit depth
        // if (screen_bpp == 16)
        if (true)
        {
            // cool, 16 bit mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV1_ATTR_RGB16);

            // test if this is a textured poly, if so override the color to WHITE,
            // so we get maximum reflection in lighting stage
            if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_TEXTURE)
                obj->plist[curr_poly].color = RGB16Bit565(255, 255, 255);
            else
                obj->plist[curr_poly].color = RGB16Bit565(materials[poly_material[curr_poly]].color.r,
                                                          materials[poly_material[curr_poly]].color.g,
                                                          materials[poly_material[curr_poly]].color.b);
            //Write_Error("\nPolygon 16-bit");
        } // end
        else
        {
            // 8 bit mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV1_ATTR_8BITCOLOR);

            // test if this is a textured poly, if so override the color to WHITE,
            // so we get maximum reflection in lighting stage
            if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_TEXTURE)
                // obj->plist[curr_poly].color = RGBto8BitIndex(255, 255, 255, palette, 0);
                obj->plist[curr_poly].color = RGB16Bit565(255, 255, 255);

            else
                obj->plist[curr_poly].color = RGB16Bit565(materials[poly_material[curr_poly]].color.r, materials[poly_material[curr_poly]].color.g, materials[poly_material[curr_poly]].color.b);
            // obj->plist[curr_poly].color = RGBto8BitIndex(materials[poly_material[curr_poly]].color.r,
            //                                              materials[poly_material[curr_poly]].color.g,
            //                                              materials[poly_material[curr_poly]].color.b,
            //                                              palette, 0);

            //Write_Error("\nPolygon 8-bit, index=%d", obj->plist[curr_poly].color);
        } // end else

        // now set all the shading flags
        // figure out which shader to use
        if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_CONSTANT)
        {
            // set shading mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_SHADE_MODE_CONSTANT);
        } // end if
        else if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_FLAT)
        {
            // set shading mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_SHADE_MODE_FLAT);
        } // end if
        else if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_GOURAUD)
        {
            // set shading mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_SHADE_MODE_GOURAUD);

            // going to need vertex normals!
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[0]].attr, VERTEX4DTV1_ATTR_NORMAL);
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[1]].attr, VERTEX4DTV1_ATTR_NORMAL);
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[2]].attr, VERTEX4DTV1_ATTR_NORMAL);

        } // end if
        else if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_FASTPHONG)
        {
            // set shading mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_SHADE_MODE_FASTPHONG);

            // going to need vertex normals!
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[0]].attr, VERTEX4DTV1_ATTR_NORMAL);
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[1]].attr, VERTEX4DTV1_ATTR_NORMAL);
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[2]].attr, VERTEX4DTV1_ATTR_NORMAL);
        } // end if
        else
        {
            // set shading mode to default flat
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_SHADE_MODE_FLAT);

        } // end if

        if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_TEXTURE)
        {
            // set shading mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_SHADE_MODE_TEXTURE);

            // apply texture to this polygon
            obj->plist[curr_poly].texture = obj->texture;

            // set texture coordinate attributes
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[0]].attr, VERTEX4DTV1_ATTR_TEXTURE);
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[1]].attr, VERTEX4DTV1_ATTR_TEXTURE);
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[2]].attr, VERTEX4DTV1_ATTR_TEXTURE);

        } // end if

        // set the material mode to ver. 1.0 emulation (for now only!!!)
        SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_DISABLE_MATERIAL);

    } // end for curr_poly

    // local object materials have been added to database, update total materials in system
    num_materials += num_materials_object;

    // now fix up all texture coordinates
    if (obj->texture)
    {
        for (int tvertex = 0; tvertex < num_texture_vertices; tvertex++)
        {
            // step 1: scale the texture coordinates by the texture size
            int texture_size = obj->texture->width;

            // scale 0..1 to 0..texture_size-1
            obj->tlist[tvertex].x *= (texture_size - 1);
            obj->tlist[tvertex].y *= (texture_size - 1);

            // now test for vertex transformation flags
            if (vertex_flags & VERTEX_FLAGS_INVERT_TEXTURE_U)
            {
                obj->tlist[tvertex].x = (texture_size - 1) - obj->tlist[tvertex].x;
            } // end if

            if (vertex_flags & VERTEX_FLAGS_INVERT_TEXTURE_V)
            {
                obj->tlist[tvertex].y = (texture_size - 1) - obj->tlist[tvertex].y;
            } // end if

            if (vertex_flags & VERTEX_FLAGS_INVERT_SWAP_UV)
            {
                float temp;
                SWAP(obj->tlist[tvertex].x, obj->tlist[tvertex].y, temp);
            } // end if

        } // end for

    } // end if there was a texture loaded for this object

#ifdef DEBUG_ON
    for (curr_material = 0; curr_material < num_materials; curr_material++)
    {
        //Write_Error("\nMaterial %d", curr_material);

        //Write_Error("\nint  state    = %d", materials[curr_material].state);
        //Write_Error("\nint  id       = %d", materials[curr_material].id);
        //Write_Error("\nchar name[64] = %s", materials[curr_material].name);
        //Write_Error("\nint  attr     = %d", materials[curr_material].attr);
        //Write_Error("\nint r         = %d", materials[curr_material].color.r);
        //Write_Error("\nint g         = %d", materials[curr_material].color.g);
        //Write_Error("\nint b         = %d", materials[curr_material].color.b);
        //Write_Error("\nint alpha     = %d", materials[curr_material].color.a);
        //Write_Error("\nint color     = %d", materials[curr_material].attr);
        //Write_Error("\nfloat ka      = %f", materials[curr_material].ka);
        //Write_Error("\nkd            = %f", materials[curr_material].kd);
        //Write_Error("\nks            = %f", materials[curr_material].ks);
        //Write_Error("\npower         = %f", materials[curr_material].power);
        //Write_Error("\nchar texture_file = %s\n", materials[curr_material].texture_file);
    } // end for curr_material
#endif

    // now that we know the correct number of polygons, we must allocate memory for them
    // and fix up the object, this is a hack, but the file formats are so stupid by not
    // all starting with NUM_VERTICES, NUM_POLYGONS -- that would make everyone's life
    // easier!

#if 0

// step 1: allocate memory for the polygons
POLY4DV2_PTR plist_temp = NULL;

// allocate memory for polygon list, the correct number of polys was overwritten
// into the object during parsing, so we can use the num_polys field
if (!(plist_temp = (POLY4DV2_PTR)malloc(sizeof(POLY4DV2)*obj->num_polys)))
   return(0);

// step 2:  now copy the polygons into the correct list
memcpy((void *)plist_temp, (void *)obj->plist, sizeof(POLY4DV2));

// step 3: now free the old memory and fix the pointer
free(obj->plist);

// now fix the pointer
obj->plist = plist_temp;

#endif

    // compute the polygon normal lengths
    Compute_OBJECT4DV2_Poly_Normals(obj);

    // compute vertex normals for any gouraud shaded polys
    Compute_OBJECT4DV2_Vertex_Normals(obj);

    // return success
    return (1);

} // end Load_OBJECT4DV2_COB
