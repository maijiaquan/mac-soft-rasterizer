#pragma once

#include "ds.h"
char *Extract_Filename_From_Path(char *filepath, char *filename);
int Load_OBJECT4DV2_COB(OBJECT4DV2_PTR obj,	// pointer to object
						char *filename,		   // filename of Caligari COB file
						VECTOR4D_PTR scale,	// initial scaling factors
						VECTOR4D_PTR pos,	  // initial position
						VECTOR4D_PTR rot,	  // initial rotations
						int vertex_flags = 0); // flags to re-order vertices
