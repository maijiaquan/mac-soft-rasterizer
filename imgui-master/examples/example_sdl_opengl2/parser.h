#pragma once
#include "predefine.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include<stdlib.h>

// parser class ///////////////////////////////////////////////
class CPARSERV1
{
public:
	// constructor /////////////////////////////////////////////////
	CPARSERV1();

	// destructor ///////////////////////////////////////////////////
	~CPARSERV1();

	// reset file system ////////////////////////////////////////////
	int Reset();

	// open file /////////////////////////////////////////////////////
	int Open(char *filename);

	// close file ////////////////////////////////////////////////////
	int Close();

	// get line //////////////////////////////////////////////////////
	char *Getline(int mode);

	// sets the comment string ///////////////////////////////////////
	int SetComment(char *string);

	// find pattern in line //////////////////////////////////////////
	int Pattern_Match(char *string, char *pattern, ...);

	// VARIABLE DECLARATIONS /////////////////////////////////////////

public:
	FILE *fstream;					  // file pointer
	char buffer[PARSER_BUFFER_SIZE];  // line buffer
	int length;						  // length of current line
	int num_lines;					  // number of lines processed
	char comment[PARSER_MAX_COMMENT]; // single line comment string

	// pattern matching parameter storage, easier that variable arguments
	// anything matched will be stored here on exit from the call to pattern()
	char pstrings[PATTERN_MAX_ARGS][PATTERN_BUFFER_SIZE]; // any strings
	int num_pstrings;

	float pfloats[PATTERN_MAX_ARGS]; // any floats
	int num_pfloats;

	int pints[PATTERN_MAX_ARGS]; // any ints
	int num_pints;

}; // end CLASS CPARSERV1 //////////////////////////////////////////////



int ReplaceChars(char *string_in, char *string_out, char *replace_chars, char rep_char, int case_on = 1);
char *StringLtrim(char *string);
char *StringRtrim(char *string);
float IsFloat(char *fstring);
int IsInt(char *istring);