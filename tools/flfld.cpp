/*  Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    C.Pain@Imperial.ac.uk
    
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation,
    version 2.1 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
    USA
*/

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "fmangle.h"

using namespace std;

// keep this filehandle scoped to this file
static ifstream fldfile;

extern "C" {
 void fl_fldopen_fc(char *, int *, int *);
 void fl_fldread_fc(char *, int *, char *, int *, int *, int *, int *);
 void fl_fldclose_fc();
}

//
// Convert word to uppercase.
//
void uppercase(string& word){
  for(string::iterator c = word.begin(); c != word.end(); ++c){
    if( (*c >= 'a')&&(*c <= 'z') ){
      *c = (*c - 'a') + 'A';
    }
  }
}

//
// Split a line up into tokens
//
void Tokenize(const string& str, vector<string>& tokens, const string& delimiters = " "){
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}

void fl_fldopen_fc(char *ffile, int *ffile_len, int *err){
  string filename(ffile, *ffile_len);
  fldfile.open(filename.c_str(), ifstream::in);
  
  if( fldfile.fail() )
    *err = -1;
  else
    *err = 0;
  
  return;
}
void fl_fldread_fc(char *fld1, int *fld1_len, char *fld2, int *fld2_len,
		 int *values, int *nvalues, int *err){
  string line;
  vector<string> tokens;
  
  *err = 0;

  // Check if there is something to read.
  if( !fldfile.is_open() ){
    *err = -1;
    return;
  }
  
  // Try to get a line ignoring bacnk lines and comments 
  for(;;){
    if( getline(fldfile, line) ){
      Tokenize(line, tokens, " ");
      unsigned tcnt = tokens.size();
      
      // If this is a blank line then move to the next one
      if( tcnt < 1) 
	continue;

      // Skip if this is a comment
      if( (tokens[0][0] == '#') || (tokens[0][0] == '@') )
	continue;
      
      // Quit if this is a foobar line
      if( tcnt < 3 ){
	*err = -2;
	return;
      }
      
      // Sort out the line
      // scalar/vector/tensor
      *fld1_len = tokens[0].size();
      uppercase( tokens[0] );
      strcpy(fld1, tokens[0].c_str());
      
      *fld2_len = tokens[1].size();
      strcpy(fld2, tokens[1].c_str());
      
      *nvalues = tcnt - 2;
      unsigned vcnt = 0;
      for(unsigned i=2; i<tcnt; i++)
	values[vcnt++] = atoi( tokens[i].c_str() );
      
      // We have what we want so lets get out of here.
      break;
    }else{
      *err = -3;
      return;
    }
  }
  
  return;
}

void fl_fldclose_fc(){
  fldfile.close();
}










