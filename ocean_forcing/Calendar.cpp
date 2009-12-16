/*  Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Dr Gerard Gorman
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    g.gorman@imperial.ac.uk
    
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

#include "Calendar.h" 

using namespace std;

Calendar::Calendar(){
#ifdef HAVE_LIBUDUNITS
  if(!utIsInit()){
    // The definitions in the units file are read into memory.
    int err = utInit(NULL);
    if(err==UT_ENOFILE){
      cerr<<"ERROR ("<<__FILE__<<"): The units file doesn’t exist.\n";
    }else if(err==UT_ESYNTAX){
      cerr<<"ERROR ("<<__FILE__<<"): The units file contains a syntax  error.\n";
    }else if(err==UT_EUNKNOWN){
      cerr<<"ERROR ("<<__FILE__<<"): The units  file  contains  an  unknown specification.\n";
    }else if(err==UT_EIO){
      cerr<<"ERROR ("<<__FILE__<<"): An I/O error occurred while accessing the units file.\n";
    }else if(err==UT_EALLOC){
      cerr<<"ERROR ("<<__FILE__<<"): A memory allocation failure occurred.\n";
    }
    if(err)
      exit(err);
  }
#endif
  SetTransformation("seconds since 1-1-1 0:0:0", 
                    "seconds since 1-1-1 0:0:0",
                    "none");
}

Calendar::Calendar(const Calendar& in){
  *this = in;
}

const Calendar& Calendar::operator=(const Calendar &in){
#ifdef HAVE_LIBUDUNITS
  from_unit = in.from_unit;
  to_unit = in.to_unit;
#endif
  calendar = in.calendar;
  
  return *this;
}

Calendar::~Calendar(){}

int Calendar::SetTransformation(string _from_unit, string _to_unit, string _calendar){
  int err=0;

  // Reserved for future use...
  // calendar = _calendar;

#ifdef HAVE_LIBUDUNITS
  err = utScan(_from_unit.c_str(), &from_unit);
  if(err==UT_ENOINIT){
    cerr<<"ERROR ("<<__FILE__<<"): UDUNITS has not been initialized.\n";
  }else if(err==UT_EINVALID){
    cerr<<"ERROR ("<<__FILE__<<"): The unit argument is a null pointer.\n";
  }else if(err==UT_EUNKNOWN){
    cerr<<"ERROR ("<<__FILE__<<"): The specification contains an unknown unit.\n";
  }else if(err==UT_ESYNTAX){
    cerr<<"ERROR ("<<__FILE__<<"): The specification contains a syntax error.\n";
  }
  if(err)
    exit(-1);

  err = utScan(_to_unit.c_str(), &to_unit);
  if(err==UT_ENOINIT){
    cerr<<"ERROR ("<<__FILE__<<"): UDUNITS has not been initialized.\n";
  }else if(err==UT_EINVALID){
    cerr<<"ERROR ("<<__FILE__<<"): The unit argument is a null pointer.\n";
  }else if(err==UT_EUNKNOWN){
    cerr<<"ERROR ("<<__FILE__<<"): The specification contains an unknown unit.\n";
  }else if(err==UT_ESYNTAX){
    cerr<<"ERROR ("<<__FILE__<<"): The specification contains a syntax error.\n";
  }
  if(err)
    exit(-1);
  
  err = utConvert(&from_unit, &to_unit, &slope, &intercept);
  if(err==UT_ENOINIT){
    cerr<<"ERROR ("<<__FILE__<<"): The package has not been initialized.\n";
  }else if(err==UT_EINVALID){
    cerr<<"ERROR ("<<__FILE__<<"): The unit structures is invalid.\n";
  }else if(err==UT_ECONVERT){
    cerr<<"ERROR ("<<__FILE__<<"): The units are not convertible.\n"
        <<"from_unit = "<<_from_unit<<endl
        <<"to_unit = "<<_to_unit<<endl;
  }
  if(err)
    exit(-1);

#endif
  return err;
}

int Calendar::Convert(double from_value, double& to_value){
  to_value = slope*from_value + intercept;
  
  return 0;
}
