/* Copyright (C) 2006 Imperial College London and others.

 Please see the AUTHORS file in the main source directory for a full list
 of copyright holders.

 Dr Gerard J Gorman
 Applied Modelling and Computation Group
 Department of Earth Science and Engineering
 Imperial College London

 g.gorman@imperial.ac.uk

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
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
#ifndef DEBUG_H
#define DEBUG_H

#include "confdefs.h"
#include <string>
#include <iostream>

#ifndef NDEBUG
#define ECHO(X)    std::cerr <<"("<<__FILE__<<","<<__LINE__<<"): "<<X<<std::endl
#define CHECK(X)   std::cerr <<"("<<__FILE__<<","<<__LINE__<<"): "<<#X<<" = "<<X<<std::endl
#define WARNING(X) std::cerr <<"("<<__FILE__<<","<<__LINE__<<"): WARNING - "<<X<<std::endl
#else 
#define ECHO(X)
#define CHECK(X)
#define WARNING(X)
#endif

#define ERROR(X) std::cerr <<"("<<__FILE__<<","<<__LINE__<<"): ERROR - "<<X<<std::endl

std::string mpiErrorMessage(const int err);

#endif
