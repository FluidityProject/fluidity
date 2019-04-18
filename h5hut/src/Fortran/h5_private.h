/*
  Copyright (c) 2006-2016, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#ifndef __FORTRAN_H5_PRIVATE_H
#define __FORTRAN_H5_PRIVATE_H

#include <stdlib.h>
#include <string.h>

#include "h5core/h5_types.h"

static inline char*
h5_strdupfor2c (
	const char* s,
	const ssize_t len
	) {
        // :FIXME: error handling
	char* dup = (char*)malloc (len + 1);
	strncpy (dup, s, len);
	dup[len] = '\0';
	for (int i = len-1; i >= 0; i--) {
		if (dup[i] == ' ') dup[i] = '\0';
		else break;
	}
	return dup;
}

static inline char*
h5_strc2for (
	char* const str,
	const ssize_t l_str
	) {
	size_t len = strlen (str);
	memset (str+len, ' ', l_str-len);

	return str;
}

static inline h5_file_t
h5_filehandlefor2c (
	const h5_int64_t* ptr
	) {
	return (h5_file_t)*ptr;
}

static inline
int strlenf (
        const char* s,
        int len
        ) {
        if (len == 0) return 0;
        while (s[--len] == ' ');
        return ++len;
}

#endif
