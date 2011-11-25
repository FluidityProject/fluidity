static char adSid[]="$Id: adBuffer.c $";

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "adStack.h"

/************ MEASUREMENT OF PUSH/POP TRAFFIC *************/

static int mmftraffic = 0 ;
static int mmftrafficM = 0 ;

void addftraffic(int n) {
  mmftraffic = mmftraffic+n ;
  while (mmftraffic >= 1000000) {
    mmftraffic = mmftraffic-1000000 ;
    ++mmftrafficM ;
  }
  while (mmftraffic < 0) {
    mmftraffic = mmftraffic+1000000 ;
    --mmftraffic ;
  }
}

void printtraffic() {
  printctraffic_() ;
  printf(" F Traffic: ") ;
  printbigbytes(mmftrafficM, 1000000, mmftraffic) ;
  printf(" bytes\n") ;
}

/************************** integer*4 ************************/
static int adi4buf[512] ;
static int adi4ibuf = 0 ;
static int adi4lbuf[512] ;
static int adi4ilbuf = -1 ;
static int adi4inlbuf = 0 ;

void pushinteger4(int x) {
  addftraffic(4) ;
  if (adi4ilbuf != -1) {
    adi4ilbuf = -1 ;
    adi4inlbuf = 0 ;
  }
  if (adi4ibuf >= 511) {
    adi4buf[511] = x ;
    pushNarray((char*)adi4buf, 512*4) ;
    addftraffic(-512*4) ;
    adi4ibuf = 0 ;
  } else {
    adi4buf[adi4ibuf] = x ;
    ++adi4ibuf ;
  }
}

void lookinteger4(int *x) {
  if (adi4ilbuf == -1) {
    adi4ilbuf = adi4ibuf ;
    resetadlookstack_() ;
  }
  if (adi4ilbuf <= 0) {
    lookNarray((char*)adi4lbuf, 512*4) ;
    adi4inlbuf = 1 ;
    adi4ilbuf = 511 ;
    *x = adi4lbuf[511] ;
  } else {
    --adi4ilbuf ;
    if (adi4inlbuf)
      *x = adi4lbuf[adi4ilbuf] ;
    else
      *x = adi4buf[adi4ilbuf] ;
  }
}

void popinteger4(int *x) {
  if (adi4ilbuf != -1) {
    adi4ilbuf = -1 ;
    adi4inlbuf = 0 ;
  }
  if (adi4ibuf <= 0) {
    popNarray((char*)adi4buf, 512*4) ;
    adi4ibuf = 511 ;
    *x = adi4buf[511] ;
  } else {
    --adi4ibuf ;
    *x = adi4buf[adi4ibuf] ;
  }
}

/*************************** bits *************************/
static unsigned int adbitbuf = 0 ;
static int adbitibuf = 1 ;
static unsigned int adbitlbuf = 0 ;
static int adbitilbuf = -1 ;
static int adbitinlbuf = 0 ;

void pushbit(int bit) {
  if (adbitilbuf != -1) {
    adbitilbuf = -1 ;
    adbitinlbuf = 0 ;
  }
  adbitbuf<<=1 ;
  if (bit) adbitbuf++ ;
  if (adbitibuf>=32) {
    pushinteger4(adbitbuf) ;
    adbitbuf = 0 ;
    adbitibuf = 1 ;
  } else
    adbitibuf++ ;
}

int lookbit() {
  if (adbitilbuf==-1) {
    adbitilbuf=adbitibuf ;
    adbitlbuf = adbitbuf ;
  }
  if (adbitilbuf<=1) {
    lookinteger4(&adbitlbuf) ;
    adbitilbuf = 32 ;
  } else
    adbitilbuf-- ;
  int bit = adbitlbuf%2 ;
  adbitlbuf>>=1 ;
  return bit ;
}

int popbit() {
  if (adbitilbuf != -1) {
    adbitilbuf = -1 ;
    adbitinlbuf = 0 ;
  }
  if (adbitibuf<=1) {
    popinteger4(&adbitbuf) ;
    adbitibuf = 32 ;
  } else
    adbitibuf-- ;
  int bit = adbitbuf%2 ;
  adbitbuf>>=1 ;
  return bit ;
}

/************************* controls ***********************/

void pushcontrol1b(int cc) {
  pushbit(cc) ;
}

void popcontrol1b(int *cc) {
  *cc = (popbit()?1:0) ;
}

void lookcontrol1b(int *cc) {
  *cc = (lookbit()?1:0) ;
}

void pushcontrol2b(int cc) {
  pushbit(cc%2) ;
  cc >>= 1 ;
  pushbit(cc) ;
}

void popcontrol2b(int *cc) {
  *cc = (popbit()?2:0) ;
  if (popbit()) (*cc)++ ;
}

void lookcontrol2b(int *cc) {
  *cc = (lookbit()?2:0) ;
  if (lookbit()) (*cc)++ ;
}

void pushcontrol3b(int cc) {
  pushbit(cc%2) ;
  cc >>= 1 ;
  pushbit(cc%2) ;
  cc >>= 1 ;
  pushbit(cc) ;
}

void popcontrol3b(int *cc) {
  *cc = (popbit()?2:0) ;
  if (popbit()) (*cc)++ ;
  (*cc) <<= 1 ;
  if (popbit()) (*cc)++ ;
}

void lookcontrol3b(int *cc) {
  *cc = (lookbit()?2:0) ;
  if (lookbit()) (*cc)++ ;
  (*cc) <<= 1 ;
  if (lookbit()) (*cc)++ ;
}

void pushcontrol4b(int cc) {
  pushbit(cc%2) ;
  cc >>= 1 ;
  pushbit(cc%2) ;
  cc >>= 1 ;
  pushbit(cc%2) ;
  cc >>= 1 ;
  pushbit(cc) ;
}

void popcontrol4b(int *cc) {
  *cc = (popbit()?2:0) ;
  if (popbit()) (*cc)++ ;
  (*cc) <<= 1 ;
  if (popbit()) (*cc)++ ;
  (*cc) <<= 1 ;
  if (popbit()) (*cc)++ ;
}

void lookcontrol4b(int *cc) {
  *cc = (lookbit()?2:0) ;
  if (lookbit()) (*cc)++ ;
  (*cc) <<= 1 ;
  if (lookbit()) (*cc)++ ;
  (*cc) <<= 1 ;
  if (lookbit()) (*cc)++ ;
}

void pushcontrol5b(int cc) {
  pushbit(cc%2) ;
  cc >>= 1 ;
  pushbit(cc%2) ;
  cc >>= 1 ;
  pushbit(cc%2) ;
  cc >>= 1 ;
  pushbit(cc%2) ;
  cc >>= 1 ;
  pushbit(cc) ;
}

void popcontrol5b(int *cc) {
  *cc = (popbit()?2:0) ;
  if (popbit()) (*cc)++ ;
  (*cc) <<= 1 ;
  if (popbit()) (*cc)++ ;
  (*cc) <<= 1 ;
  if (popbit()) (*cc)++ ;
  (*cc) <<= 1 ;
  if (popbit()) (*cc)++ ;
}

void lookcontrol5b(int *cc) {
  *cc = (lookbit()?2:0) ;
  if (lookbit()) (*cc)++ ;
  (*cc) <<= 1 ;
  if (lookbit()) (*cc)++ ;
  (*cc) <<= 1 ;
  if (lookbit()) (*cc)++ ;
  (*cc) <<= 1 ;
  if (lookbit()) (*cc)++ ;
}

/************************** real*4 ************************/
static float adr4buf[512] ;
static int adr4ibuf = 0 ;
static float adr4lbuf[512] ;
static int adr4ilbuf = -1 ;
static int adr4inlbuf = 0 ;

void pushreal4(float x) {
  addftraffic(4) ;
  if (adr4ilbuf != -1) {
    adr4ilbuf = -1 ;
    adr4inlbuf = 0 ;
  }
  if (adr4ibuf >= 511) {
    adr4buf[511] = x ;
    pushNarray((char*)adr4buf, 512*4) ;
    addftraffic(-512*4) ;
    adr4ibuf = 0 ;
  } else {
    adr4buf[adr4ibuf] = x ;
    ++adr4ibuf ;
  }
}

void lookreal4(float *x) {
  if (adr4ilbuf == -1) {
    adr4ilbuf = adr4ibuf ;
    resetadlookstack_() ;
  }
  if (adr4ilbuf <= 0) {
    lookNarray((char*)adr4lbuf, 512*4) ;
    adr4inlbuf = 1 ;
    adr4ilbuf = 511 ;
    *x = adr4lbuf[511] ;
  } else {
    --adr4ilbuf ;
    if (adr4inlbuf)
      *x = adr4lbuf[adr4ilbuf] ;
    else
      *x = adr4buf[adr4ilbuf] ;
  }
}

void popreal4(float *x) {
  if (adr4ilbuf != -1) {
    adr4ilbuf = -1 ;
    adr4inlbuf = 0 ;
  }
  if (adr4ibuf <= 0) {
    popNarray((char*)adr4buf, 512*4) ;
    adr4ibuf = 511 ;
    *x = adr4buf[511] ;
  } else {
    --adr4ibuf ;
    *x = adr4buf[adr4ibuf] ;
  }
}

/************************** real*8 ************************/
static double adr8buf[512] ;
static int adr8ibuf = 0 ;
static double adr8lbuf[512] ;
static int adr8ilbuf = -1 ;
static int adr8inlbuf = 0 ;

void pushreal8(double x) {
  addftraffic(8) ;
  if (adr8ilbuf != -1) {
    adr8ilbuf = -1 ;
    adr8inlbuf = 0 ;
  }
  if (adr8ibuf >= 511) {
    adr8buf[511] = x ;
    pushNarray((char*)adr8buf, 512*8) ;
    addftraffic(-4096) ;
    adr8ibuf = 0 ;
  } else {
    adr8buf[adr8ibuf] = x ;
    ++adr8ibuf ;
  }
}

void lookreal8(double *x) {
  if (adr8ilbuf == -1) {
    adr8ilbuf = adr8ibuf ;
    resetadlookstack_() ;
  }
  if (adr8ilbuf <= 0) {
    lookNarray((char*)adr8lbuf, 512*8) ;
    adr8inlbuf = 1 ;
    adr8ilbuf = 511 ;
    *x = adr8lbuf[511] ;
  } else {
    --adr8ilbuf ;
    if (adr8inlbuf)
      *x = adr8lbuf[adr8ilbuf] ;
    else
      *x = adr8buf[adr8ilbuf] ;
  }
}

void popreal8(double *x) {
  if (adr8ilbuf != -1) {
    adr8ilbuf = -1 ;
    adr8inlbuf = 0 ;
  }
  if (adr8ibuf <= 0) {
    popNarray((char*)adr8buf, 512*8) ;
    adr8ibuf = 511 ;
    *x = adr8buf[511] ;
  } else {
    --adr8ibuf ;
    *x = adr8buf[adr8ibuf] ;
  }
}

/********* PRINTING THE SIZE OF STACKS AND BUFFERS ********/

void printbuffertop() {
  int size = 0 ;
  size += adi4ibuf*4 ;
  size += adr4ibuf*4 ;
  size += adr8ibuf*8 ;
  printf("Buffer size:%i bytes i.e. %i Kbytes\n",
         size, (int) (size/1024.0)) ;
}

void showallstacks() {
  int i ;
  printf("BIT STACK      : %x == %i\n",adbitbuf,adbitbuf) ; 
  printf("INTEGER*4 BUFFER[%i]:",adi4ibuf) ;
  for (i=0 ; i<adi4ibuf ; ++i) printf(" %i",adi4buf[i]) ;
  printf("\n") ;
  printf("REAL*8 BUFFER[%i]:",adr8ibuf) ;
  for (i=0 ; i<adr8ibuf ; ++i) printf(" %d",(int) adr8buf[i]) ;
  printf("\n") ;
  printf("REAL*4 BUFFER[%i]:",adr4ibuf) ;
  for (i=0 ; i<adr4ibuf ; ++i) printf(" %f",adr4buf[i]) ;
  printf("\n") ;
  showrecentcstack_() ;
}

/**********************************************************
 *        HOW TO CREATE PUSH* LOOK* POP* SUBROUTINES
 *              YET FOR OTHER DATA TYPES
 * Duplicate and uncomment the commented code below.
 * In the duplicated code, replace:
 *   ctct -> C type name (e.g. float double, int...)
 *   TTTT -> BASIC TAPENADE TYPE NAME
 *     (in character, boolean, integer, real, complex, pointer,...)
 *   z7   -> LETTER-SIZE FOR TYPE
 *     (in s,         b,       i,       r,    c,       p,      ...)
 *   7    -> TYPE SIZE IN BYTES
 * Don't forget to insert the corresponding lines in
 * procedure printbuffertop(), otherwise the contribution of
 * this new type to buffer occupation will not be seen.
 * (not very important anyway...)
 **********************************************************/

/************************** TTTT*7 ************************/
/*
static ctct adz7buf[512] ;
static int adz7ibuf = 0 ;
static ctct adz7lbuf[512] ;
static int adz7ilbuf = -1 ;
static int adz7inlbuf = 0 ;

void pushTTTT7(ctct x) {
  addftraffic(7) ;
  if (adz7ilbuf != -1) {
    adz7ilbuf = -1 ;
    adz7inlbuf = 0 ;
  }
  if (adz7ibuf >= 511) {
    adz7buf[511] = x ;
    pushNarray((char*)adz7buf, 512*7) ;
    addftraffic(-512*7) ;
    adz7ibuf = 0 ;
  } else {
    adz7buf[adz7ibuf] = x ;
    ++adz7ibuf ;
  }
}

void lookTTTT7(ctct *x) {
  if (adz7ilbuf == -1) {
    adz7ilbuf = adz7ibuf ;
    resetadlookstack_() ;
  }
  if (adz7ilbuf <= 0) {
    lookNarray((char*)adz7lbuf, 512*7) ;
    adz7inlbuf = 1 ;
    adz7ilbuf = 511 ;
    *x = adz7lbuf[511] ;
  } else {
    --adz7ilbuf ;
    if (adz7inlbuf)
      *x = adz7lbuf[adz7ilbuf] ;
    else
      *x = adz7buf[adz7ilbuf] ;
  }
}

void popTTTT7(ctct *x) {
  if (adz7ilbuf != -1) {
    adz7ilbuf = -1 ;
    adz7inlbuf = 0 ;
  }
  if (adz7ibuf <= 0) {
    popNarray((char*)adz7buf, 512*7) ;
    adz7ibuf = 511 ;
    *x = adz7buf[511] ;
  } else {
    --adz7ibuf ;
    *x = adz7buf[adz7ibuf] ;
  }
}
*/

/**********************************************************/
