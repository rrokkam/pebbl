// This file is massive KLUDGE just to see if we can get a compile

#define ACRO_VALIDATING 1


#ifndef _INCLUDE_ACRO_CONFIG_H
#define _INCLUDE_ACRO_CONFIG_H 1
 
#ifndef ACRO_BUILD_CPU_X86 
#define ACRO_BUILD_CPU_X86  1 
#endif

#define ACRO_BUILD_LINUX 1

/* define if the compiler supports exceptions */
#ifndef ACRO_HAVE_EXCEPTIONS 
#define ACRO_HAVE_EXCEPTIONS   
#endif

/* define if the compiler supports the explicit keyword */
#ifndef ACRO_HAVE_EXPLICIT 
#define ACRO_HAVE_EXPLICIT   
#endif

/* Define to 1 if you have the `getrusage' function. */
/* #ifndef ACRO_HAVE_GETRUSAGE  */
/* #define ACRO_HAVE_GETRUSAGE  1  */
/* #endif */

/* Define to 1 if you have the <inttypes.h> header file. */
// #ifndef ACRO_HAVE_INTTYPES_H 
// #define ACRO_HAVE_INTTYPES_H  1 
// #endif

/* Define if you have LAPACK library. */
/* #undef ACRO_HAVE_LAPACK */

/* define if the compiler supports member templates */
#ifndef ACRO_HAVE_MEMBER_TEMPLATES 
#define ACRO_HAVE_MEMBER_TEMPLATES   
#endif

/* Define to 1 if you have the <memory.h> header file. */
#ifndef ACRO_HAVE_MEMORY_H 
#define ACRO_HAVE_MEMORY_H  1 
#endif

/* define that mpi is being used */
#define ACRO_HAVE_MPI 1
#define UTILIB_HAVE_MPI 1

/* define if the compiler implements namespaces */
#ifndef ACRO_HAVE_NAMESPACES 
#define ACRO_HAVE_NAMESPACES   
#endif

/* define if the compiler has stringstream */
#ifndef ACRO_HAVE_SSTREAM 
#define ACRO_HAVE_SSTREAM   
#endif

/* define if the compiler supports ISO C++ standard library */
#ifndef ACRO_HAVE_STD 
#define ACRO_HAVE_STD   
#endif

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef ACRO_HAVE_STDINT_H 
#define ACRO_HAVE_STDINT_H  1 
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef ACRO_HAVE_STDLIB_H 
#define ACRO_HAVE_STDLIB_H  1 
#endif

/* Define to 1 if you have the `strerror' function. */
#ifndef ACRO_HAVE_STRERROR 
#define ACRO_HAVE_STRERROR  1 
#endif

/* Define to 1 if you have the <strings.h> header file. */
/* #ifndef ACRO_HAVE_STRINGS_H  */
/* #define ACRO_HAVE_STRINGS_H  1  */
/* #endif */

/* Define to 1 if you have the <string.h> header file. */
#ifndef ACRO_HAVE_STRING_H 
#define ACRO_HAVE_STRING_H  1 
#endif

/* Define to 1 if you have the <sys/stat.h> header file. */
/* #ifndef ACRO_HAVE_SYS_STAT_H  */
/* #define ACRO_HAVE_SYS_STAT_H  1  */
/* #endif */

/* Define to 1 if you have the <sys/types.h> header file. */
/* #ifndef ACRO_HAVE_SYS_TYPES_H  */
/* #define ACRO_HAVE_SYS_TYPES_H  1  */
/* #endif */

/* Define to 1 if you have the <unistd.h> header file. */
/* #ifndef ACRO_HAVE_UNISTD_H  */
/* #define ACRO_HAVE_UNISTD_H  1  */
/* #endif */

/* Define to 1 if you have the <values.h> header file. */
/* #undef ACRO_HAVE_VALUES_H */

/* software host will be cygwin */
// #ifndef ACRO_HOST_CYGWIN 
// #define ACRO_HOST_CYGWIN  1 
// #endif

/* software host is GNU */
/* #undef ACRO_HOST_GNU */

/* software host will be linux */
/* #undef ACRO_HOST_LINUX */

/* software host will be mingw */
/* #undef ACRO_HOST_MINGW */

/* software host will be solaris */
/* #undef ACRO_HOST_SOLARIS */

/* Name of package */
#ifndef ACRO_PACKAGE 
#define ACRO_PACKAGE  "acro" 
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef ACRO_PACKAGE_BUGREPORT 
#define ACRO_PACKAGE_BUGREPORT  "acro-help@software.sandia.gov" 
#endif

/* Define to the full name of this package. */
#ifndef ACRO_PACKAGE_NAME 
#define ACRO_PACKAGE_NAME  "acro" 
#endif

/* Define to the full name and version of this package. */
#ifndef ACRO_PACKAGE_STRING 
#define ACRO_PACKAGE_STRING  "acro VOTD" 
#endif

/* Define to the one symbol short name of this package. */
#ifndef ACRO_PACKAGE_TARNAME 
#define ACRO_PACKAGE_TARNAME  "acro" 
#endif

/* Define to the version of this package. */
#ifndef ACRO_PACKAGE_VERSION 
#define ACRO_PACKAGE_VERSION  "VOTD" 
#endif

/* Define to 1 if you have the ANSI C header files. */
/* #ifndef ACRO_STDC_HEADERS  */
/* #define ACRO_STDC_HEADERS  1  */
/* #endif */

/* software target will be cygwin */
// #ifndef ACRO_TARGET_CYGWIN 
// #define ACRO_TARGET_CYGWIN  1 
// #endif

/* software target will be linux */
/* #undef ACRO_TARGET_LINUX */

/* software target will be mingw */
/* #undef ACRO_TARGET_MINGW */

/* software target will be solaris */
/* #undef ACRO_TARGET_SOLARIS */

/* Define if want to build with ampl enabled */
/* #undef ACRO_USING_AMPL */

/* Define if want to build with appspack enabled */
/* #undef ACRO_USING_APPSPACK */

/* define that clp is being used */
// #ifndef ACRO_USING_CLP 
// #define ACRO_USING_CLP  1 
// #endif

/* define that cobyla is being used */
// #ifndef ACRO_USING_COBYLA 
// #define ACRO_USING_COBYLA  1 
// #endif

/* Define if want to build with coin enabled */
/* #undef ACRO_USING_COIN */

/* Define if want to build with colin enabled */
/* #undef ACRO_USING_COLIN */

/* Define if want to build with coliny enabled */
/* #undef ACRO_USING_COLINY */

/* Define if want to build with dscpack enabled */
/* #undef ACRO_USING_DSCPACK */

/* Define if want to build with exact enabled */
// #ifndef ACRO_USING_EXACT 
// #define ACRO_USING_EXACT   
// #endif

/* Define if want to build with filib enabled */
/* #undef ACRO_USING_FILIB */

/* Define if want to build with glpk enabled */
/* #undef ACRO_USING_GLPK */

/* Define if want to build with gnlp enabled */
/* #undef ACRO_USING_GNLP */

/* Define if want to build with ipopt enabled */
/* #undef ACRO_USING_IPOPT */

/* Define if want to build with mtl enabled */
/* #undef ACRO_USING_MTL */

/* Define if want to build with optpp enabled */
/* #undef ACRO_USING_OPTPP */

/* Define if want to build with parpcx enabled */
/* #undef ACRO_USING_PARPCX */

/* Define if want to build with pebbl enabled */
/* #undef ACRO_USING_PEBBL */

/* Define if want to build with pico enabled */
/* #undef ACRO_USING_PICO */

/* define that plgo is being used */
// #ifndef ACRO_USING_PLGO 
// #define ACRO_USING_PLGO  1 
// #endif

/* Define if want to build with soplex enabled */
/* #undef ACRO_USING_SOPLEX */

/* Define if want to build with 3po enabled */
/* #undef ACRO_USING_THREEPO */

/* Define if want to build with tmf enabled */
/* #undef ACRO_USING_TMF */

/* Define if want to build with tracecache enabled */
/* #undef ACRO_USING_TRACECACHE */

/* Define if want to build with trilinos enabled */
/* #undef ACRO_USING_TRILINOS */

/* Define if want to build with utilib enabled */
#ifndef ACRO_USING_UTILIB 
#define ACRO_USING_UTILIB   
#endif

/* turn on code validation tests */
/* #undef ACRO_VALIDATING */

/* Version number of package */
// #ifndef ACRO_VERSION 
// #define ACRO_VERSION  "VOTD" 
// #endif

/* Define to 1 if your processor stores words with the most significant byte
   first (like Motorola and SPARC, unlike Intel and VAX). */
/* #undef ACRO_WORDS_BIGENDIAN */

/* define whether checksum function is included in utilib */
// #ifndef ACRO_YES_CHECKSUM 
// #define ACRO_YES_CHECKSUM   
// #endif

/* define whether CommonIO is included in utilib */
#ifndef ACRO_YES_COMMONIO 
#define ACRO_YES_COMMONIO   
#endif

/* define whether memdebug is included in utilib */
/* #undef ACRO_YES_MEMDEBUG */

/* Define to `unsigned' if <sys/types.h> does not define. */
/* #undef _acro_size_t */

// #ifdef ACRO_USING_UTILIB
// #include <pebbl/utilib/utilib_config.h>
// #include <acro_config_bool.h>
// #include <acro_config_explicit.h>
// #endif

 
/* once: _INCLUDE_ACRO_CONFIG_H */
#endif


#ifndef _UTILIB_UTILIB_CONFIG_H
#define _UTILIB_UTILIB_CONFIG_H 1
 
/* utilib/utilib_config.h. Generated automatically at end of configure. */
/* utilib/config.h.  Generated by configure.  */
/* utilib/config.h.in.  Generated from configure.ac by autoheader.  */


/*  _________________________________________________________________________
 *
 *  UTILIB: A utility library for developing portable C++ codes.
 *  Copyright (c) 2008 Sandia Corporation.
 *  This software is distributed under the BSD License.
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *  For more information, see the README file in the top UTILIB directory.
 *  _________________________________________________________________________
 */

/**
 * \file utilib_config.h
 *
 * A header file for UTILIB configuration options that is generated by
 * autoconf.
 */


/* software build cpu is sparc */
/* #undef UTILIB_BUILD_CPU_SPARC */

/* software build cpu is x86 */
#ifndef UTILIB_BUILD_CPU_X86 
#define UTILIB_BUILD_CPU_X86  1 
#endif

/* software build cpu is 64 bit x86 */
/* #undef UTILIB_BUILD_CPU_X86_64 */

/* software build os is cygwin */
#ifndef UTILIB_BUILD_CYGWIN 
#define UTILIB_BUILD_CYGWIN  1 
#endif

/* software build os is linux */
/* #undef UTILIB_BUILD_LINUX */

/* software build os is solaris */
/* #undef UTILIB_BUILD_SOLARIS */

/* Define to 1 if you have the `atexit' function. */
#ifndef UTILIB_HAVE_ATEXIT 
#define UTILIB_HAVE_ATEXIT  1 
#endif

/* Define to 1 if you have the `clock' function. */
#define UTILIB_HAVE_CLOCK 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#ifndef UTILIB_HAVE_DLFCN_H 
#define UTILIB_HAVE_DLFCN_H  1 
#endif

/* define if the compiler supports exceptions */
#ifndef UTILIB_HAVE_EXCEPTIONS 
#define UTILIB_HAVE_EXCEPTIONS   
#endif

/* define if the compiler supports the explicit keyword */
#ifndef UTILIB_HAVE_EXPLICIT 
#define UTILIB_HAVE_EXPLICIT   
#endif

/* Define to 1 if you have the <float.h> header file. */
#ifndef UTILIB_HAVE_FLOAT_H 
#define UTILIB_HAVE_FLOAT_H  1 
#endif

/* Define to 1 if you have the `ftime' function. */
/* #undef UTILIB_HAVE_FTIME */

/* Define to 1 if you have the `getcwd' function. */
#ifndef UTILIB_HAVE_GETCWD 
#define UTILIB_HAVE_GETCWD  1 
#endif

/* Define to 1 if you have the `getrusage' function. */
#ifndef UTILIB_HAVE_GETRUSAGE 
#define UTILIB_HAVE_GETRUSAGE  1 
#endif 

/* Define to 1 if you have the `gettimeofday' function. */
#ifndef UTILIB_HAVE_GETTIMEOFDAY
#define UTILIB_HAVE_GETTIMEOFDAY  1
#endif

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef UTILIB_HAVE_INTTYPES_H 
#define UTILIB_HAVE_INTTYPES_H  1 
#endif

/* Define to 1 if you have the `m' library (-lm). */
#ifndef UTILIB_HAVE_LIBM 
#define UTILIB_HAVE_LIBM  1 
#endif

/* Define to 1 if you have the <limits.h> header file. */
#ifndef UTILIB_HAVE_LIMITS_H 
#define UTILIB_HAVE_LIMITS_H  1 
#endif

/* Define to 1 if you have the `localtime' function. */
#ifndef UTILIB_HAVE_LOCALTIME
#define UTILIB_HAVE_LOCALTIME 1
#endif

/* Define to 1 if long double works and has more range or precision than
   double. */
#ifndef UTILIB_HAVE_LONG_DOUBLE 
#define UTILIB_HAVE_LONG_DOUBLE  1 
#endif

#if 0
/* define whether libm contains lround() */
#ifndef UTILIB_HAVE_LROUND 
#define UTILIB_HAVE_LROUND   
#endif
#endif

/* define if the compiler supports member templates */
#ifndef UTILIB_HAVE_MEMBER_TEMPLATES 
#define UTILIB_HAVE_MEMBER_TEMPLATES   
#endif

#define UTILIB_NO_MEMBER_TEMPLATE_FRIENDS 1

/* Define to 1 if you have the <memory.h> header file. */
#ifndef UTILIB_HAVE_MEMORY_H 
#define UTILIB_HAVE_MEMORY_H  1 
#endif

/* define that mpi is being used */
/* #undef UTILIB_HAVE_MPI */

/* define if the compiler implements namespaces */
#ifndef UTILIB_HAVE_NAMESPACES
#define UTILIB_HAVE_NAMESPACES
#endif

/* Define to 1 if you have the `nrand48' function. */
/* #ifndef UTILIB_HAVE_NRAND48 */
/* #define UTILIB_HAVE_NRAND48  1 */
/* #endif */

/* Define to 1 if the system has the type `ptrdiff_t'. */
#ifndef UTILIB_HAVE_PTRDIFF_T 
#define UTILIB_HAVE_PTRDIFF_T  1 
#endif

/* Define to 1 if the system has the type `size_t'. */
#ifndef UTILIB_HAVE_SIZE_T 
#define UTILIB_HAVE_SIZE_T  1 
#endif

/* define if the compiler has stringstream */
#ifndef UTILIB_HAVE_SSTREAM 
#define UTILIB_HAVE_SSTREAM   
#endif

/* define if the compiler supports ISO C++ standard library */
#ifndef UTILIB_HAVE_STD 
#define UTILIB_HAVE_STD  1
#endif

/* Define to 1 if stdbool.h conforms to C99. */
#ifndef UTILIB_HAVE_STDBOOL_H 
#define UTILIB_HAVE_STDBOOL_H  1 
#endif

/* Define to 1 if you have the <stddef.h> header file. */
#ifndef UTILIB_HAVE_STDDEF_H 
#define UTILIB_HAVE_STDDEF_H  1 
#endif

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef UTILIB_HAVE_STDINT_H 
#define UTILIB_HAVE_STDINT_H  1 
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef UTILIB_HAVE_STDLIB_H 
#define UTILIB_HAVE_STDLIB_H  1 
#endif

/* Define to 1 if you have the `strcasecmp' function. */
#ifndef UTILIB_HAVE_STRCASECMP 
#define UTILIB_HAVE_STRCASECMP  1 
#endif

/* Define to 1 if you have the `strchr' function. */
#ifndef UTILIB_HAVE_STRCHR 
#define UTILIB_HAVE_STRCHR  1 
#endif

/* Define to 1 if you have the `stricmp' function. */
/* #undef UTILIB_HAVE_STRICMP */

/* Define to 1 if you have the <strings.h> header file. */
/* #ifndef UTILIB_HAVE_STRINGS_H */
/* #define UTILIB_HAVE_STRINGS_H  1 */ 
/* #endif */

/* Define to 1 if you have the <string.h> header file. */
#ifndef UTILIB_HAVE_STRING_H 
#define UTILIB_HAVE_STRING_H  1 
#endif

/* Define to 1 if you have the `sysconf' function. */
#ifndef UTILIB_HAVE_SYSCONF 
#define UTILIB_HAVE_SYSCONF  1 
#endif

/* Define to 1 if you have the <sys/resource.h> header file. */
/* #ifndef UTILIB_HAVE_SYS_RESOURCE_H */
/* #define UTILIB_HAVE_SYS_RESOURCE_H  1 */
/* #endif */

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef UTILIB_HAVE_SYS_STAT_H 
#define UTILIB_HAVE_SYS_STAT_H  1 
#endif

/* Define to 1 if you have the <sys/timeb.h> header file. */
#ifndef UTILIB_HAVE_SYS_TIMEB_H
#define UTILIB_HAVE_SYS_TIMEB_H  1
#endif

/* Define to 1 if you have the <sys/time.h> header file. */
/* #ifndef UTILIB_HAVE_SYS_TIME_H */
/* #define UTILIB_HAVE_SYS_TIME_H  1 */
/* #endif */

/* Define to 1 if you have the <sys/types.h> header file. */
/* #ifndef UTILIB_HAVE_SYS_TYPES_H */
/* #define UTILIB_HAVE_SYS_TYPES_H  1 */
/* #endif */

/* Define to 1 if you have the `times' function. */
/* #undef UTILIB_HAVE_TIMES */

/* Define to 1 if you have the <unistd.h> header file. */
/* #ifndef UTILIB_HAVE_UNISTD_H */
/* #define UTILIB_HAVE_UNISTD_H  1 */
/* #endif */

/* Define to 1 if you have the <values.h> header file. */
/* #undef UTILIB_HAVE_VALUES_H */

/* Define to 1 if the system has the type `_Bool'. */
#ifndef UTILIB_HAVE__BOOL 
#define UTILIB_HAVE__BOOL  1 
#endif

// /* software host will be cygwin */
// #ifndef UTILIB_HOST_CYGWIN 
// #define UTILIB_HOST_CYGWIN  1 
// #endif

/* software host is GNU */
/* #undef UTILIB_HOST_GNU */

/* software host will be linux */
/* #undef UTILIB_HOST_LINUX */

/* software host will be mingw */
/* #undef UTILIB_HOST_MINGW */

/* software host will be solaris */
/* #undef UTILIB_HOST_SOLARIS */

/* Name of package */
#ifndef UTILIB_PACKAGE 
#define UTILIB_PACKAGE  "utilib" 
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef UTILIB_PACKAGE_BUGREPORT 
#define UTILIB_PACKAGE_BUGREPORT  "acro-help@sandia.gov" 
#endif

/* Define to the full name of this package. */
#ifndef UTILIB_PACKAGE_NAME 
#define UTILIB_PACKAGE_NAME  "utilib" 
#endif

// Define to the full name and version of this package. 
#ifndef UTILIB_PACKAGE_STRING 
#define UTILIB_PACKAGE_STRING  "utilib VOTD" 
#endif

/* Define to the one symbol short name of this package. */
#ifndef UTILIB_PACKAGE_TARNAME 
#define UTILIB_PACKAGE_TARNAME  "utilib" 
#endif

/* Define to the version of this package. */
#ifndef UTILIB_PACKAGE_VERSION 
#define UTILIB_PACKAGE_VERSION  "VOTD" 
#endif

/* Define to 1 if you have the ANSI C header files. */
/* #undef UTILIB_STDC_HEADERS */

/* software target will be cygwin */
#ifndef UTILIB_TARGET_CYGWIN 
#define UTILIB_TARGET_CYGWIN  1 
#endif

/* software target will be linux */
/* #undef UTILIB_TARGET_LINUX */

/* software target will be mingw */
/* #undef UTILIB_TARGET_MINGW */

/* software target will be solaris */
/* #undef UTILIB_TARGET_SOLARIS */

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
/* #ifndef UTILIB_TIME_WITH_SYS_TIME */
/* #define UTILIB_TIME_WITH_SYS_TIME  1 */
/* #endif */

/* turn on code validation tests */
/* #undef UTILIB_VALIDATING */

/* Version number of package */
#ifndef UTILIB_VERSION 
#define UTILIB_VERSION  "VOTD" 
#endif

/* define whether checksum function is included in utilib */
#ifndef UTILIB_YES_CHECKSUM 
#define UTILIB_YES_CHECKSUM   
#endif

/* define whether CommonIO is included in utilib */
#ifndef UTILIB_YES_COMMONIO 
#define UTILIB_YES_COMMONIO   
#endif

/* define whether linpack will be included */
/* #undef UTILIB_YES_LINPACK */

/* define whether memdebug is included in utilib */
/* #undef UTILIB_YES_MEMDEBUG */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef _utilib_const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef _utilib_inline */
#endif

/* once: _UTILIB_UTILIB_CONFIG_H */
#endif