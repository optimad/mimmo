/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
 \*---------------------------------------------------------------------------*/

#ifndef __MIMMO_COMPILER_HPP__
#define __MIMMO_COMPILER_HPP__

/**
 * Unreachable macro. Useful for suppressing "control reaches end of non-void
 * function" warnings.
 */
#if __GNUC__ >= 4 && __GNUC_MINOR__ >= 5
#define MIMMO_UNREACHABLE(str)    \
do {                         \
    assert(!str);            \
    __builtin_unreachable(); \
} while (0)
#elif (defined(__clang__) && defined(__has_builtin))
# if __has_builtin(__builtin_unreachable)
#  define MIMMO_UNREACHABLE(str)  \
do {                         \
    assert(!str);            \
    __builtin_unreachable(); \
} while (0)
# endif
#endif

#ifndef MIMMO_UNREACHABLE
#define MIMMO_UNREACHABLE(str)
#endif

/*!
 * Unused macro.
 */
#define MIMMO_UNUSED(variable)     \
do {                  \
    (void)(variable); \
} while (0)


/*!
 * Deprecated macro
 */
#if defined __GNUC__
#   define MIMMO_DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined __clang__
#   define MIMMO_DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
    #define MIMMO_DEPRECATED(func) __declspec(deprecated) func
#else
#   pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#   define MIMMO_DEPRECATED(func) func
#endif

#endif
