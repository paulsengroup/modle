// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// Source: https://www.fluentcpp.com/2019/08/30/how-to-disable-a-warning-in-cpp/
// GCC to MSVC codes: https://github.com/srz-zumix/awesome-cpp-warning

// clang-format off

// Defines for GCC and Clang
#if defined(__GNUC__) || defined(__clang__)
    #define DO_PRAGMA(X)                     _Pragma(#X)                                    // NOLINT(cppcoreguidelines-macro-usage)
    #define DISABLE_WARNING_PUSH             DO_PRAGMA(GCC diagnostic push)
    #define DISABLE_WARNING_POP              DO_PRAGMA(GCC diagnostic pop)
    #define DISABLE_WARNING(warningName)     DO_PRAGMA(GCC diagnostic ignored warningName)  // NOLINT(cppcoreguidelines-macro-usage)

    #define DISABLE_WARNING_SIGN_CONVERSION  DISABLE_WARNING("-Wsign-conversion")
    #define DISABLE_WARNING_SIGN_COMPARE     DISABLE_WARNING("-Wsign-compare")
    #define DISABLE_WARNING_CONVERSION       DISABLE_WARNING("-Wconversion")
    #define DISABLE_WARNING_DOUBLE_PROMOTION DISABLE_WARNING("-Wdouble-promotion")
    #define DISABLE_WARNING_PADDED           DISABLE_WARNING("-Wpadded")
    #define DISABLE_WARNING_FLOAT_EQUAL      DISABLE_WARNING("-Wfloat-equal")
    #define DISABLE_WARNING_SWITCH_ENUM      DISABLE_WARNING("-Wswitch-enum")
    #define DISABLE_WARNING_SHADOW           DISABLE_WARNING("-Wshadow")
    #define DISABLE_WARNING_ATTRIBUTES       DISABLE_WARNING("-Wattributes")
    #define DISABLE_WARNING_UNUSED_VARIABLE  DISABLE_WARNING("-Wunused-variable")
    #define DISABLE_WARNING_NULL_DEREF       DISABLE_WARNING("-Wnull-dereference")
#endif

// Defines specific to Clang
#ifdef __clang__
    #define DISABLE_WARNING_BOOL_COMPARE           DISABLE_WARNING("-Wtautological-constant-out-of-range-compare")
    #define DISABLE_WARNING_SHORTEN_64_TO_32       DISABLE_WARNING("-Wshorten-64-to-32")
    #define DISABLE_WARNING_USED_BUT_MARKED_UNUSED DISABLE_WARNING("-Wused-but-marked-unused")
    #define DISABLE_WARNING_RANGE_LOOP_ANALYSIS    DISABLE_WARNING("-Wrange-loop-analysis")
    #define DISABLE_WARNING_DUPLICATED_BRANCHES
    #define DISABLE_WARNING_MAYBE_UNINITIALIZED

    #if defined(__APPLE__) || (__clang_major__ == 11)
        #define DISABLE_WARNING_IMPL_INT_TO_FLOAT  DISABLE_WARNING("-Wimplicit-const-int-float-conversion")
    #else
        #define DISABLE_WARNING_IMPL_INT_TO_FLOAT
    #endif

    #if __has_warning("-Wunused-but-set-parameter")
        #define DISABLE_WARNING_UNUSED_PARAMETER   DISABLE_WARNING("-Wunused-but-set-parameter")
    #else
        #define DISABLE_WARNING_UNUSED_PARAMETER
    #endif
#endif

// Defines specific to GCC
#if defined(__GNUC__) && !defined(__clang__)
    #define DISABLE_WARNING_BOOL_COMPARE           DISABLE_WARNING("-Wbool-compare")
    #define DISABLE_WARNING_IMPL_INT_TO_FLOAT      DISABLE_WARNING("-Wfloat-conversion")
    #define DISABLE_WARNING_MAYBE_UNINITIALIZED    DISABLE_WARNING("-Wmaybe-uninitialized")
    #define DISABLE_WARNING_DUPLICATED_BRANCHES    DISABLE_WARNING("-Wduplicated-branches")
    #define DISABLE_WARNING_UNUSED_PARAMETER       DISABLE_WARNING("-Wunused-but-set-parameter")
    #define DISABLE_WARNING_SHORTEN_64_TO_32
    #define DISABLE_WARNING_USED_BUT_MARKED_UNUSED
    #define DISABLE_WARNING_RANGE_LOOP_ANALYSIS
#endif

// Defines for unknown/unsupported compilers
#if !defined(_MSC_VER) && !defined(__GNUC__) && !defined(__clang__)
    #define DISABLE_WARNING
    #define DISABLE_WARNING_PUSH
    #define DISABLE_WARNING_POP

    #define DISABLE_WARNING_BOOL_COMPARE
    #define DISABLE_WARNING_CONVERSION
    #define DISABLE_WARNING_DOUBLE_PROMOTION
    #define DISABLE_WARNING_IMPL_INT_TO_FLOAT
    #define DISABLE_WARNING_SHORTEN_64_TO_32
    #define DISABLE_WARNING_SIGN_COMPARE
    #define DISABLE_WARNING_SIGN_CONVERSION
    #define DISABLE_WARNING_PADDED
    #define DISABLE_WARNING_USED_BUT_MARKED_UNUSED
    #define DISABLE_WARNING_FLOAT_EQUAL
    #define DISABLE_WARNING_SWITCH_ENUM
    #define DISABLE_WARNING_SHADOW
    #define DISABLE_WARNING_MAYBE_UNINITIALIZED
    #define DISABLE_WARNING_ATTRIBUTES
    #define DISABLE_WARNING_DUPLICATED_BRANCHES
    #define DISABLE_WARNING_UNUSED_VARIABLE
    #define DISABLE_WARNING_UNUSED_PARAMETER
    #define DISABLE_WARNING_RANGE_LOOP_ANALYSIS
    #define DISABLE_WARNING_NULL_DEREF
#endif

// clang-format on
