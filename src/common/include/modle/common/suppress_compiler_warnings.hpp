#pragma once

// Source: https://www.fluentcpp.com/2019/08/30/how-to-disable-a-warning-in-cpp/
// clang-format off

#if defined(_MSC_VER)
    #define DISABLE_WARNING_PUSH           __pragma(warning( push ))
    #define DISABLE_WARNING_POP            __pragma(warning( pop ))
    #define DISABLE_WARNING(warningNumber) __pragma(warning( disable : warningNumber ))

    // #define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER    DISABLE_WARNING(4100)
    // #define DISABLE_WARNING_UNREFERENCED_FUNCTION            DISABLE_WARNING(4505)

#elif defined(__GNUC__) || defined(__clang__)
  #define DO_PRAGMA(X) _Pragma(#X)
  #define DISABLE_WARNING_PUSH                   DO_PRAGMA(GCC diagnostic push)
  #define DISABLE_WARNING_POP                    DO_PRAGMA(GCC diagnostic pop)
  #define DISABLE_WARNING(warningName)           DO_PRAGMA(GCC diagnostic ignored warningName)

  #define DISABLE_WARNING_SIGN_CONVERSION        DISABLE_WARNING("-Wsign-conversion")
  #define DISABLE_WARNING_SIGN_COMPARE           DISABLE_WARNING("-Wsign-compare")
  #define DISABLE_WARNING_CONVERSION             DISABLE_WARNING("-Wconversion")
  #define DISABLE_WARNING_ARITH_CONVERSION       DISABLE_WARNING("-Warith-conversion")
  #define DISABLE_WARNING_DOUBLE_PROMOTION       DISABLE_WARNING("-Wdouble-promotion")
  #define DISABLE_WARNING_UNUSED_FUNCTION        DISABLE_WARNING("-Wunused-function")
  #define DISABLE_WARNING_UNUSED_PARAMETER       DISABLE_WARNING("-Wunused-parameter")
  #define DISABLE_WARNING_SHADOW                 DISABLE_WARNING("-Wshadow")

#if defined(__clang__)
    #define DISABLE_WARNING_C_VLA                DISABLE_WARNING("-Wvla-extension")
    #define DISABLE_WARNING_BOOL_COMPARE
    #define DISABLE_WARNING_USELESS_CAST
    #define DISABLE_WARNING_DUPLICATED_BANCHES

#if defined(__APPLE__) || (__clang_major__ < 11)
      #define DISABLE_WARNING_IMPL_INT_TO_FLOAT
#else
      #define DISABLE_WARNING_IMPL_INT_TO_FLOAT  DISABLE_WARNING("-Wimplicit-const-int-float-conversion")
#endif

#else
    #define DISABLE_WARNING_C_VLA                DISABLE_WARNING("-Wvla")
    #define DISABLE_WARNING_BOOL_COMPARE         DISABLE_WARNING("-Wbool-compare")
    #define DISABLE_WARNING_IMPL_INT_TO_FLOAT    DISABLE_WARNING("-Wfloat-conversion")
    #define DISABLE_WARNING_USELESS_CAST         DISABLE_WARNING("-Wuseless-cast")
    #define DISABLE_WARNING_DUPLICATED_BANCHES   DISABLE_WARNING("-Wduplicated-branches")
#endif
#else
    // Defines for unknown/unsupported compilers
    #define DISABLE_WARNING_PUSH
    #define DISABLE_WARNING_POP
    #define DISABLE_WARNING_SIGN_CONVERSION
    #define DISABLE_WARNING_SIGN_COMPARE
    #define DISABLE_WARNING_CONVERSION
    #define DISABLE_WARNING_ARITH_CONVERSION
    #define DISABLE_WARNING_DOUBLE_PROMOTION
    #define DISABLE_WARNING_UNUSED_FUNCTION
    #define DISABLE_WARNING_UNUSED_PARAMETER
    #define DISABLE_WARNING_SHADOW
    #define DISABLE_WARNING_C_VLA
    #define DISABLE_WARNING_BOOL_COMPARE
    #define DISABLE_WARNING_USELESS_CAST
    #define DISABLE_WARNING_IMPL_INT_TO_FLOAT

#endif
// clang-format on