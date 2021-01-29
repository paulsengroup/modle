#pragma once

// Source: https://www.fluentcpp.com/2019/08/30/how-to-disable-a-warning-in-cpp/
// clang-format off

#if defined(_MSC_VER)
    #define DISABLE_WARNING_PUSH           __pragma(warning( push ))
    #define DISABLE_WARNING_POP            __pragma(warning( pop ))
    #define DISABLE_WARNING(warningNumber) __pragma(warning( disable : warningNumber ))

    #define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER    DISABLE_WARNING(4100)
    #define DISABLE_WARNING_UNREFERENCED_FUNCTION            DISABLE_WARNING(4505)
    // other warnings you want to deactivate...

#elif defined(__GNUC__) || defined(__clang__)
  #define DO_PRAGMA(X) _Pragma(#X)
  #define DISABLE_WARNING_PUSH           DO_PRAGMA(GCC diagnostic push)
  #define DISABLE_WARNING_POP            DO_PRAGMA(GCC diagnostic pop)
  #define DISABLE_WARNING(warningName)   DO_PRAGMA(GCC diagnostic ignored warningName)

  #define DISABLE_WARNING_SIGN_CONVERSION    DISABLE_WARNING("-Wsign-conversion")
  #define DISABLE_WARNING_SIGN_COMPARE       DISABLE_WARNING("-Wsign-compare")
  #define DISABLE_WARNING_CONVERSION         DISABLE_WARNING("-Wconversion")
  #define DISABLE_WARNING_ARITH_CONVERSION   DISABLE_WARNING("-Warith-conversion")
  #define DISABLE_WARNING_DOUBLE_PROMOTION   DISABLE_WARNING("-Wdouble-promotion")
  #define DISABLE_WARNING_UNUSED_FUNCTION    DISABLE_WARNING("-Wunused-function")
  #define DISABLE_WARNING_UNUSED_PARAMETER   DISABLE_WARNING("-Wunused-parameter")
  #define DISABLE_WARNING_SHADOW             DISABLE_WARNING("-Wshadow")
#if defined(__clang__)
    #define DISABLE_WARNING_C_VLA              DISABLE_WARNING("-Wvla-extension")
    #define DISABLE_WARNING_BOOL_COMPARE       (void)0;
    #define DISABLE_WARNING_IMPL_INT_TO_FLOAT  DISABLE_WARNING("-Wimplicit-const-int-float-conversion")
    #define DISABLE_WARNING_USELESS_CAST       (void)0;
#else
    #define DISABLE_WARNING_C_VLA             DISABLE_WARNING("-Wvla")
    #define DISABLE_WARNING_BOOL_COMPARE      DISABLE_WARNING("-Wbool-compare")
    #define DISABLE_WARNING_IMPL_INT_TO_FLOAT DISABLE_WARNING("-Wfloat-conversion")
    #define DISABLE_WARNING_USELESS_CAST       DISABLE_WARNING("-Wuseless-cast")
#endif
#else
    #define DISABLE_WARNING_PUSH
    #define DISABLE_WARNING_POP
    #define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER
    #define DISABLE_WARNING_UNREFERENCED_FUNCTION
    // other warnings you want to deactivate...

#endif
// clang-format on
