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

#else
    #define DISABLE_WARNING_PUSH
    #define DISABLE_WARNING_POP
    #define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER
    #define DISABLE_WARNING_UNREFERENCED_FUNCTION
    // other warnings you want to deactivate...

#endif
// clang-format on