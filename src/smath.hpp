/**
 * Description:  Single-header math library
 * Author:       Alicia Amarilla (smushy) 
 * File Created: November 23, 2022 
*/
#pragma once
#include <cstdint>
#include <immintrin.h>

// TODO(alicia): SIMD!
#ifndef SMUSHY_TYPE_ALIAS
#define SMUSHY_TYPE_ALIAS 1

// TYPEDEF -----------------------------------------------------------------------------------------
// Pointer-sized unsigned integer
typedef uintptr_t usize;
// Pointer-sized signed integer
typedef intptr_t  isize;

// 8-bit unsigned integer
typedef uint8_t   u8;
// 16-bit unsigned integer
typedef uint16_t u16;
// 32-bit unsigned integer
typedef uint32_t u32;
// 64-bit unsigned integer
typedef uint64_t u64;
// 32-bit unsigned integer for representing boolean values
typedef uint32_t bool32;

// 8-bit signed integer
typedef int8_t   i8;
// 16-bit signed integer
typedef int16_t i16;
// 32-bit signed integer
typedef int32_t i32;
// 64-bit signed integer
typedef int64_t i64;

// 32-bit floating-point number
typedef float  f32;
// 64-bit floating-point number
typedef double f64;

// Float 32-bit constants
namespace F32 {
    // Largest finite f32 value
    inline const constexpr f32 MAX = 3.40282347E+38f;
    // Smallest finite f32 value
    inline const constexpr f32 MIN = -3.40282347E+38f;
    // Not a number
    inline const f32 NaN = 0.0f / 0.0f;
    // Smallest positive f32 value
    inline const constexpr f32 MIN_POS = 1.17549435E-38f;
    // Positive infinity
    inline const f32 POS_INFINITY = 1.0f / 0.0f;
    // Positive infinity
    inline const f32 NEG_INFINITY = -1.0f / 0.0f;
    // Pi constant
    inline const constexpr f32 PI = 3.141592741f;
    // Tau constant
    inline const constexpr f32 TAU = 2.0f * PI;
};

// Float 64-bit constants
namespace F64 {
    // Largest finite f64 value
    inline const constexpr f64 MAX = 1.7976931348623157E+308;
    // Smallest finite f64 value
    inline const constexpr f64 MIN = -1.7976931348623157E+308;
    // Not a number
    inline const f64 NaN = 0.0 / 0.0;
    // Smallest positive f32 value
    inline const constexpr f64 MIN_POS = 2.2250738585072014E-308;
    // Positive infinity
    inline const f64 POS_INFINITY = 1.0 / 0.0;
    // Positive infinity
    inline const f64 NEG_INFINITY = -1.0 / 0.0;
    // Pi constant
    inline const constexpr f64 PI = 3.14159265358979323846;
    // Tau constant
    inline const constexpr f64 TAU = 2.0 * PI;
};

// Unsigned integer 8-bit constants
namespace U8 {
    // Largest u8 value
    inline const constexpr u8 MAX = 255;
    // Smallest u8 value
    inline const constexpr u8 MIN = 0;
};

// Unsigned integer 16-bit constants
namespace U16 {
    // Largest u16 value
    inline const constexpr u16 MAX = 65535;
    // Smallest u16 value
    inline const constexpr u16 MIN = 0;
};

// Unsigned integer 32-bit constants
namespace U32 {
    // Largest u32 value
    inline const constexpr u32 MAX = 4294967295;
    // Smallest u32 value
    inline const constexpr u32 MIN = 0;
};

// Unsigned integer 64-bit constants
namespace U64 {
    // Largest u64 value
    inline const constexpr u64 MAX = 18446744073709551615ULL;
    // Smallest u64 value
    inline const constexpr u64 MIN = 0ULL;
};

// Integer 8-bit constants
namespace I8 {
    // Largest i8 value
    inline const constexpr i8 MAX = 127;
    // Smallest i8 value
    inline const constexpr i8 MIN = -128;
};

// Integer 16-bit constants
namespace I16 {
    // Largest i16 value
    inline const constexpr i16 MAX = 32767;
    // Smallest i16 value
    inline const constexpr i16 MIN = -32768;
};

// Integer 32-bit constants
namespace I32 {
    // Largest i32 value
    inline const constexpr i32 MAX = 2147483647;
    // Smallest i32 value
    inline const constexpr i32 MIN = -2147483648;
};

// Integer 64-bit constants
namespace I64 {
    // Largest i64 value
    inline const constexpr i64 MAX = 9223372036854775807LL;
    // Smallest i64 value
    inline const constexpr i64 MIN = -9223372036854775807 - 1;
};

#endif

namespace smath {

// NOTE(alicia): FUNCTIONS ---------------------------------------------------------------------------------------

// forward-declarations

inline f32 sqrt( f32 x );
inline f64 sqrt( f64 x );
inline constexpr f32 abs( f32 x );
inline constexpr f64 abs( f64 x );
inline constexpr f32 sign( f32 x );
inline constexpr f64 sign( f64 x );

// sine of x
inline f32 sin( f32 x ) {
    // TODO(alicia): temp?
    return __builtin_sinf(x);
}
// sine of x
inline f64 sin( f64 x ) {
    // TODO(alicia): temp?
    return __builtin_sin(x);
}
// cosine of x
inline f32 cos( f32 x ) {
    // TODO(alicia): temp?
    return __builtin_cosf(x);
}
// cosine of x
inline f64 cos( f64 x ) {
    // TODO(alicia): temp?
    return __builtin_cos(x);
}
// tangent of x
inline f32 tan( f32 x ) {
    // TODO(alicia): temp?
    return __builtin_tanf(x);
}
// tangent of x
inline f64 tan( f64 x ) {
    // TODO(alicia): temp?
    return __builtin_tan(x);
}
/// arc-sine 
inline f32 asin( f32 x ) {
    // TODO(alicia): temp?
    return __builtin_asinf(x);
}
/// arc-sine 
inline f64 asin( f64 x ) {
    // TODO(alicia): temp?
    return __builtin_asin(x);
}
/// arc-sine, returns pi/2 * sign(x) instead of NaN when outside -1 to 1 range
inline f32 asinNoNaN( f32 x ) {
    if( abs( x ) >= 1.0f ) {
        return ( F32::PI / 2.0f ) * sign(x);
    } else {
        return asin(x);
    }
}
/// arc-sine, returns pi/2 * sign(x) instead of NaN when outside -1 to 1 range
inline f64 asinNoNaN( f64 x ) {
    if( abs( x ) >= 1.0 ) {
        return ( F64::PI / 2.0 ) * sign(x);
    } else {
        return asin(x);
    }
}
/// arc-cosine
inline f32 acos( f32 x ) {
    // TODO(alicia): temp?
    return __builtin_acosf( x );
}
/// arc-cosine
inline f64 acos( f64 x ) {
    // TODO(alicia): temp?
    return __builtin_acos( x );
}
// 2 argument arc-tangent
inline f32 atan2( f32 a, f32 b ) {
    // TODO(alicia): temp?
    return __builtin_atan2f( a, b );
}
// 2 argument arc-tangent
inline f64 atan2( f64 a, f64 b ) {
    // TODO(alicia): temp?
    return __builtin_atan2( a, b );
}

// degrees to radians
inline constexpr f32 toRad( f32 deg ) {
    return deg * ( F32::PI / 180.0f );
}
// radians to degrees
inline constexpr f32 toDeg( f32 rad ) {
    return rad * ( 180.0f / F32::PI );
}

inline i64 rand( i64 seed ) {
    // linear congruential generator :)
    // TODO(alicia): do more research on rng's

    i64 a = 122;
    i64 m = 6454;
    i64 b = 455;

    static i64 inc = 0;
    i64 result = (a * seed + inc) % m;
    inc = (inc + b) % m;

    return result;
}

/// square root
inline f32 sqrt( f32 x ) {
    // NOTE(alicia): SSE

    __m128 temp = _mm_set_ss(x);
    temp = _mm_sqrt_ss( temp );
    return _mm_cvtss_f32( temp );
}
// square root
inline f64 sqrt( f64 x ) {
    // NOTE(alicia): SSE

    __m128d temp = _mm_set_sd( x );
    temp = _mm_sqrt_pd( temp );
    return _mm_cvtsd_f64( temp );
}

// raise to the power of 2
inline constexpr i8 sqr( i8 x ) {
    return x * x;
}
// raise to the power of 2
inline constexpr i16 sqr( i16 x ) {
    return x * x;
}
// raise to the power of 2
inline constexpr i32 sqr( i32 x ) {
    return x * x;
}
// raise to the power of 2
inline constexpr i64 sqr( i64 x ) {
    return x * x;
}
// raise to the power of 2
inline constexpr u8 sqr( u8 x ) {
    return x * x;
}
// raise to the power of 2
inline constexpr u16 sqr( u16 x ) {
    return x * x;
}
// raise to the power of 2
inline constexpr u32 sqr( u32 x ) {
    return x * x;
}
// raise to the power of 2
inline constexpr u64 sqr( u64 x ) {
    return x * x;
}
// raise to the power of 2
inline constexpr f32 sqr( f32 x ) {
    return x * x;
}
// raise to the power of 2
inline constexpr f64 sqr( f64 x ) {
    return x * x;
}

/// @brief Clamp an int between min and max values
/// @param value value to clamp
/// @param min minimum range, inclusive
/// @param max maximum range, inclusive
/// @return clamped value
inline constexpr i8 clamp( i8 value, i8 min, i8 max ) {
    const i8 t = value < min ? min : value;
    return t > max ? max : t;
}
/// @brief Clamp an int between min and max values
/// @param value value to clamp
/// @param min minimum range, inclusive
/// @param max maximum range, inclusive
/// @return clamped value
inline constexpr i16 clamp( i16 value, i16 min, i16 max ) {
    const i16 t = value < min ? min : value;
    return t > max ? max : t;
}
/// @brief Clamp an int between min and max values
/// @param value value to clamp
/// @param min minimum range, inclusive
/// @param max maximum range, inclusive
/// @return clamped value
inline constexpr i32 clamp( i32 value, i32 min, i32 max ) {
    const i32 t = value < min ? min : value;
    return t > max ? max : t;
}
/// @brief Clamp an int between min and max values
/// @param value value to clamp
/// @param min minimum range, inclusive
/// @param max maximum range, inclusive
/// @return clamped value
inline constexpr i64 clamp( i64 value, i64 min, i64 max ) {
    const i64 t = value < min ? min : value;
    return t > max ? max : t;
}
/// @brief Clamp an int between min and max values
/// @param value value to clamp
/// @param min minimum range, inclusive
/// @param max maximum range, inclusive
/// @return clamped value
inline constexpr u8 clamp( u8 value, u8 min, u8 max ) {
    const u8 t = value < min ? min : value;
    return t > max ? max : t;
}
/// @brief Clamp an int between min and max values
/// @param value value to clamp
/// @param min minimum range, inclusive
/// @param max maximum range, inclusive
/// @return clamped value
inline constexpr u16 clamp( u16 value, u16 min, u16 max ) {
    const u16 t = value < min ? min : value;
    return t > max ? max : t;
}
/// @brief Clamp an int between min and max values
/// @param value value to clamp
/// @param min minimum range, inclusive
/// @param max maximum range, inclusive
/// @return clamped value
inline constexpr u32 clamp( u32 value, u32 min, u32 max ) {
    const u32 t = value < min ? min : value;
    return t > max ? max : t;
}
/// @brief Clamp an int between min and max values
/// @param value value to clamp
/// @param min minimum range, inclusive
/// @param max maximum range, inclusive
/// @return clamped value
inline constexpr u64 clamp( u64 value, u64 min, u64 max ) {
    const u64 t = value < min ? min : value;
    return t > max ? max : t;
}
/// @brief Clamp an int between min and max values
/// @param value value to clamp
/// @param min minimum range, inclusive
/// @param max maximum range, inclusive
/// @return clamped value
inline constexpr f32 clamp( f32 value, f32 min, f32 max ) {
    const f32 t = value < min ? min : value;
    return t > max ? max : t;
}
/// @brief Clamp an int between min and max values
/// @param value value to clamp
/// @param min minimum range, inclusive
/// @param max maximum range, inclusive
/// @return clamped value
inline constexpr f64 clamp( f64 value, f64 min, f64 max ) {
    const f64 t = value < min ? min : value;
    return t > max ? max : t;
}

/// @brief Linear interpolation
/// @param a minimum value
/// @param b maximum value
/// @param t fraction
/// @return blend between a and b, based on fraction t
inline f32 lerp( f32 a, f32 b, f32 t ) {
    return ( 1.0f - t ) * a + b * t;
}
/// @brief Linear interpolation
/// @param a minimum value
/// @param b maximum value
/// @param t fraction
/// @return blend between a and b, based on fraction t
inline f64 lerp( f64 a, f64 b, f64 t ) {
    return ( 1.0 - t ) * a + b * t;
}
/// @brief Linear interpolation, t is clamped between 0.0-1.0
/// @param a minimum value
/// @param b maximum value
/// @param t 0.0-1.0 fraction
/// @return blend between a and b, based on fraction t
inline f32 clampedLerp( f32 a, f32 b, f32 t ) {
    return lerp( a, b, clamp( t, 0.0f, 1.0f ) );
}
/// @brief Linear interpolation, t is clamped between 0.0-1.0
/// @param a minimum value
/// @param b maximum value
/// @param t 0.0-1.0 fraction
/// @return blend between a and b, based on fraction t
inline f64 clampedLerp( f64 a, f64 b, f64 t ) {
    return lerp( a, b, clamp( t, 0.0, 1.0 ) );
}
/// @brief Inverse linear interpolation
/// @param a minimum
/// @param b maximum
/// @param v value between
/// @return fraction that value occupies between a and b
inline f32 invLerp( f32 a, f32 b, f32 v ) {
    return ( v - a ) / ( b - a );
}
/// @brief Inverse linear interpolation
/// @param a minimum
/// @param b maximum
/// @param v value between
/// @return fraction that value occupies between a and b
inline f64 invLerp( f64 a, f64 b, f64 v ) {
    return ( v - a ) / ( b - a );
}
/// @brief Remap value from input range to output range
/// @param imin input range minimum
/// @param imax input range maximum
/// @param omin output range minimum
/// @param omax output range maximum
/// @param v value to remap
/// @return remapped value
inline f32 remap( f32 imin, f32 imax, f32 omin, f32 omax, f32 v ) {
    const f32 t = invLerp( imin, imax, v );
    return lerp( omin, omax, t );
}
/// @brief Remap value from input range to output range
/// @param imin input range minimum
/// @param imax input range maximum
/// @param omin output range minimum
/// @param omax output range maximum
/// @param v value to remap
/// @return remapped value
inline f64 remap( f64 imin, f64 imax, f64 omin, f64 omax, f64 v ) {
    const f64 t = invLerp( imin, imax, v );
    return lerp( omin, omax, t );
}

// normalize integer from -1.0-1.0
inline constexpr f32 normalize( i8 value ) {
    if( value > 0 ) {
        return (f32)value / (f32)I8::MAX;
    } else {
        return (f32)value / -(f32)I8::MIN;
    }
}
// normalize integer from -1.0-1.0
inline constexpr f32 normalize( i16 value ) {
    if( value > 0 ) {
        return (f32)value / (f32)I16::MAX;
    } else {
        return (f32)value / -(f32)I16::MIN;
    }
}
// normalize integer from -1.0-1.0
inline constexpr f32 normalize( i32 value ) {
    if( value > 0 ) {
        return (f32)value / (f32)I32::MAX;
    } else {
        return (f32)value / -(f32)I32::MIN;
    }
}
// normalize integer from -1.0-1.0
inline constexpr f32 normalize( i64 value ) {
    if( value > 0 ) {
        return (f32)value / (f32)I64::MAX;
    } else {
        return (f32)value / -(f32)I64::MIN;
    }
}
// normalize unsigned integer from 0.0-1.0
inline constexpr f32 normalize( u8 value ) {
    return (f32)value / (f32)U8::MAX;
}
// normalize unsigned integer from 0.0-1.0
inline constexpr f32 normalize( u16 value ) {
    return (f32)value / (f32)U16::MAX;
}
// normalize unsigned integer from 0.0-1.0
inline constexpr f32 normalize( u32 value ) {
    return (f32)value / (f32)U32::MAX;
}
// normalize unsigned integer from 0.0-1.0
inline constexpr f32 normalize( u64 value ) {
    return (f32)value / (f32)U64::MAX;
}

// sign of value
inline constexpr i8 sign( i8 value ) {
    return ( value > 0 ) - ( value < 0 );
}
// sign of value
inline constexpr i16 sign( i16 value ) {
    return ( value > 0 ) - ( value < 0 );
}
// sign of value
inline constexpr i32 sign( i32 value ) {
    return ( value > 0 ) - ( value < 0 );
}
// sign of value
inline constexpr i64 sign( i64 value ) {
    return ( value > 0 ) - ( value < 0 );
}
// sign of value
inline constexpr f32 sign( f32 value ) {
    return ( value > 0.0f ) - ( value < 0.0f );
}
// sign of value
inline constexpr f64 sign( f64 value ) {
    return ( value > 0.0 ) - ( value < 0.0 );
}

// absolute value
inline constexpr i8 abs( i8 value ) {
    return value * sign(value);
}
// absolute value
inline constexpr i16 abs( i16 value ) {
    return value * sign(value);
}
// absolute value
inline constexpr i32 abs( i32 value ) {
    return value * sign(value);
}
// absolute value
inline constexpr i64 abs( i64 value ) {
    return value * sign(value);
}
// absolute value
inline constexpr f32 abs( f32 value ) {
    return value * sign(value);
}
// absolute value
inline constexpr f64 abs( f64 value ) {
    return value * sign(value);
}
/// @brief multiply 4 floats by 4 floats
/// @param x1,y1,z1,w1 a
/// @param x2,y2,z2,w2 b
/// @param out f32[4] a * b
inline void mul(
    f32 x1, f32 y1, f32 z1, f32 w1,
    f32 x2, f32 y2, f32 z2, f32 w2,
    f32* out
) {
    // NOTE(alicia): SSE
    __m128 a = _mm_set_ps( w1, z1, y1, x1 );
    __m128 b = _mm_set_ps( w2, z2, y2, x2 );
    _mm_storeu_ps( out, _mm_mul_ps( a, b ) );
}
/// @brief multiply 4 floats by scalar
/// @param x,y,z,w a
/// @param scalar b
/// @param out f32[4] a * b
inline void mul(
    f32 x, f32 y, f32 z, f32 w,
    f32 scalar,
    f32* out
) {
    // NOTE(alicia): SSE
    __m128 a = _mm_set_ps( w, z, y, x );
    __m128 b = _mm_set1_ps( scalar );
    _mm_storeu_ps( out, _mm_mul_ps( a, b ) );
}

// NOTE(alicia): TYPES -------------------------------------------------------------------------------------------

const f32 VEC_CMP_THRESHOLD = 0.0001f;

struct vec2;
struct vec3;
struct vec4;
struct quat;
struct mat3;
struct mat4;
// TODO(alicia): MAT3, MAT4

vec2 operator+( const vec2& lhs, const vec2& rhs );
vec2 operator-( const vec2& lhs, const vec2& rhs );
vec2 operator*( const vec2& lhs, f32 rhs );
vec2 operator*( f32 lhs, const vec2& rhs );
vec2 operator/( const vec2& lhs, f32 rhs );
// normalize vector 
inline vec2 normalize( const vec2& v );
// clamp vector length up to max, inclusive
inline vec2 clamp( const vec2& v, f32 max );
// compare two vectors
inline bool cmp( const vec2& lhs, const vec2& rhs );
// dot product
inline f32 dot( const vec2& lhs, const vec2& rhs );
// scale vectors component-wise
inline vec2 mul( const vec2& lhs, const vec2& rhs );
// square magnitude
inline f32 sqrMag( const vec2& v );
// magnitude
inline f32 mag( const vec2& v );
// signed angle between two vectors
inline f32 angle( const vec2& lhs, const vec2& rhs );
// unsigned angle between two vectors
inline f32 unsignedAngle( const vec2& lhs, const vec2& rhs );
/// @brief Linear interpolation
/// @param a minimum value
/// @param b maximum value
/// @param t fraction
/// @return blend between a and b, based on fraction t
inline vec2 lerp( const vec2& a, const vec2& b, f32 t );
/// @brief Linear interpolation, t is clamped between 0.0-1.0
/// @param a minimum value
/// @param b maximum value
/// @param t 0.0-1.0 fraction
/// @return blend between a and b, based on fraction t
inline vec2 clampedLerp( const vec2& a, const vec2& b, f32 t );
// 2-component vector
struct vec2 {
    union {
        struct { f32 x, y; };
        struct { f32 u, v; };
    };

    vec2() {}
    vec2( f32 scalar ) : x(scalar), y(scalar) {}
    vec2( f32 x, f32 y ) : x(x), y(y) {}

    // get pointer to struct as f32 
    f32* ptr() { return (f32*)this; }

    vec2& operator-() { return *this *= -1.0f; }
    bool operator==( const vec2& rhs ) { return cmp( *this, rhs ); }
    bool operator!=( const vec2& rhs ) { return !(*this == rhs); }

    vec2& operator+=( const vec2& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( 0.0f, 0.0f, this->x, this->y );
        __m128 b = _mm_set_ps( 0.0f, 0.0f, rhs.x, rhs.y );
        f32 result[2];
        _mm_storeu_ps( result, _mm_add_ps( a, b ) );

        this->x = result[1];
        this->y = result[0];
        return *this;
    }
    vec2& operator-=( const vec2& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( 0.0f, 0.0f, this->x, this->y );
        __m128 b = _mm_set_ps( 0.0f, 0.0f, rhs.x, rhs.y );
        f32 result[2];
        _mm_storeu_ps( result, _mm_sub_ps( a, b ) );

        this->x = result[1];
        this->y = result[0];
        return *this;
    }
    vec2& operator*=( const f32& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( 0.0f, 0.0f, this->x, this->y );
        __m128 b = _mm_set_ps( 0.0f, 0.0f, rhs, rhs );
        f32 result[2];
        _mm_storeu_ps( result, _mm_mul_ps( a, b ) );

        this->x = result[1];
        this->y = result[0];
        return *this;
    }
    vec2& operator/=( const f32& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( 0.0f, 0.0f, this->x, this->y );
        __m128 b = _mm_set_ps( 0.0f, 0.0f, rhs, rhs );
        f32 result[2];
        _mm_storeu_ps( result, _mm_div_ps( a, b ) );

        this->x = result[1];
        this->y = result[0];
        return *this;
    }

    // ( 1.0f, 1.0f ) 
    static vec2 one()   { return {  1.0f,  1.0f }; }
    // ( 0.0f, 0.0f ) 
    static vec2 zero()  { return {  0.0f,  0.0f }; }
    // ( -1.0f, 0.0f ) 
    static vec2 left()  { return { -1.0f,  0.0f }; }
    // ( 1.0f, 0.0f ) 
    static vec2 right() { return {  1.0f,  0.0f }; }
    // ( 0.0f, 1.0f ) 
    static vec2 up()    { return {  0.0f,  1.0f }; }
    // ( 0.0f, -1.0f ) 
    static vec2 down()  { return {  0.0f, -1.0f }; }
};
inline vec2 operator+( const vec2& lhs, const vec2& rhs ) {
    return vec2(lhs) += rhs;
}
inline vec2 operator-( const vec2& lhs, const vec2& rhs ) {
    return vec2(lhs) -= rhs;
}
inline vec2 operator*( const vec2& lhs, f32 rhs ) {
    return vec2(lhs) *= rhs;
}
inline vec2 operator*( f32 lhs, const vec2& rhs ) {
    return vec2(rhs) *= lhs;
}
inline vec2 operator/( const vec2& lhs, f32 rhs ) {
    return vec2(lhs) /= rhs;
}
inline vec2 normalize( const vec2& v ) {
    f32 m = mag( v );
    if( m != 0.0f ) {
        return v / m;
    } else {
        return vec2::zero();
    }
}
inline bool cmp( const vec2& lhs, const vec2& rhs ) {
    return sqrMag(lhs - rhs) < VEC_CMP_THRESHOLD;
}
inline vec2 clamp( const vec2& v, f32 max ) {
    f32 magnitude = mag(v);
    if( magnitude > max ) {
        vec2 result = v / magnitude;
        result *= max;
        return result;
    } else {
        return v;
    }
}
inline f32 dot( const vec2& lhs, const vec2& rhs ) {
    // NOTE(alicia): SSE
    // TODO(alicia): further optimizations? FMA?

    __m128 a = _mm_set_ps( 0.0f, 0.0f, lhs.x, lhs.y );
    __m128 b = _mm_set_ps( 0.0f, 0.0f, rhs.x, rhs.y );
    f32 result[2];
    _mm_storeu_ps( result, _mm_mul_ps( a, b ) );

    return result[0] + result[1];
}
inline vec2 mul( const vec2& lhs, const vec2& rhs ) {
    // NOTE(alicia): SSE

    __m128 a = _mm_set_ps( 0.0f, 0.0f, lhs.x, lhs.y );
    __m128 b = _mm_set_ps( 0.0f, 0.0f, rhs.x, rhs.y );
    f32 result[2];
    _mm_storeu_ps( result, _mm_mul_ps( a, b ) );

    return { result[1], result[0] };
}
inline f32 sqrMag( const vec2& v ) {
    // NOTE(alicia): SSE
    // TODO(alicia): further optimizations? FMA?

    __m128 a = _mm_set_ps( 0.0f, 0.0f, v.x, v.y );
    f32 result[2];
    _mm_storeu_ps( result, _mm_mul_ps( a, a ) );

    return result[0] + result[1];
}
inline f32 mag( const vec2& v ) {
    return sqrt( sqrMag( v ) );
}
inline f32 angle( const vec2& lhs, const vec2& rhs ) {
    return acos( dot( lhs, rhs ) );
}
inline f32 unsignedAngle( const vec2& lhs, const vec2& rhs ) {
    return abs( angle( lhs, rhs ) );
}
inline vec2 lerp( const vec2& a, const vec2& b, f32 t ) {
    return ( 1.0f - t ) * a + b * t;
}
inline vec2 clampedLerp( const vec2& a, const vec2& b, f32 t ) {
    return lerp( a, b, clamp( t, 0.0f, 1.0f ) );
}

vec3 operator+( const vec3& lhs, const vec3& rhs );
vec3 operator-( const vec3& lhs, const vec3& rhs );
vec3 operator*( const vec3& lhs, f32 rhs );
vec3 operator*( f32 lhs, const vec3& rhs );
vec3 operator/( const vec3& lhs, f32 rhs );
// normalize vector 
inline vec3 normalize( const vec3& v );
// clamp vector length up to max, inclusive
inline vec3 clamp( const vec3& v, f32 max );
// compare two vectors
inline bool cmp( const vec3& lhs, const vec3& rhs );
// dot product
inline f32 dot( const vec3& lhs, const vec3& rhs );
// cross product
inline vec3 cross( const vec3& lhs, const vec3& rhs );
// scale vectors component-wise
inline vec3 mul( const vec3& lhs, const vec3& rhs );
// square magnitude
inline f32 sqrMag( const vec3& v );
// magnitude
inline f32 mag( const vec3& v );
// signed angle between two vectors
inline f32 angle( const vec3& lhs, const vec3& rhs );
// unsigned angle between two vectors
inline f32 unsignedAngle( const vec3& lhs, const vec3& rhs );
// reflect vector off of given normal
inline vec3 reflect( const vec3& direction, const vec3& normal );
/// @brief Linear interpolation
/// @param a minimum value
/// @param b maximum value
/// @param t fraction
/// @return blend between a and b, based on fraction t
inline vec3 lerp( const vec3& a, const vec3& b, f32 t );
/// @brief Linear interpolation, t is clamped between 0.0-1.0
/// @param a minimum value
/// @param b maximum value
/// @param t 0.0-1.0 fraction
/// @return blend between a and b, based on fraction t
inline vec3 clampedLerp( const vec3& a, const vec3& b, f32 t );
// 3-component vector
struct vec3 {
    union {
        struct { f32 x, y, z; };
        struct { f32 r, g, b; };
    };

    vec3() {}
    vec3( f32 scalar ) : x(scalar), y(scalar), z(scalar) {}
    vec3( const vec2& v ) : x(v.x), y(v.y), z(0.0f) {}
    vec3( f32 x, f32 y, f32 z ) : x(x), y(y), z(z) {}

    // get pointer to struct as f32 
    f32* ptr() { return (f32*)this; }

    vec3& operator-() { return *this *= -1.0f; }
    bool operator==( const vec3& rhs ) { return cmp( *this, rhs ); }
    bool operator!=( const vec3& rhs ) { return !(*this == rhs); }
    vec3& operator+=( const vec3& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( 0.0f, this->x, this->y, this->z );
        __m128 b = _mm_set_ps( 0.0f, rhs.x, rhs.y, rhs.z );
        f32 result[3];
        _mm_storeu_ps( result, _mm_add_ps( a, b ) );

        this->x = result[2];
        this->y = result[1];
        this->z = result[0];
        return *this;
    }
    vec3& operator-=( const vec3& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( 0.0f, this->x, this->y, this->z );
        __m128 b = _mm_set_ps( 0.0f, rhs.x, rhs.y, rhs.z );
        f32 result[3];
        _mm_storeu_ps( result, _mm_sub_ps( a, b ) );

        this->x = result[2];
        this->y = result[1];
        this->z = result[0];
        return *this;
    }
    vec3& operator*=( const f32& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( 0.0f, this->x, this->y, this->z );
        __m128 b = _mm_set_ps( 0.0f, rhs, rhs, rhs );
        f32 result[3];
        _mm_storeu_ps( result, _mm_mul_ps( a, b ) );

        this->x = result[2];
        this->y = result[1];
        this->z = result[0];
        return *this;
    }
    vec3& operator/=( const f32& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( 0.0f, this->x, this->y, this->z );
        __m128 b = _mm_set_ps( 0.0f, rhs, rhs, rhs );
        f32 result[3];
        _mm_storeu_ps( result, _mm_div_ps( a, b ) );

        this->x = result[2];
        this->y = result[1];
        this->z = result[0];
        return *this;
    }

    // ( 1.0f, 1.0f, 1.0f ) 
    static vec3 one()     { return {  1.0f,  1.0f,  1.0f }; }
    // ( 0.0f, 0.0f, 0.0f ) 
    static vec3 zero()    { return {  0.0f,  0.0f,  0.0f }; }
    // ( -1.0f, 0.0f, 0.0f ) 
    static vec3 left()    { return { -1.0f,  0.0f,  0.0f }; }
    // ( 1.0f, 0.0f, 0.0f ) 
    static vec3 right()   { return {  1.0f,  0.0f,  0.0f }; }
    // ( 0.0f, 1.0f, 0.0f )
    static vec3 up()      { return {  0.0f,  1.0f,  0.0f }; }
    // ( 0.0f, -1.0f, 0.0f )
    static vec3 down()    { return {  0.0f, -1.0f,  0.0f }; }
    // ( 0.0f, 0.0f, 1.0f )
    static vec3 forward() { return {  0.0f,  0.0f,  1.0f }; }
    // ( 0.0f, 0.0f, -1.0f )
    static vec3 back()    { return {  0.0f, 0.0f,  -1.0f }; }
};
inline vec3 operator+( const vec3& lhs, const vec3& rhs ) {
    return vec3(lhs) += rhs;
}
inline vec3 operator-( const vec3& lhs, const vec3& rhs ) {
    return vec3(lhs) -= rhs;
}
inline vec3 operator*( const vec3& lhs, f32 rhs ) {
    return vec3(lhs) *= rhs;
}
inline vec3 operator*( f32 lhs, const vec3& rhs ) {
    return vec3(rhs) *= lhs;
}
inline vec3 operator/( const vec3& lhs, f32 rhs ) {
    return vec3(lhs) /= rhs;
}
inline vec3 normalize( const vec3& v ) {
    f32 m = mag( v );
    if( m != 0.0f ) {
        return v / m;
    } else {
        return vec3::zero();
    }
}
inline bool cmp( const vec3& lhs, const vec3& rhs ) {
    return sqrMag(lhs - rhs) < VEC_CMP_THRESHOLD;
}
inline vec3 clamp( const vec3& v, f32 max ) {
    f32 magnitude = mag(v);
    if( magnitude > max ) {
        vec3 result = v / magnitude;
        result *= max;
        return result;
    } else {
        return v;
    }
}
inline f32 dot( const vec3& lhs, const vec3& rhs ) {
    // NOTE(alicia): SSE
    // TODO(alicia): further optimizations? FMA?

    __m128 a = _mm_set_ps( 0.0f, lhs.x, lhs.y, lhs.z );
    __m128 b = _mm_set_ps( 0.0f, rhs.x, rhs.y, rhs.z );
    f32 result[3];
    _mm_storeu_ps( result, _mm_mul_ps( a, b ) );

    return result[0] + result[1] + result[2];
}
inline vec3 cross( const vec3& lhs, const vec3& rhs ) {
    // TODO(alicia): SIMD!
    return {
        ( lhs.y * rhs.z ) - ( lhs.z * rhs.y ),
        ( lhs.z * rhs.x ) - ( lhs.x * rhs.z ),
        ( lhs.x * rhs.y ) - ( lhs.y * rhs.x )
    };
}
inline vec3 mul( const vec3& lhs, const vec3& rhs ) {
    // NOTE(alicia): SSE

    __m128 a = _mm_set_ps( 0.0f, lhs.x, lhs.y, lhs.z );
    __m128 b = _mm_set_ps( 0.0f, rhs.x, rhs.y, rhs.z );
    f32 result[3];
    _mm_storeu_ps( result, _mm_mul_ps( a, b ) );

    return {result[2], result[1], result[0]};
}
inline f32 sqrMag( const vec3& v ) {
    // NOTE(alicia): SSE
    // TODO(alicia): further optimizations? FMA?

    __m128 a = _mm_set_ps( 0.0f, v.x, v.y, v.z );
    f32 result[3];
    _mm_storeu_ps( result, _mm_mul_ps( a, a ) );

    return result[0] + result[1] + result[2];
}
inline f32 mag( const vec3& v ) {
    return sqrt( sqrMag( v ) );
}
inline f32 angle( const vec3& lhs, const vec3& rhs ) {
    return acos( dot( lhs, rhs ) );
}
inline f32 unsignedAngle( const vec3& lhs, const vec3& rhs ) {
    return abs( angle( lhs, rhs ) );
}
inline vec3 reflect( const vec3& direction, const vec3& normal ) {
    return ( normal - direction ) * ( 2.0f * dot( direction, normal ) );
}
inline vec3 lerp( const vec3& a, const vec3& b, f32 t ) {
    return ( 1.0f - t ) * a + b * t;
}
inline vec3 clampedLerp( const vec3& a, const vec3& b, f32 t ) {
    return lerp( a, b, clamp( t, 0.0f, 1.0f ) );
}

vec4 operator+( const vec4& lhs, const vec4& rhs );
vec4 operator-( const vec4& lhs, const vec4& rhs );
vec4 operator*( const vec4& lhs, f32 rhs );
vec4 operator*( f32 lhs, const vec4& rhs );
vec4 operator/( const vec4& lhs, f32 rhs );
// normalize vector 
inline vec4 normalize( const vec4& v );
// clamp vector length up to max, inclusive
inline vec4 clamp( const vec4& v, f32 max );
// compare two vectors
inline bool cmp( const vec4& lhs, const vec4& rhs );
// dot product
inline f32 dot( const vec4& lhs, const vec4& rhs );
// scale vectors component-wise
inline vec4 mul( const vec4& lhs, const vec4& rhs );
// square magnitude
inline f32 sqrMag( const vec4& v );
// magnitude
inline f32 mag( const vec4& v );
/// @brief Linear interpolation
/// @param a minimum value
/// @param b maximum value
/// @param t fraction
/// @return blend between a and b, based on fraction t
inline vec4 lerp( const vec4& a, const vec4& b, f32 t );
/// @brief Linear interpolation, t is clamped between 0.0-1.0
/// @param a minimum value
/// @param b maximum value
/// @param t 0.0-1.0 fraction
/// @return blend between a and b, based on fraction t
inline vec4 clampedLerp( const vec4& a, const vec4& b, f32 t );
// 4-component vector
struct vec4 {
    union {
        struct { f32 x, y, z, w; };
        struct { f32 r, g, b, a; };
    };

    vec4() {}
    vec4( f32 scalar ) : x(scalar), y(scalar), z(scalar), w(scalar) {}
    vec4( const vec2& v ) : x(v.x), y(v.y), z(0.0f), w(0.0f) {}
    vec4( const vec3& v ) : x(v.x), y(v.y), z(v.z), w(1.0f) {}
    vec4( f32 x, f32 y, f32 z, f32 w ) : x(x), y(y), z(z), w(w) {}

    // get pointer to struct as f32 
    f32* ptr() { return (f32*)this; }

    vec4& operator-() { return *this *= -1.0f; }
    bool operator==( const vec4& rhs ) { return cmp( *this, rhs ); }
    bool operator!=( const vec4& rhs ) { return !(*this == rhs); }
    vec4& operator+=( const vec4& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( this->w, this->z, this->y, this->x );
        __m128 b = _mm_set_ps( rhs.w, rhs.z, rhs.y, rhs.x );
        _mm_storeu_ps( this->ptr(), _mm_add_ps( a, b ) );
        return *this;
    }
    vec4& operator-=( const vec4& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( this->w, this->z, this->y, this->x );
        __m128 b = _mm_set_ps( rhs.w, rhs.z, rhs.y, rhs.x );
        _mm_storeu_ps( this->ptr(), _mm_sub_ps( a, b ) );
        return *this;
    }
    vec4& operator*=( const f32& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( this->w, this->z, this->y, this->x );
        __m128 b = _mm_set_ps( rhs, rhs, rhs, rhs );
        _mm_storeu_ps( this->ptr(), _mm_mul_ps( a, b ) );

        return *this;
    }
    vec4& operator/=( const f32& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( this->w, this->z, this->y, this->x );
        __m128 b = _mm_set_ps( rhs, rhs, rhs, rhs );
        _mm_storeu_ps( this->ptr(), _mm_div_ps( a, b ) );

        return *this;
    }

    // ( 1.0f, 1.0f, 1.0f, 1.0f ) 
    static vec4 one()   { return {  1.0f,  1.0f,  1.0f,  1.0f }; }
    // ( 0.0f, 0.0f, 0.0f, 0.0f ) 
    static vec4 zero()  { return {  0.0f,  0.0f,  0.0f,  0.0f }; }
};
inline vec4 operator+( const vec4& lhs, const vec4& rhs ) {
    return vec4(lhs) += rhs;
}
inline vec4 operator-( const vec4& lhs, const vec4& rhs ) {
    return vec4(lhs) -= rhs;
}
inline vec4 operator*( const vec4& lhs, f32 rhs ) {
    return vec4(lhs) *= rhs;
}
inline vec4 operator*( f32 lhs, const vec4& rhs ) {
    return vec4(rhs) *= lhs;
}
inline vec4 operator/( const vec4& lhs, f32 rhs ) {
    return vec4(lhs) /= rhs;
}
inline vec4 normalize( const vec4& v ) {
    f32 m = mag( v );
    if( m != 0.0f ) {
        return v / m;
    } else {
        return vec4::zero();
    }
}
inline bool cmp( const vec4& lhs, const vec4& rhs ) {
    return sqrMag(lhs - rhs) < VEC_CMP_THRESHOLD;
}
inline vec4 clamp( const vec4& v, f32 max ) {
    f32 magnitude = mag(v);
    if( magnitude > max ) {
        vec4 result = v / magnitude;
        result *= max;
        return result;
    } else {
        return v;
    }
}
inline f32 dot( const vec4& lhs, const vec4& rhs ) {
    // NOTE(alicia): SSE
    // TODO(alicia): further optimizations? FMA?

    __m128 a = _mm_set_ps( lhs.x, lhs.y, lhs.z, lhs.w );
    __m128 b = _mm_set_ps( rhs.x, rhs.y, rhs.z, rhs.w );
    f32 result[4];
    _mm_storeu_ps( result, _mm_mul_ps( a, b ) );

    return result[0] + result[1] + result[2] + result[3];
}
inline vec4 mul( const vec4& lhs, const vec4& rhs ) {
    // NOTE(alicia): SSE

    __m128 a = _mm_set_ps( lhs.w, lhs.z, lhs.y, lhs.x );
    __m128 b = _mm_set_ps( rhs.w, rhs.z, rhs.y, rhs.x );
    vec4 result = {};
    _mm_storeu_ps( result.ptr(), _mm_mul_ps( a, b ) );

    return result;
}
inline f32 sqrMag( const vec4& v ) {
    // NOTE(alicia): SSE
    // TODO(alicia): further optimizations? FMA?

    __m128 a = _mm_set_ps( v.x, v.y, v.z, v.w );
    f32 result[4];
    _mm_storeu_ps( result, _mm_mul_ps( a, a ) );

    return result[0] + result[1] + result[2] + result[3];
}
inline f32 mag( const vec4& v ) {
    return sqrt( sqrMag( v ) );
}
inline vec4 lerp( const vec4& a, const vec4& b, f32 t ) {
    return ( 1.0f - t ) * a + b * t;
}
inline vec4 clampedLerp( const vec4& a, const vec4& b, f32 t ) {
    return lerp( a, b, clamp( t, 0.0f, 1.0f ) );
}

quat operator+( const quat& lhs, const quat& rhs );
quat operator-( const quat& lhs, const quat& rhs );
quat operator*( const quat& lhs, const quat& rhs );
vec3 operator*( const quat& lhs, const vec3& rhs );
quat operator*( const quat& lhs, f32 rhs );
quat operator*( f32 lhs, const quat& rhs );
quat operator/( const quat& lhs, f32 rhs );
// normalize quaternion 
inline quat normalize( const quat& q );
// compare two quaternions
inline bool cmp( const quat& lhs, const quat& rhs );
// dot product
inline f32 dot( const quat& lhs, const quat& rhs );
// square magnitude
inline f32 sqrMag( const quat& q );
// magnitude
inline f32 mag( const quat& q );
// conjugate of quaternion
inline quat conjugate( const quat& q );
// invert quaternion
inline quat inverse( const quat& q );
/// @brief Linear interpolation
/// @param a minimum value
/// @param b maximum value
/// @param t fraction
/// @return blend between a and b, based on fraction t
inline quat slerp( const quat& a, const quat& b, f32 t );
/// @brief Linear interpolation
/// @param a minimum value
/// @param b maximum value
/// @param t fraction
/// @return blend between a and b, based on fraction t
inline quat clampedSlerp( const quat& a, const quat& b, f32 t );
// quaternion rotation
struct quat {
    union {
        struct { f32 w, x, y, z; };
        struct { f32 a, b, c, d; };
    };

    quat() {}
    quat( const vec4& v ) : w(v.w), x(v.x), y(v.y), z(v.z) {}
    quat( f32 w, f32 x, f32 y, f32 z ) : w(w), x(x), y(y), z(z) {}

    // get pointer to struct as f32 
    f32* ptr() { return (f32*)this; }

    quat& operator-() { return *this *= -1.0f; }
    bool operator==( const quat& rhs ) { return cmp( *this, rhs ); }
    bool operator!=( const quat& rhs ) { return !(*this == rhs); }
    quat& operator+=( const quat& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( this->z, this->y, this->x, this->w );
        __m128 b = _mm_set_ps( rhs.z, rhs.y, rhs.x, rhs.w );
        _mm_storeu_ps( this->ptr(), _mm_add_ps( a, b ) );
        return *this;
    }
    quat& operator-=( const quat& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( this->z, this->y, this->x, this->w );
        __m128 b = _mm_set_ps( rhs.z, rhs.y, rhs.x, rhs.w );
        _mm_storeu_ps( this->ptr(), _mm_sub_ps( a, b ) );
        return *this;
    }
    quat& operator*=( const quat& rhs ) {
        // TODO(alicia): SIMD
        this->w = ( this->w * rhs.w ) - ( (this->x * rhs.x) + (this->y * rhs.y) + (this->z * rhs.z) );
        this->x = ( this->w * rhs.x ) + ( (rhs.w * this->x) + ( (this->y * rhs.z) - (this->z * rhs.y) ) );
        this->y = ( this->w * rhs.y ) + ( (rhs.w * this->y) + ( (this->z * rhs.x) - (this->x * rhs.z) ) );
        this->x = ( this->w * rhs.z ) + ( (rhs.w * this->z) + ( (this->x * rhs.y) - (this->y * rhs.x) ) );
        return *this;
    }
    quat& operator*=( const f32& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( this->z, this->y, this->x, this->w );
        __m128 b = _mm_set_ps( rhs, rhs, rhs, rhs );
        _mm_storeu_ps( this->ptr(), _mm_mul_ps( a, b ) );
        return *this;
    }
    quat& operator/=( const f32& rhs ) {
        // NOTE(alicia): SSE

        __m128 a = _mm_set_ps( this->z, this->y, this->x, this->w );
        __m128 b = _mm_set_ps( rhs, rhs, rhs, rhs );
        _mm_storeu_ps( this->ptr(), _mm_div_ps( a, b ) );
        return *this;
    }

    /// @brief convert to angle-axis
    /// @param angle [out] angle
    /// @param axis [out] axis
    void toAngleAxis( f32& angle, vec3& axis ) {
        angle = 2.0f * acos(w);
        f32 invW2sqrt = sqrt(1.0f - sqr(w));

        // NOTE(alicia): SSE
        __m128 a = _mm_set_ps( 1.0f, z, y, x );
        __m128 b = _mm_set1_ps( invW2sqrt );
        f32 result[3];
        _mm_storeu_ps( result, _mm_div_ps( a, b ) );

        axis = {
            result[0],
            result[1],
            result[2]
        };
    }
    /// @brief convert to euler angles
    /// @param euler [out] euler angles
    void toEuler( vec3& euler ) {
        // TODO(alicia): SIMD!

        euler.x = atan2( 2.0f * ( w * x + y * z ), 1.0f - 2.0f * ( x * x + y * y ) );
        euler.y = asinNoNaN( 2.0f * ( ( w * y ) - ( z * x ) ) );
        euler.z = atan2( 2.0f * ( w * z + x * y ), 1.0f - 2.0f * ( y * y + z * z ) );
    }
    // identity quaternion
    static quat identity() { return { 1.0f, 0.0f, 0.0f, 0.0f }; }
    // construct quaternion from angle-axis
    static quat angleAxis( f32 theta, const vec3& axis ) {
        f32 halfTheta = theta / 2.0f;
        f32 s = sin( halfTheta );

        // NOTE(alicia): SSE
        __m128 a = _mm_set_ps( 0.0f, axis.z, axis.y, axis.x );
        __m128 b = _mm_set1_ps( s );
        f32 result[3];
        _mm_storeu_ps( result, _mm_mul_ps(a, b) );

        return {
            cos( halfTheta ),
            result[0],
            result[1],
            result[2]
        };
    }
    // construct quaternion from euler angles
    static quat euler( f32 x, f32 y, f32 z ) {
        // TODO(alicia): SIMD?

        f32 halfX = x / 2.0f;
        f32 halfY = y / 2.0f;
        f32 halfZ = z / 2.0f;

        f32 xSin = sin( halfX );
        f32 xCos = cos( halfX );

        f32 ySin = sin( halfY );
        f32 yCos = cos( halfY );

        f32 zSin = sin( halfZ );
        f32 zCos = cos( halfZ );

        return {
            ( xCos * yCos * zCos ) + ( xSin * ySin * zSin ),

            ( xSin * yCos * zCos ) + ( xCos * ySin * zSin ),
            ( xCos * ySin * zCos ) + ( xSin * yCos * zSin ),
            ( xCos * yCos * zSin ) + ( xSin * ySin * zCos )
        };
    }
    // construct quaternion from euler angles
    static quat euler( const vec3& euler ) {
        return quat::euler( euler.x, euler.y, euler.z );
    }
};
inline quat operator+( const quat& lhs, const quat& rhs ) {
    return quat(lhs) += rhs;
}
inline quat operator-( const quat& lhs, const quat& rhs ) {
    return quat(lhs) -= rhs;
}
quat operator*( const quat& lhs, const quat& rhs ) {
    return quat(lhs) *= rhs;
}
vec3 operator*( const quat& lhs, const vec3& rhs ) {
    quat p = { 0.0f, rhs.x, rhs.y, rhs.z };
    quat result = lhs * p * conjugate( lhs );
    return { result.x, result.y, result.z };
}
inline quat operator*( const quat& lhs, f32 rhs ) {
    return quat(lhs) *= rhs;
}
inline quat operator*( f32 lhs, const quat& rhs ) {
    return quat(rhs) *= lhs;
}
inline quat operator/( const quat& lhs, f32 rhs ) {
    return quat(lhs) /= rhs;
}
inline quat normalize( const quat& q ) {
    f32 m = mag( q );
    if( m != 0.0f ) {
        return q / m;
    } else {
        return {};
    }
}
inline bool cmp( const quat& lhs, const quat& rhs ) {
    return sqrMag(lhs - rhs) < VEC_CMP_THRESHOLD;
}
inline f32 dot( const quat& lhs, const quat& rhs ) {
    // NOTE(alicia): SSE
    // TODO(alicia): further optimizations? FMA?

    __m128 a = _mm_set_ps( lhs.z, lhs.y, lhs.x, lhs.w );
    __m128 b = _mm_set_ps( rhs.z, rhs.y, rhs.x, rhs.w );
    f32 result[4];
    _mm_storeu_ps( result, _mm_mul_ps( a, b ) );

    return result[0] + result[1] + result[2] + result[3];
}
inline f32 sqrMag( const quat& q ) {
    // NOTE(alicia): SSE
    // TODO(alicia): further optimizations? FMA?

    __m128 a = _mm_set_ps( q.z, q.y, q.x, q.w );
    f32 result[4];
    _mm_storeu_ps( result, _mm_mul_ps( a, a ) );

    return result[0] + result[1] + result[2] + result[3];
}
inline f32 mag( const quat& q ) {
    return sqrt( sqrMag( q ) );
}
inline quat conjugate( const quat& q ) {
    quat r = q;
    r   = -r;
    r.w = -r.w;
    return r;
}
inline quat inverse( const quat& q ) {
    return conjugate( q ) / sqrMag( q );
}
inline quat slerp( const quat& a, const quat& b, f32 t ) {
    f32 dotProd = dot( a, b );
    f32 lambda  = t / 2.0f;
    f32 theta   = abs( acos( dotProd ) );

    f32 thetaSin = sin( theta );
    f32 coeff1 = sin( ( 1.0f - theta ) * theta ) / thetaSin;
    f32 coeff2 = sin( lambda * theta ) / thetaSin;

    // NOTE(alicia): SSE

    __m128 a1 = _mm_set_ps( a.z, a.y, a.x, a.w );
    __m128 a2 = _mm_set_ps( b.z, b.y, b.x, b.w );
    __m128 b1 = _mm_set1_ps( coeff1 );
    __m128 b2 = _mm_set1_ps( coeff2 );

    __m128 c1 = _mm_mul_ps( a1, b1 );
    __m128 c2 = _mm_mul_ps( a2, b2 );

    quat result = {};
    _mm_storeu_ps( result.ptr(), _mm_add_ps( c1, c2 ) );

    return normalize( result );
}
inline quat clampedSlerp( const quat& a, const quat& b, f32 t ) {
    return slerp( a, b, clamp( t, 0.0f, 1.0f ) );
}

} // namespace smath