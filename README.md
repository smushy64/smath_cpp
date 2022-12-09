# smath (c++)
smath is a header-only math library for my personal C++ projects.

- [src](src)
- [to-do list](TODO.md)

## About

This is **not** a general-purpose math library. It is specifically written for my purposes only.

Feel free to use it in your projects if it fits your needs, no attribution needed :)

Currently uses only SSE instructions for vector arithmetic.

Only tested with c++20, gcc compiler.

## How to use

Download [smath.hpp](src/smath.hpp) and include it in your project.

Requires
- `cstdint`     should already have it as every compiler includes this header
- `immintrin.h` should already have it as every compiler includes this header

## Example
```cpp

#include "smath.hpp"

int main(void) {
    smath::vec2 v1 = smath::vec2::left() * 2.0f;
    smath::vec2 v2 = smath::vec2::up();
    smath::vec2 v3 = smath::vec2( 0.25f, 0.5f );
    
    smath::vec2 addition = v1 + v2;
    // result = { -2.0f, 1.0f }

    smath::vec2 subtraction = v1 - v3;
    // result = { -2.25f, -0.5f }

    f32 dotProd = smath::dot( v1, v3 );
    // result = -0.5f
    
    return 0;
}

```