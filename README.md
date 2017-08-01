# chicken-fftw3
Thin CHICKEN wrapper for fftw3

## Usage

There are two submodules containing procedures working on f32vectors and f64vectors

```scheme
(require-library fftw)
; Load the f32 or f64 interface
(import fftw.f32)
; or
(import fftw.f64)
; You can also load both at the same time
(import (prefix fftw.f32 f32:))
(import (prefix fftw.f64 f64:))
```
