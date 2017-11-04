(include "fftw.scm")

(module fftw-f64 = fftw-impl
  (import scheme chicken foreign)
  (use srfi-4)

  (define fvector-length f64vector-length)

  (define-foreign-type c64vector (nonnull-scheme-pointer "fftw_complex")
    f64vector->blob/shared)

  (define c-fftw-destroy-plan
    (foreign-lambda void fftw_destroy_plan c-pointer))
  (define c-fftw-alignment-of
    (foreign-lambda int fftw_alignment_of scheme-pointer))

  (define c-fftw-execute-c2c
    (foreign-lambda void fftw_execute_dft c-pointer c64vector c64vector))
  (define c-fftw-execute-r2c
    (foreign-lambda void fftw_execute_dft_r2c c-pointer f64vector c64vector))
  (define c-fftw-execute-c2r
    (foreign-lambda void fftw_execute_dft_c2r c-pointer c64vector f64vector))
  (define c-fftw-execute-r2r
    (foreign-lambda void fftw_execute_r2r c-pointer f64vector f64vector))

  (define c-fftw-plan-c2c
    (foreign-lambda c-pointer fftw_plan_dft int s32vector c64vector c64vector int int))
  (define c-fftw-plan-r2c
    (foreign-lambda c-pointer fftw_plan_dft_r2c int s32vector f64vector c64vector int))
  (define c-fftw-plan-c2r
    (foreign-lambda c-pointer fftw_plan_dft_c2r int s32vector c64vector f64vector int))
  (define c-fftw-plan-r2r
    (foreign-lambda c-pointer fftw_plan_r2r int s32vector f64vector f64vector c-pointer int)))
