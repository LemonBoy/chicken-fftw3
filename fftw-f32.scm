(include "fftw.scm")

(module fftw-f32 = fftw-impl
  (import scheme chicken foreign)
  (use srfi-4)

  (define fvector-length f32vector-length)

  (define-foreign-type c32vector (nonnull-scheme-pointer "fftwf_complex")
    f32vector->blob/shared)

  (define c-fftw-destroy-plan
    (foreign-lambda void fftwf_destroy_plan c-pointer))
  (define c-fftw-alignment-of
    (foreign-lambda int fftwf_alignment_of scheme-pointer))

  (define c-fftw-execute-c2c
    (foreign-lambda void fftwf_execute_dft c-pointer c32vector c32vector))
  (define c-fftw-execute-r2c
    (foreign-lambda void fftwf_execute_dft_r2c c-pointer f32vector c32vector))
  (define c-fftw-execute-c2r
    (foreign-lambda void fftwf_execute_dft_c2r c-pointer c32vector f32vector))
  (define c-fftw-execute-r2r
    (foreign-lambda void fftwf_execute_r2r c-pointer f32vector f32vector))

  (define c-fftw-plan-c2c
    (foreign-lambda c-pointer fftwf_plan_dft int s32vector c32vector c32vector int int))
  (define c-fftw-plan-r2c
    (foreign-lambda c-pointer fftwf_plan_dft_r2c int s32vector f32vector c32vector int))
  (define c-fftw-plan-c2r
    (foreign-lambda c-pointer fftwf_plan_dft_c2r int s32vector c32vector f32vector int))
  (define c-fftw-plan-r2r
    (foreign-lambda c-pointer fftwf_plan_r2r int s32vector f32vector f32vector c-pointer int)))
