(include "fftw.scm")

(module fftw-f32 = fftw-impl
  (import scheme chicken foreign)
  (use srfi-4)

  (define fvector-length f32vector-length)

  (define-foreign-type c32vector (nonnull-scheme-pointer "fftwf_complex")
    f32vector->blob/shared)

  (define c-fftw-execute
    (foreign-lambda void fftwf_execute plan))
  (define c-fftw-destroy-plan
    (foreign-lambda void fftwf_destroy_plan plan))

  (define c-fftw-plan-c2c
    (foreign-lambda plan fftwf_plan_dft int s32vector c32vector c32vector int int))
  (define c-fftw-plan-r2c
    (foreign-lambda plan fftwf_plan_dft_r2c int s32vector f32vector c32vector int))
  (define c-fftw-plan-c2r
    (foreign-lambda plan fftwf_plan_dft_c2r int s32vector c32vector f32vector int))
  (define c-fftw-plan-r2r
    (foreign-lambda plan fftwf_plan_r2r int s32vector f32vector f32vector c-pointer int)))
