(include "fftw.scm")

(module fftw-f64 = fftw-impl
  (import scheme chicken foreign)
  (use srfi-4)

  (define fvector-length f64vector-length)

  (define-foreign-type c64vector (nonnull-scheme-pointer "fftw_complex")
    f64vector->blob/shared)

  (define c-fftw-execute
    (foreign-lambda void fftw_execute plan))
  (define c-fftw-destroy-plan
    (foreign-lambda void fftw_destroy_plan plan))

  (define c-fftw-execute-c2c
    (foreign-lambda void fftw_execute_dft plan c64vector c64vector))
  (define c-fftw-execute-r2c
    (foreign-lambda void fftw_execute_dft_r2c plan f64vector c64vector))
  (define c-fftw-execute-c2r
    (foreign-lambda void fftw_execute_dft_c2r plan c64vector f64vector))
  (define c-fftw-execute-r2r
    (foreign-lambda void fftw_execute_r2r plan f64vector f64vector))

  (define c-fftw-plan-c2c
    (foreign-lambda plan fftw_plan_dft int s32vector c64vector c64vector int int))
  (define c-fftw-plan-r2c
    (foreign-lambda plan fftw_plan_dft_r2c int s32vector f64vector c64vector int))
  (define c-fftw-plan-c2r
    (foreign-lambda plan fftw_plan_dft_c2r int s32vector c64vector f64vector int))
  (define c-fftw-plan-r2r
    (foreign-lambda plan fftw_plan_r2r int s32vector f64vector f64vector c-pointer int)))
