(module fftw
  (fft! rfft! ifft! irfft!
   fft* rfft* ifft* irfft*)
  (import scheme chicken foreign)

(use srfi-4)

(foreign-declare "#include <fftw3.h>")

(define-foreign-type plan c-pointer)

(define-foreign-type cvector (nonnull-scheme-pointer "fftw_complex")
  (lambda (x)
    (let* ((blob (f64vector->blob/shared x))
           (elem (fx/ (##sys#size blob) 8)))
      (when (odd? elem) (error "odd!" elem))
      blob)))

(define-foreign-variable fftw-estimate   int "FFTW_ESTIMATE")
(define-foreign-variable fftw-measure    int "FFTW_MEASURE")
(define-foreign-variable fftw-patient    int "FFTW_PATIENT")
(define-foreign-variable fftw-exhaustive int "FFTW_EXHAUSTIVE")

(define c-fftw-execute
  (foreign-lambda void fftw_execute plan))
(define c-fftw-destroy-plan
  (foreign-lambda void fftw_destroy_plan plan))
(define c-fftw-cleanup
  (foreign-lambda void fftw_cleanup))

; complex->complex
(define c-fftw-plan-c2c
  (foreign-lambda plan fftw_plan_dft int s32vector cvector cvector int int))

; real->complex
(define c-fftw-plan-r2c
  (foreign-lambda plan fftw_plan_dft_r2c int s32vector f64vector cvector int))

; complex->real
(define c-fftw-plan-c2r
  (foreign-lambda plan fftw_plan_dft_c2r int s32vector cvector f64vector int))

(define-syntax wrap-complex-trasform
  (er-macro-transformer
    (lambda (x r c)
      (let* ((name     (cadr (strip-syntax x)))
             (forward? (caddr x))
             (execute? (cadddr x)))
        `(define ,(if execute? (symbol-append name '!) name)
           (lambda (dim in out #!optional flags)
             (let* ((rank       (length dim))
                    (total-dim  (foldl fx* 1 dim))
                    (min-length (fx* 2 total-dim)))
               (unless (fx> total-dim 1)
                 (error "dim"))
               (unless (fx>= (f64vector-length out) min-length)
                 (error "out-size"))
               (unless (fx>= (f64vector-length in) min-length)
                 (error "in-size"))
               (let ((plan
                       (c-fftw-plan-c2c rank (list->s32vector dim) in out
                                        ,(if forward? -1 +1)
                                        (or flags fftw-estimate))))
                 ,(if execute?
                      `(begin
                         (c-fftw-execute plan)
                         (c-fftw-destroy-plan plan))
                      `(set-finalizer! plan c-fftw-destroy-plan))))))))))

(wrap-complex-trasform fft   #t #t)
(wrap-complex-trasform ifft  #f #t)
(wrap-complex-trasform fft*  #t #f)
(wrap-complex-trasform ifft* #f #f)

(define-syntax wrap-real-trasform
  (er-macro-transformer
    (lambda (x r c)
      (let* ((name     (cadr (strip-syntax x)))
             (forward? (caddr x))
             (execute? (cadddr x))
             (re-vec   (if forward? 'in 'out))
             (im-vec   (if forward? 'out 'in)))
        `(define ,(if execute? (symbol-append name '!) name)
           (lambda (dim in out #!optional flags)
             (let* ((rank      (length dim))
                    (total-dim (foldl fx* 1 dim))
                    (re-length (f64vector-length ,re-vec))
                    (im-length (f64vector-length ,im-vec)))
               (unless (fx> total-dim 1)
                 (error "dim"))
               ; The Hermitian symmetry allows us to spare some space by
               ; only storing half of the spectrum
               ; (N/2 + 1) * 2 = N + 2
               (unless (fx>= im-length (fx+ 2 total-dim))
                 (error "out-size"))
               ; The dimension refers to the size of the real vector
               (unless (fx>= re-length total-dim)
                 (error "in-size"))
               (let ((plan
                       ,(if forward?
                            ; real->complex
                            `(c-fftw-plan-r2c rank (list->s32vector dim)
                                              ,re-vec ,im-vec
                                              (or flags fftw-estimate))
                            ; complex->real
                            `(c-fftw-plan-c2r rank (list->s32vector dim)
                                              ,im-vec ,re-vec
                                              (or flags fftw-estimate)))))
                 ,(if execute?
                      `(begin
                         (c-fftw-execute plan)
                         (c-fftw-destroy-plan plan))
                      `(set-finalizer! plan c-fftw-destroy-plan))))))))))

(wrap-real-trasform rfft   #t #t)
(wrap-real-trasform irfft  #f #t)
(wrap-real-trasform rfft*  #t #f)
(wrap-real-trasform irfft* #f #f)
)
