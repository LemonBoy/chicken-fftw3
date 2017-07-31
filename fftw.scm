(module fftw
  (fft! rfft! ifft! irfft! dct! dst!
   plan-fft plan-rfft plan-ifft plan-irfft plan-dct plan-dst
   execute-plan
   fftw-estimate fftw-measure fftw-patient fftw-exhaustive fftw-preserve-input)
  (import scheme chicken foreign)

(use srfi-4)

(foreign-declare "#include <fftw3.h>")

(define-foreign-type plan c-pointer)

(define-foreign-type cvector (nonnull-scheme-pointer "fftw_complex")
  (lambda (x)
    (let* ((blob (f64vector->blob/shared x))
           (elem (fx/ (##sys#size blob) 8)))
      (when (odd? elem) (error "invalid f64vector size - must be even"))
      blob)))

(define-foreign-variable c-fftw-estimate int       "FFTW_ESTIMATE")
(define-foreign-variable c-fftw-measure int        "FFTW_MEASURE")
(define-foreign-variable c-fftw-patient int        "FFTW_PATIENT")
(define-foreign-variable c-fftw-exhaustive int     "FFTW_EXHAUSTIVE")
(define-foreign-variable c-fftw-preserve-input int "FFTW_PRESERVE_INPUT")

(define-foreign-variable fftw-dctI  int "FFTW_REDFT00")
(define-foreign-variable fftw-dstI  int "FFTW_RODFT00")

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

; real->real
(define c-fftw-plan-r2r
  (foreign-lambda plan fftw_plan_r2r int s32vector f64vector f64vector c-pointer int))

; flags
(define fftw-estimate       c-fftw-estimate)
(define fftw-measure        c-fftw-measure)
(define fftw-patient        c-fftw-patient)
(define fftw-exhaustive     c-fftw-exhaustive)
(define fftw-preserve-input c-fftw-preserve-input)

(define execute-plan c-fftw-execute)

(define (kind->enum base x)
  (case x
    ((I)   base)
    ((II)  (fx+ base 1))
    ((III) (fx+ base 2))
    ((IV)  (fx+ base 3))
    (else
      (error "invalid transform kind - must be I II III or IV" x))))

(define-syntax wrap-real-eo-transform
  (er-macro-transformer
    (lambda (x r c)
      (let* ((name  (cadr (strip-syntax x)))
             (even? (caddr x)))
        `(begin
           (define ,(symbol-append 'plan- name)
             (lambda (kind in out #!optional dim flags)
               (let* ((flags      (or flags c-fftw-estimate))
                      (dim        (or dim (list (f64vector-length in))))
                      (rank       (length dim))
                      (total-dim  (foldl fx* 1 dim))
                      (base       ,(if even? 'fftw-dctI 'fftw-dstI))
                      (kind       (kind->enum base kind)))
                 (unless (fx> total-dim 1)
                   (error "invalid transform size" dim))
                 (unless (fx>= (f64vector-length out) total-dim)
                   (error "output vector length is too short"))
                 (unless (fx>= (f64vector-length in) total-dim)
                   (error "input vector length is too short"))
                 (set-finalizer!
                   (c-fftw-plan-r2r rank (list->s32vector dim) in out
                                    #$(make-s32vector rank kind)
                                    flags)
                   c-fftw-destroy-plan))))
           (define ,(symbol-append name '!)
             (lambda (kind in out #!optional dim)
               (c-fftw-execute
                 (,(symbol-append 'plan- name) kind in out dim #f)))))))))

(wrap-real-eo-transform dct #t)
(wrap-real-eo-transform dst #f)

(define-syntax wrap-complex-trasform
  (er-macro-transformer
    (lambda (x r c)
      (let* ((name     (cadr (strip-syntax x)))
             (forward? (caddr x)))
        `(begin
           (define ,(symbol-append 'plan- name)
             (lambda (in out #!optional dim flags)
               (let* ((flags      (or flags c-fftw-estimate))
                      (dim        (or dim (list (f64vector-length in))))
                      (rank       (length dim))
                      (total-dim  (foldl fx* 1 dim))
                      (min-length (fx* 2 total-dim)))
                 (unless (fx> total-dim 1)
                   (error "invalid transform size" dim))
                 (unless (fx>= (f64vector-length out) min-length)
                   (error "output vector length is too short"))
                 (unless (fx>= (f64vector-length in) min-length)
                   (error "input vector length is too short"))
                 (set-finalizer!
                   (c-fftw-plan-c2c rank (list->s32vector dim) in out
                                    ,(if forward? -1 +1) flags)
                   c-fftw-destroy-plan))))
           (define ,(symbol-append name '!)
             (lambda (in out #!optional dim)
               (c-fftw-execute
                 (,(symbol-append 'plan- name) in out dim #f)))))))))

(wrap-complex-trasform fft  #t)
(wrap-complex-trasform ifft #f)

(define-syntax wrap-real-trasform
  (er-macro-transformer
    (lambda (x r c)
      (let* ((name     (cadr (strip-syntax x)))
             (forward? (caddr x))
             (re-vec   (if forward? 'in 'out))
             (cpl-vec  (if forward? 'out 'in)))
        `(begin
           (define ,(symbol-append 'plan- name)
             (lambda (in out #!optional dim flags)
               (let* ((flags      (or flags c-fftw-estimate))
                      (dim        (or dim (list (f64vector-length in))))
                      (rank       (length dim))
                      (total-dim  (foldl fx* 1 dim))
                      (re-length  (f64vector-length ,re-vec))
                      (cpl-length (f64vector-length ,cpl-vec)))
                 (unless (fx> total-dim 1)
                   (error "invalid transform size" dim))
                 ; The Hermitian symmetry allows us to spare some space by
                 ; only storing half of the spectrum
                 ; (N/2 + 1) * 2 = N + 2
                 (unless (fx>= cpl-length (fx+ 2 total-dim))
                   (error "complex vector length is too short"))
                 ; The dimension refers to the size of the real vector
                 (unless (fx>= re-length total-dim)
                   (error "real vector length is too short"))
                 (set-finalizer!
                   (,(if forward? 'c-fftw-plan-r2c 'c-fftw-plan-c2r)
                     rank (list->s32vector dim) in out flags)
                   c-fftw-destroy-plan))))
           (define ,(symbol-append name '!)
             (lambda (in out dim)
               (c-fftw-execute
                 (,(symbol-append 'plan- name) in out dim #f)))))))))

(wrap-real-trasform rfft  #t)
(wrap-real-trasform irfft #f)
)
