(foreign-declare "#include <fftw3.h>")

(define-foreign-variable c-fftw-estimate int "FFTW_ESTIMATE")
(define-foreign-variable c-fftw-measure int "FFTW_MEASURE")
(define-foreign-variable c-fftw-patient int "FFTW_PATIENT")
(define-foreign-variable c-fftw-exhaustive int "FFTW_EXHAUSTIVE")
(define-foreign-variable c-fftw-preserve-input int "FFTW_PRESERVE_INPUT")
(define-foreign-variable c-fftw-unaligned int "FFTW_UNALIGNED")
(define-foreign-variable c-fftw-conserve-memory int "FFTW_CONSERVE_MEMORY")

(define-foreign-variable fftw-dctI int "FFTW_REDFT00")
(define-foreign-variable fftw-dstI int "FFTW_RODFT00")

(define-foreign-type plan c-pointer)

(functor
  (fftw-impl
    (I (c-fftw-destroy-plan
        c-fftw-execute c-fftw-execute-c2c c-fftw-execute-c2r c-fftw-execute-r2c
        c-fftw-execute-r2r
        c-fftw-plan-c2c c-fftw-plan-r2c c-fftw-plan-c2r c-fftw-plan-r2r)))
  (fft! rfft! ifft! irfft! dct! dst!
   plan-fft plan-rfft plan-ifft plan-irfft plan-dct plan-dst
   execute-plan
   fftw-estimate fftw-measure fftw-patient fftw-exhaustive fftw-preserve-input
   fftw-unaligned fftw-conserve-memory)

  (import chicken scheme foreign I)
  (use srfi-4)

  (define (kind->enum base x)
    (case x
      ((I)   base)
      ((II)  (fx+ base 1))
      ((III) (fx+ base 2))
      ((IV)  (fx+ base 3))
      (else
        (error "invalid transform kind - must be I II III or IV" x))))

  ; exported flags

  (define fftw-estimate c-fftw-estimate)
  (define fftw-measure c-fftw-measure)
  (define fftw-patient c-fftw-patient)
  (define fftw-exhaustive c-fftw-exhaustive)
  (define fftw-preserve-input c-fftw-preserve-input)
  (define fftw-unaligned c-fftw-unaligned)
  (define fftw-conserve-memory c-fftw-conserve-memory)

  ; exported procedures

  (define execute-plan c-fftw-execute)

  (define-syntax wrap-real-eo-transform
    (er-macro-transformer
      (lambda (x r c)
        (let* ((name  (cadr (strip-syntax x)))
               (even? (caddr x)))
          `(begin
             (define ,(symbol-append 'plan- name)
               (lambda (kind in out #!optional dim flags)
                 (let* ((flags      (or flags c-fftw-estimate))
                        (dim        (or dim (list (fvector-length in))))
                        (rank       (length dim))
                        (total-dim  (foldl fx* 1 dim))
                        (base       ,(if even? 'fftw-dctI 'fftw-dstI))
                        (kind       (kind->enum base kind)))
                   (unless (fx> total-dim 1)
                     (error "invalid transform size" dim))
                   (unless (fx>= (fvector-length out) total-dim)
                     (error "output vector length is too short"))
                   (unless (fx>= (fvector-length in) total-dim)
                     (error "input vector length is too short"))
                   (set-finalizer!
                     (c-fftw-plan-r2r rank (list->s32vector dim) in out
                                      #$(make-s32vector rank kind)
                                      flags)
                     c-fftw-destroy-plan))))
             (define ,(symbol-append name '!)
               (lambda (kind in out #!optional dim)
                 (c-fftw-execute-r2r
                   (,(symbol-append 'plan- name) kind in out dim #f)
                   in out))))))))

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
                        (dim        (or dim (list (fvector-length in))))
                        (rank       (length dim))
                        (total-dim  (foldl fx* 1 dim))
                        (min-length (fx* 2 total-dim)))
                   (unless (fx> total-dim 1)
                     (error "invalid transform size" dim))
                   (unless (fx>= (fvector-length out) min-length)
                     (error "output vector length is too short"))
                   (unless (fx>= (fvector-length in) min-length)
                     (error "input vector length is too short"))
                   (set-finalizer!
                     (c-fftw-plan-c2c rank (list->s32vector dim) in out
                                      ,(if forward? -1 +1) flags)
                     c-fftw-destroy-plan))))
             (define ,(symbol-append name '!)
               (lambda (in out #!optional dim)
                 (c-fftw-execute-c2c
                   (,(symbol-append 'plan- name) in out dim #f)
                   in out))))))))

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
                        (dim        (or dim (list (fvector-length in))))
                        (rank       (length dim))
                        (total-dim  (foldl fx* 1 dim))
                        (re-length  (fvector-length ,re-vec))
                        (cpl-length (fvector-length ,cpl-vec)))
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
                 (,(if forward? 'c-fftw-execute-r2c 'c-fftw-execute-c2r)
                   (,(symbol-append 'plan- name) in out dim #f)
                   in out))))))))

  (wrap-real-trasform rfft  #t)
  (wrap-real-trasform irfft #f)
)
