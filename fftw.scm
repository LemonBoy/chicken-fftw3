#>
#include <fftw3.h>
<#

(define-foreign-variable c-fftw-estimate int "FFTW_ESTIMATE")
(define-foreign-variable c-fftw-measure int "FFTW_MEASURE")
(define-foreign-variable c-fftw-patient int "FFTW_PATIENT")
(define-foreign-variable c-fftw-exhaustive int "FFTW_EXHAUSTIVE")
(define-foreign-variable c-fftw-destroy-input int "FFTW_DESTROY_INPUT")
(define-foreign-variable c-fftw-preserve-input int "FFTW_PRESERVE_INPUT")
(define-foreign-variable c-fftw-unaligned int "FFTW_UNALIGNED")

(functor
  (fftw-impl
    (I (c-fftw-destroy-plan c-fftw-alignment-of c-fftw-execute-c2c
         c-fftw-execute-c2r c-fftw-execute-r2c c-fftw-execute-r2r
         c-fftw-plan-c2c c-fftw-plan-r2c c-fftw-plan-c2r c-fftw-plan-r2r)))
  (plan-fft plan-rfft plan-ifft plan-irfft plan-dct plan-dst
    fftw-estimate fftw-measure fftw-patient fftw-exhaustive)

  (import chicken scheme foreign I)
  (use srfi-4)

  (define-foreign-variable FFTW_REDFT00 int "FFTW_REDFT00")
  (define-foreign-variable FFTW_REDFT10 int "FFTW_REDFT10")
  (define-foreign-variable FFTW_REDFT01 int "FFTW_REDFT01")
  (define-foreign-variable FFTW_REDFT11 int "FFTW_REDFT11")
  (define-foreign-variable FFTW_RODFT00 int "FFTW_RODFT00")
  (define-foreign-variable FFTW_RODFT10 int "FFTW_RODFT10")
  (define-foreign-variable FFTW_RODFT01 int "FFTW_RODFT01")
  (define-foreign-variable FFTW_RODFT11 int "FFTW_RODFT11")

  (define +even-flags+
    (list FFTW_REDFT00 FFTW_REDFT10 FFTW_REDFT01 FFTW_REDFT11))
  (define +odd-flags+
    (list FFTW_RODFT00 FFTW_RODFT10 FFTW_RODFT01 FFTW_RODFT11))

  (define (kind->enum even? x loc)
    (##sys#check-range x 1 4 loc)
    (list-ref (if even? +even-flags+ +odd-flags+) (fx- x 1)))

  (define (vector-alignment v)
    (c-fftw-alignment-of (##sys#slot v 1)))

  (define (flags->bitmap args)
    (let ((strategy (or (get-keyword strategy: args) c-fftw-estimate)))
      (##sys#check-exact strategy strategy:)
      (apply bitwise-ior
        (cons strategy
          (map (lambda (k+v)
                 (if (get-keyword (car k+v) args) (cdr k+v) 0))
               `((unaligned?:      . ,c-fftw-unaligned)
                 (preserve-input?: . ,c-fftw-preserve-input)
                 (destroy-input?:  . ,c-fftw-destroy-input)))))))

  ; exported flags

  (define fftw-estimate c-fftw-estimate)
  (define fftw-measure c-fftw-measure)
  (define fftw-patient c-fftw-patient)
  (define fftw-exhaustive c-fftw-exhaustive)

  ; exported procedures

  (define-syntax wrap-real-eo-transform
    (er-macro-transformer
      (lambda (x r c)
        (let* ((name (strip-syntax (cadr x)))
               (even? (strip-syntax (caddr x)))
               (proc-name (symbol-append 'plan- name)))
          `(define (,proc-name kind in out #!optional dim . flags)
             (let ((dim (or dim (list (fvector-length in))))
                   (flags (flags->bitmap flags))
                   (total-dim (foldl fx* 1 dim)))
               ; make sure the parameters are correct
               (cond
                   ((< total-dim 1)
                    (error "invalid transform size" dim))
                   ((> total-dim (fvector-length out))
                    (error "output vector length is too short"))
                   ((> total-dim (fvector-length in))
                    (error "input vector length is too short")))
               (let ((ensure-aligned? (eq? 0 (fxand c-fftw-unaligned flags)))
                     (in-place? (eq? in out))
                     (rank (length dim))
                     (align-in  (vector-alignment in))
                     (align-out (vector-alignment out))
                     (enum (kind->enum ,even? kind ',proc-name))
                     (plan #f))
                 ; create the plan
                 (set! plan
                   (c-fftw-plan-r2r rank (list->s32vector dim) in out
                     #$(make-s32vector rank enum) flags))
                 ; FFTW caches the plans, delay the destruction
                 (set-finalizer! plan c-fftw-destroy-plan)
                 (lambda (in out)
                   (cond
                       ((> total-dim (fvector-length out))
                        (error "output vector length is too short"))
                       ((> total-dim (fvector-length in))
                        (error "input vector length is too short"))
                       ((not (eq? in-place? (eq? in out)))
                        (error "planning mismatch"))
                       ((and ensure-aligned?
                             (eq? align-in  (vector-alignment in))
                             (eq? align-out (vector-alignment out)))
                        (error "alignment mismatch")))
                   (c-fftw-execute-r2r plan in out)))))))))

  (wrap-real-eo-transform dct #t)
  (wrap-real-eo-transform dst #f)

  (define-syntax wrap-complex-trasform
    (er-macro-transformer
      (lambda (x r c)
        (let* ((name (strip-syntax (cadr x)))
               (forward? (strip-syntax (caddr x)))
               (expt (if forward? -1 +1))
               (proc-name (symbol-append 'plan- name)))
          `(define (,proc-name kind in out #!optional dim . flags)
             (let ((dim (or dim (list (fvector-length in))))
                   (flags (flags->bitmap flags))
                   (total-dim (foldl fx* 1 dim)))
               ; make sure the parameters are correct
               (cond
                   ((< total-dim 1)
                    (error "invalid transform size" dim))
                   ((> (fx* 2 total-dim) (fvector-length out))
                    (error "output vector length is too short"))
                   ((> (fx* 2 total-dim) (fvector-length in))
                    (error "input vector length is too short")))
               (let ((ensure-aligned? (eq? 0 (fxand c-fftw-unaligned flags)))
                     (in-place? (eq? in out))
                     (rank (length dim))
                     (align-in  (vector-alignment in))
                     (align-out (vector-alignment out))
                     (plan #f))
                 ; create the plan
                 (set! plan
                   (c-fftw-plan-c2c rank (list->s32vector dim) in out ,expt flags))
                 ; FFTW caches the plans, delay the destruction
                 (set-finalizer! plan c-fftw-destroy-plan)
                 (lambda (in out)
                   (cond
                       ((> (fx* 2 total-dim) (fvector-length out))
                        (error "output vector length is too short"))
                       ((> (fx* 2 total-dim) (fvector-length in))
                        (error "input vector length is too short"))
                       ((not (eq? in-place? (eq? in out)))
                        (error "planning mismatch"))
                       ((and ensure-aligned?
                             (eq? align-in  (vector-alignment in))
                             (eq? align-out (vector-alignment out)))
                        (error "alignment mismatch")))
                   (c-fftw-execute-c2c plan in out)))))))))

  (wrap-complex-trasform fft  #t)
  (wrap-complex-trasform ifft #f)

  (define-syntax wrap-real-trasform
    (er-macro-transformer
      (lambda (x r c)
        (let* ((name (strip-syntax (cadr x)))
               (forward? (strip-syntax (caddr x)))
               (re-vec (if forward? 'in 'out))
               (cpl-vec (if forward? 'out 'in))
               (proc-name (symbol-append 'plan- name))
               (exec-name (if forward?
                              'c-fftw-execute-r2c
                              'c-fftw-execute-c2r))
               (plan-name (if forward?
                              'c-fftw-plan-r2c
                              'c-fftw-plan-c2r)))
          `(define (,proc-name kind in out #!optional dim . flags)
             (let ((dim (or dim (list (fvector-length in))))
                   (flags (flags->bitmap flags))
                   (total-dim (foldl fx* 1 dim)))
               ; make sure the parameters are correct
               (cond
                   ((< total-dim 1)
                    (error "invalid transform size" dim))
                   ((> total-dim (fvector-length ,re-vec))
                    (error "real vector length is too short"))
                   ((> (fx+ 2 total-dim) (fvector-length ,cpl-vec))
                    (error "complex vector length is too short")))
               (let ((ensure-aligned? (eq? 0 (fxand c-fftw-unaligned flags)))
                     (in-place? (eq? in out))
                     (rank (length dim))
                     (align-in  (vector-alignment in))
                     (align-out (vector-alignment out))
                     (plan #f))
                 ; create the plan
                 (set! plan
                   (,plan-name rank (list->s32vector dim) in out flags))
                 ; FFTW caches the plans, delay the destruction
                 (set-finalizer! plan c-fftw-destroy-plan)
                 ; The Hermitian symmetry allows us to spare some space by
                 ; only storing half of the spectrum
                 ; (N/2 + 1) * 2 = N + 2
                 (lambda (in out)
                   (cond
                       ((> (fx+ 2 total-dim) (fvector-length ,cpl-vec))
                        (error "complex vector length is too short"))
                       ((> total-dim (fvector-length ,re-vec))
                        (error "real vector length is too short"))
                       ((not (eq? in-place? (eq? in out)))
                        (error "planning mismatch"))
                       ((and ensure-aligned?
                             (eq? align-in  (vector-alignment in))
                             (eq? align-out (vector-alignment out)))
                        (error "alignment mismatch")))
                   (,exec-name plan in out)))))))))

  (wrap-real-trasform rfft  #t)
  (wrap-real-trasform irfft #f)
)
