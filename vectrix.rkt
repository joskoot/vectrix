#lang racket

(provide
 vectrix-ref
 vectrix-set!
 make-rectangular-vectrix
 build-rectangular-vectrix
 make-triangular-vectrix
 build-triangular-vectrix
 vectrix-map
 vectrix-for-each
 copy-vectrix
 vectrix-dimensions
 vectrix-product)

(define (vectrix-ref vectrix . indices) (xvectrix-ref vectrix (flatten indices)))
(define (vectrix-set! value vectrix . indices) (xvectrix-set! value vectrix (flatten indices)))

(define (make-rectangular-vectrix  fill . dimensions)
 (xmake-rectangular-vectrix  fill (flatten dimensions)))
(define (build-rectangular-vectrix proc . dimensions)
  (xbuild-rectangular-vectrix proc (flatten dimensions)))

(define (flatten x)
 (cond
  ((list? x) (apply append (map flatten x)))
  ((pair? x) (append (flatten (car x)) (flatten (cdr x))))
  ((box?  x) (flatten (unbox x)))
  ((vector? x) (flatten (vector->list x)))
  (else (list x))))

(define (xvectrix-ref vectrix indices)
 (if (null? indices) vectrix
  (xvectrix-ref (vector-ref vectrix (car indices)) (cdr indices))))

(define (xvectrix-set! value vectrix indices)
 (let ((index (car indices)) (indices (cdr indices)))
  (if (null? indices) (vector-set! vectrix index value)
   (xvectrix-set! value (vector-ref vectrix index) indices))))

(define (xmake-rectangular-vectrix fill dimensions)
 (if (null? dimensions) fill
  (let ((dimension (car dimensions)) (dimensions (cdr dimensions)))
   (build-vector dimension
    (lambda (index) (xmake-rectangular-vectrix fill dimensions))))))

(define (xbuild-rectangular-vectrix proc dimensions)
 (if (null? dimensions) (proc)
  (let ((dimension (car dimensions)) (dimensions (cdr dimensions)))
   (build-vector dimension
    (lambda (index)
     (xbuild-rectangular-vectrix (lambda indices (apply proc index indices))
      dimensions))))))

(define (make-triangular-vectrix fill dim depth)
 (if (zero? depth) fill
  (build-vector dim
   (lambda (index)
    (make-triangular-vectrix fill
     (add1 index) (sub1 depth))))))

(define (build-triangular-vectrix proc dim depth)
 (if (zero? depth) (proc)
  (build-vector dim
   (lambda (index)
    (build-triangular-vectrix (lambda indices (apply proc index indices))
     (add1 index) (sub1 depth))))))

(define (vectrix-map proc depth . vectrices)
 (if (zero? depth) (apply proc vectrices)
  (apply vector-map
   (lambda (index . vectrices)
    (apply vectrix-map
     (lambda args (apply proc index args)) (sub1 depth) vectrices))
   vectrices)))

(define (vectrix-for-each proc depth . vectrices)
 (if (zero? depth) (void (apply proc vectrices))
  (apply vector-for-each
   (lambda (index . vectrices)
    (apply vectrix-for-each
     (lambda args (apply proc index args)) (sub1 depth) vectrices))
   vectrices)))

(define (copy-vectrix vectrix depth)
 (if (zero? depth) vectrix
  (vector-map
   (lambda (index subvectrix) (copy-vectrix subvectrix (sub1 depth)))
   vectrix)))

(define vectrix-dimensions
 (case-lambda
  ((vectrix)
   (let vectrix-dimensions ((vectrices (list vectrix)))
    (if (not (andmap vector? vectrices)) '()
     (let ((kar (car vectrices)) (kdr (cdr vectrices)))
      (let ((len (vector-length kar)))
       (if (zero? len) '()
        (if (not (andmap (lambda (y) (= (vector-length y) len)) kdr)) '()
         (cons len (vectrix-dimensions (apply append (map vector->list vectrices)))))))))))
  ((vectrix depth)
   (let vectrix-dimensions ((vectrix vectrix) (depth depth))
    (if (zero? depth) '()
     (cons (vector-length vectrix) (vectrix-dimensions (vector-ref vectrix 0) (sub1 depth))))))))
 
(define-syntax (vectrix-product stx)
 (define (get-sum-indices result-indices all-indices)
  (define (id-member? id id-list)
   (and (not (null? id-list))
    (or (bound-identifier=? id (car id-list))
     (id-member? id (cdr id-list)))))
  (let loop ((all all-indices) (done result-indices))
   (if (null? all) '()
    (let ((a (car all)) (all (cdr all)))
     (if (id-member? a done) (loop all done)
      (cons a (loop all (cons a done))))))))
 (define (sum-index-control sum-index) #`(#,sum-index (get-dimension '#,sum-index)))
 (syntax-case stx ()
  ((vectrix-product adder multiplier (result-index ...) (vectrix index ...) ...)
   (let ((sum-indices
          (get-sum-indices (syntax->list #'(result-index ...)) (syntax->list #'(index ... ...)))))
  #`(let ((dimension-table (make-hash)) (all-indices '((index ...) ...)))
     (define (get-dimension id) (hash-ref dimension-table id))
     (define (put-dimension id dim) (hash-set! dimension-table id dim))
     (define (put-dimensions ids dims) (for-each put-dimension ids dims))
     (define add adder)
     (define mul multiplier)
     (define vectrices (list vectrix ...))
     (for-each put-dimensions '((index ...) ...)
      (list (vectrix-dimensions vectrix (length '(index ...))) ...))
     (apply build-rectangular-vectrix
      (lambda (result-index ...)
       (let ((products '()))
        (for #,(map sum-index-control sum-indices)
         (set! products (cons (mul (vectrix-ref vectrix index ...) ...) products)))
        (apply add (reverse products))))
      (map get-dimension '(result-index ...))))))))

(define (vector-for-each proc . vectors)
 (let ((n (vector-length (car vectors))))
  (for ((index (vector-length (car vectors))))
   (apply proc index (map (lambda (vector) (vector-ref vector index)) vectors)))))

(define (vector-map proc . vectors)
 (build-vector (vector-length (car vectors))
  (lambda (index)
   (apply proc index (map (lambda (vector) (vector-ref vector index)) vectors)))))

