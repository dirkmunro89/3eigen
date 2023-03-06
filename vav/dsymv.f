*> \brief \b DSYMV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER INCX,INCY,LDA,N
*       CHARACTER UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYMV  performs the matrix-vector  operation
*>
*>    y := alpha*A*x + beta*y,
*>
*> where alpha and beta are scalars, x and y are n element vectors and
*> A is an n by n symmetric matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the upper or lower
*>           triangular part of the array A is to be referenced as
*>           follows:
*>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of A
*>                                  is to be referenced.
*>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of A
*>                                  is to be referenced.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension ( LDA, N )
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n
*>           upper triangular part of the array A must contain the upper
*>           triangular part of the symmetric matrix and the strictly
*>           lower triangular part of A is not referenced.
*>           Before entry with UPLO = 'L' or 'l', the leading n by n
*>           lower triangular part of the array A must contain the lower
*>           triangular part of the symmetric matrix and the strictly
*>           upper triangular part of A is not referenced.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, n ).
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the n
*>           element vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry, BETA specifies the scalar beta. When BETA is
*>           supplied as zero then Y need not be set on input.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is DOUBLE PRECISION array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ).
*>           Before entry, the incremented array Y must contain the n
*>           element vector y. On exit, Y is overwritten by the updated
*>           vector y.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup double_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE dsymv(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
          info = 1
      ELSE IF (n.LT.0) THEN
          info = 2
      ELSE IF (lda.LT.max(1,n)) THEN
          info = 5
      ELSE IF (incx.EQ.0) THEN
          info = 7
      ELSE IF (incy.EQ.0) THEN
          info = 10
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DSYMV ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((n.EQ.0) .OR. ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
*
*     Set up the start points in  X  and  Y.
*
      IF (incx.GT.0) THEN
          kx = 1
      ELSE
          kx = 1 - (n-1)*incx
      END IF
      IF (incy.GT.0) THEN
          ky = 1
      ELSE
          ky = 1 - (n-1)*incy
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
*     First form  y := beta*y.
*
      IF (beta.NE.one) THEN
          IF (incy.EQ.1) THEN
              IF (beta.EQ.zero) THEN
                  DO 10 i = 1,n
                      y(i) = zero
   10             CONTINUE
              ELSE
                  DO 20 i = 1,n
                      y(i) = beta*y(i)
   20             CONTINUE
              END IF
          ELSE
              iy = ky
              IF (beta.EQ.zero) THEN
                  DO 30 i = 1,n
                      y(iy) = zero
                      iy = iy + incy
   30             CONTINUE
              ELSE
                  DO 40 i = 1,n
                      y(iy) = beta*y(iy)
                      iy = iy + incy
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (alpha.EQ.zero) RETURN
      IF (lsame(uplo,'U')) THEN
*
*        Form  y  when A is stored in upper triangle.
*
          IF ((incx.EQ.1) .AND. (incy.EQ.1)) THEN
              DO 60 j = 1,n
                  temp1 = alpha*x(j)
                  temp2 = zero
                  DO 50 i = 1,j - 1
                      y(i) = y(i) + temp1*a(i,j)
                      temp2 = temp2 + a(i,j)*x(i)
   50             CONTINUE
                  y(j) = y(j) + temp1*a(j,j) + alpha*temp2
   60         CONTINUE
          ELSE
              jx = kx
              jy = ky
              DO 80 j = 1,n
                  temp1 = alpha*x(jx)
                  temp2 = zero
                  ix = kx
                  iy = ky
                  DO 70 i = 1,j - 1
                      y(iy) = y(iy) + temp1*a(i,j)
                      temp2 = temp2 + a(i,j)*x(ix)
                      ix = ix + incx
                      iy = iy + incy
   70             CONTINUE
                  y(jy) = y(jy) + temp1*a(j,j) + alpha*temp2
                  jx = jx + incx
                  jy = jy + incy
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y  when A is stored in lower triangle.
*
          IF ((incx.EQ.1) .AND. (incy.EQ.1)) THEN
              DO 100 j = 1,n
                  temp1 = alpha*x(j)
                  temp2 = zero
                  y(j) = y(j) + temp1*a(j,j)
                  DO 90 i = j + 1,n
                      y(i) = y(i) + temp1*a(i,j)
                      temp2 = temp2 + a(i,j)*x(i)
   90             CONTINUE
                  y(j) = y(j) + alpha*temp2
  100         CONTINUE
          ELSE
              jx = kx
              jy = ky
              DO 120 j = 1,n
                  temp1 = alpha*x(jx)
                  temp2 = zero
                  y(jy) = y(jy) + temp1*a(j,j)
                  ix = jx
                  iy = jy
                  DO 110 i = j + 1,n
                      ix = ix + incx
                      iy = iy + incy
                      y(iy) = y(iy) + temp1*a(i,j)
                      temp2 = temp2 + a(i,j)*x(ix)
  110             CONTINUE
                  y(jy) = y(jy) + alpha*temp2
                  jx = jx + incx
                  jy = jy + incy
  120         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSYMV
*
      END
