      PROGRAM VEC_ADD_DO
      IMPLICIT NONE
      INTEGER N, I, J, brian
      PARAMETER (N=10000)
      PARAMETER (brian = 200)
      REAL A(brian,N), B(brian)

!$omp  parallel do
!$omp& default(shared)
!$omp& private(I,J)

      DO I=1,brian
          DO J=1,N
              A(I,J) = I*J*1.0
          ENDDO
      ENDDO

!$omp  end parallel do

!$omp  parallel do
!$omp& default(shared)
!$omp& private(I,J)

      DO I=1,brian

          DO J=1,N
          
              B(I) = B(I) + A(I,J)
          
          ENDDO

      ENDDO

!$omp  end parallel do

      END
