
real(4) a(3,3), b(3), c(3),ans
integer i,j

a(:,:) = 0

print *, Pi

do i=1,3
  do j=1,3

    a(i,j) =  i
    b(i) = -i

  end do
end do

!a = a/10.0
!print *, a
!print *, b

!c = a - b

!print *, dot_product(a,b)
!print *, c

!print *, a

!a(-1) = 2 

!print *, a 

print *, minval(a)

ans = erf(1.0)

print *, ans

END
