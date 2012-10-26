module vectors

  implicit none
contains

  function length(v)
    real(8)    :: length, v(:), s
    integer :: n, i
    n = size(v)
    s = 0.0
    do i=1,n
      s = s + v(i)**2
    end do
    length=sqrt(s)
  end function

  function cross(a,b)

    real(8)    :: a(3), b(3), cross(3) 
    
    cross(1) =  (a(2)*b(3)-a(3)*b(2))
    cross(2) = -(a(1)*b(3)-a(3)*b(1))
    cross(3) =  (a(1)*b(2)-a(2)*b(1))

  end function

end module
