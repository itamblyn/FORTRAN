program concat
implicit none

integer intvar, i
character(4) string
character(13) fout

intvar = 1000

do i=1000, 1010
  write(unit=string, fmt='(I4)') intvar
  fout = "file" // string // ".out"
  open(2,file=fout)
  write(2,*) "test"
  close(2)
  write(*,*) fout
  intvar = intvar + 1
end do

end program concat
