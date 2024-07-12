program test

  implicit none
  real*8 :: a(10)
  integer :: b(10)

  ! call random_number(a)
  ! call quicksort(a,1,10)
 !  print*, a
  call random_number(a)
  b = floor(a*10)
  call quicksortInt(b,1,10)
  print*, b


  contains

  recursive subroutine quicksort(a, first, last)
    implicit none
    real*8 :: a(*), x, t
    integer :: first, last
    integer :: i, j

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
      do while (a(i) < x)
          i=i+1
      end do
      do while (x < a(j))
          j=j-1
      end do
      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      i=i+1
      j=j-1
    end do
    if (first < i-1) call quicksort(a, first, i-1)
    if (j+1 < last)  call quicksort(a, j+1, last)
  end subroutine quicksort

  recursive subroutine quicksortInt(a, first, last)
    implicit none
    integer :: a(*), x, t
    integer :: first, last
    integer :: i, j

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
      do while (a(i) < x)
          i=i+1
      end do
      do while (x < a(j))
          j=j-1
      end do
      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      i=i+1
      j=j-1
    end do
    if (first < i-1) call quicksortInt(a, first, i-1)
    if (j+1 < last)  call quicksortInt(a, j+1, last)
  end subroutine quicksortInt
  


end program

