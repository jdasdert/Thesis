module sortLib

  interface bubbleSort
    module procedure ibubbleSortVec, bubbleSortVec, bubbleSortMatByRow
  end interface

  contains

  subroutine ibubbleSortVec(array)
    implicit none
    logical :: swapped
    integer:: array(:), dmmy
    integer :: i
    do
      swapped = .false.
      do i=1,size(array)-1
        if(array(i).gt.array(i+1)) then
          dmmy = array(i)
          array(i) = array(i+1)
          array(i+1) = dmmy
          swapped = .true.
        endif
      enddo
      if(.not.swapped) exit
    enddo
  end subroutine

  subroutine bubbleSortVec(array)
    implicit none
    logical :: swapped
    real*8 :: array(:), dmmy
    integer :: i
    do
      swapped = .false.
      do i=1,size(array)-1
        if(array(i).gt.array(i+1)) then
          dmmy = array(i)
          array(i) = array(i+1)
          array(i+1) = dmmy
          swapped = .true.
        endif
      enddo
      if(.not.swapped) exit
    enddo
  end subroutine

  subroutine bubbleSortMatByRow(Mat,row)
    implicit none
    real*8 :: Mat(:,:), dmmy(size(Mat,1))
    integer :: row,i
    logical :: swapped
    do
      swapped = .false.
      do i=1,size(Mat,2)-1
        if(Mat(row,i).gt.Mat(row,i+1)) then
          dmmy = Mat(:,i)
          Mat(:,i) = Mat(:,i+1)
          Mat(:,i+1) = dmmy
          swapped = .true.
        endif
      enddo
      if(.not.swapped) exit
    enddo
  end subroutine
end module