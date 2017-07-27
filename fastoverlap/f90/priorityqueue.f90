module priorityqueue
implicit none

type node
  double precision :: vector(3), upperbound, lowerbound, width
  integer :: niter, idnum
end type

type queue
  type(node), allocatable :: buf(:)
  integer                 :: n = 0
contains
  procedure :: top
  procedure :: enqueue
  procedure :: siftdown
  procedure :: destroy
end type

contains

subroutine siftdown(this, a)
  class (queue)           :: this
  integer                 :: a, parent, child
  associate (x => this%buf)
  parent = a
  do while(parent*2 <= this%n)
    child = parent*2
    if (child + 1 <= this%n) then
      if (x(child+1)%lowerbound < x(child)%lowerbound ) then
        child = child +1
      end if
    end if
    if (x(parent)%lowerbound > x(child)%lowerbound) then
      x([child, parent]) = x([parent, child])
      parent = child
    else
      exit
    end if
  end do
  end associate
end subroutine

function top(this) result (res)
  class(queue) :: this
  type(node)   :: res
  res = this%buf(1)
  this%buf(1) = this%buf(this%n)
  this%n = this%n - 1
  call this%siftdown(1)
end function

subroutine enqueue(this, lowerbound, upperbound, vector, width, niter, idnum)
  class(queue), intent(inout) :: this
  double precision            :: lowerbound, upperbound, vector(3), width
  integer                     :: niter, idnum
  type(node)                  :: x
  type(node), allocatable     :: tmp(:)
  integer                     :: i
  x%vector = vector
  x%upperbound = upperbound
  x%lowerbound = lowerbound
  x%width = width
  x%niter = niter
  x%idnum = idnum
  this%n = this%n +1
  if (.not.allocated(this%buf)) allocate(this%buf(1))
  if (size(this%buf)<this%n) then
    allocate(tmp(2*size(this%buf)))
    tmp(1:this%n-1) = this%buf
    call move_alloc(tmp, this%buf)
  end if
  this%buf(this%n) = x
  i = this%n
  do
    i = i / 2
    if (i==0) exit
    call this%siftdown(i)
  end do
end subroutine

subroutine destroy(this)
   class(queue), intent(inout) :: this

   deallocate(this%buf)
end subroutine

!subroutine putinqueue(q, lowerbound, upperbound, vector, width, niter, idnum)
!
!    implicit none
!    class(queue), intent(inout) :: q
!    double precision            :: lowerbound, upperbound, vector(3), width
!    integer                     :: niter, idnum
!
!    call enqueue(q, lowerbound, upperbound, vector, width, niter, idnum)
!
!end subroutine

end module
