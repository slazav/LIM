!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Larkin-Imry-Ma effect
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      program lim
!        call lim_init(12321, 3, 0)
!        call lim_calc(0.01D0, 1)
!        call lim_save(0, 0)
!      end

      program lim 
        implicit none
        integer i
        real*8 grad

        call lim_init(12321, 3, 0)
        grad=0.01D0
        do i=1,60
         !  call lim_init(12321, 0, 0)
          call lim_calc(grad, 1)
          call lim_save(i, 0)
          if (i.lt.30) grad=grad*dsqrt(2D0)
          if (i.ge.30) grad=grad/dsqrt(2D0)
        enddo

        do i=61,120
         !  call lim_init(12321, 0, 0)
          call lim_calc(grad, 1)
          call lim_save(i, 0)
          if (i.lt.90) grad=grad*dsqrt(2D0)
          if (i.ge.90) grad=grad/dsqrt(2D0)
        enddo

        do i=121,180
         !  call lim_init(12321, 0, 0)
          call lim_calc(grad, 1)
          call lim_save(i, 0)
          if (i.lt.150) grad=grad*dsqrt(2D0)
          if (i.ge.150) grad=grad/dsqrt(2D0)
        enddo

      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize fields A nd L (phases, 0..2pi)
! ltype -- initial value for the L field:
!   0 - constant
!   1 - same as A
!   2 - x gradient
!   3 - vortex
! atype -- initial value for the A field:
!   0 - random
!   1 - random with a hole
      subroutine lim_init(seed, ltype, atype)
        implicit none
        include 'lim.fh'
        integer seed, ltype, atype
        real*8 x,y

        call srand(seed)
        do ix=1,nx
          do iy=1,ny
            x = 2D0*dble(ix)/dble(nx)-1D0
            y = 2D0*dble(iy)/dble(ny)-1D0
            if (atype.eq.0) A(ix,iy) = dpi*rand(0)
            if (atype.eq.1) then
              A(ix,iy) = dpi*rand(0)
              if (x**2 + y**2 < 0.75D0) A(ix,iy) = 0D0
            endif

            if (ltype.eq.0) L(ix,iy) = 0D0
            if (ltype.eq.1) L(ix,iy) = A(ix,iy)
            if (ltype.eq.2) L(ix,iy) = dpi*x
            if (ltype.eq.3) L(ix,iy) = datan2(y+0.25D0,x-0.30D0)
          enddo
        enddo
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine lim_calc(grad, msg_lev)
        implicit none
        include 'lim.fh'
        real*8 grad

        integer error, lw, msg_lev
        parameter (lw = 14*nx*ny)
        real*8 w(lw)
        external calc_en

        Ug = grad
        call tn(error, nx*ny, L, E, dE, w,lw, calc_en, msg_lev)
      end

      subroutine calc_en(n,LL,E,dE)
        implicit none
        include 'lim.fh'
        integer n
        real*8 LL(nx,ny), lx,ly, gmod
        E=0D0
        do ix=1,nx
          do iy=1,ny

            ! interaction energy
            E = E + dcos(A(ix,iy)-LL(ix,iy))**2

            ! its derivative d/dLL(ix,iy)
            dE(ix,iy) = 2D0*dcos(A(ix,iy)-LL(ix,iy))*
     .                      dsin(A(ix,iy)-LL(ix,iy))

            ! gradient energy
            if (ix.lt.nx.and.iy.lt.ny) then 
              lx = gmod(LL(ix+1,iy)-LL(ix,iy), dpi)
              ly = gmod(LL(ix,iy+1)-LL(ix,iy), dpi)
              E = E + Ug*(lx**2 + ly**2)
            endif

            ! its derivative d/dLL(ix,iy)
            if (ix.lt.nx) dE(ix,iy) = dE(ix,iy)
     .           - 2D0*Ug*gmod(LL(ix+1,iy)-LL(ix,iy), dpi)
            if (ix.gt.1)  dE(ix,iy) = dE(ix,iy)
     .           - 2D0*Ug*gmod(LL(ix-1,iy)-LL(ix,iy), dpi)
            if (iy.lt.ny) dE(ix,iy) = dE(ix,iy)
     .           - 2D0*Ug*gmod(LL(ix,iy+1)-LL(ix,iy), dpi)
            if (iy.gt.1)  dE(ix,iy) = dE(ix,iy)
     .           - 2D0*Ug*gmod(LL(ix,iy-1)-LL(ix,iy), dpi)

          enddo
        enddo
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! type=0 -- phase 0..2pi
! type=1 -- phase 0..pi
      subroutine lim_save(frame, type)
        implicit none
        include 'lim.fh'
        character fname*1024
        integer type, frame

        write (fname, "(A1,I0.4,A4)") "f", frame, ".pnm"
        open(unit=5,file=trim(fname))

        ! convert to colors and print L
        call print_head(nx,ny)
        do ix=1,nx
          do iy=1,ny
            call print_col(L(ix,iy), type)
          enddo
        enddo
        close(5)
      end

      subroutine print_head(x,y)
        implicit none
        integer x,y
        write(5,'(A)') 'P6'
        write(5,'(2I5.1)') x,y
        write(5,'(A)') '255'
      end

      subroutine print_col(ph, type)
        implicit none
        real*8 ph, c, x, r,g,b, mmod
        integer  i, type
        include 'lim.fh'

        if (type.eq.0) c = 6D0*(ph/dpi - floor(ph/dpi))
        if (type.eq.1) c = 6D0*(2D0*ph/dpi - floor(2D0*ph/dpi))

        x = dmod(c, 1D0)
        i = int(c)
        if (i.eq.0) then
          r=1D0
          g=x
          b=0D0
        elseif (i.eq.1) then
          r=1D0-x
          g=1D0
          b=0D0
        elseif (i.eq.2) then
          r=0D0
          g=1D0
          b=x
        elseif (i.eq.3) then
          r=0D0
          g=1D0-x
          b=1D0
        elseif (i.eq.4) then
          r=x
          g=0D0
          b=1D0
        elseif (i.eq.5) then
          r=1D0
          g=0D0
          b=1D0-x
        endif
        write(5,'(3A1)', advance='no') int(255*r),int(255*g),int(255*b)
      end

      ! mormalize phase gradient to be -v..v
      function gmod(x, v)
        implicit none
        real*8 gmod, x, v
        gmod = (x - v*floor(x/v+0.5D0))
      end



