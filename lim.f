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
        integer i, rnd_seed, ltype, atype, msg_lev
        real*8 Ugr, Uin

        rnd_seed = 12321
        ltype = 4
        atype = 0
        msg_lev = -1 ! tn message level
        Uin=1D0

        call lim_init(rnd_seed, ltype, atype)
        Ugr=0.01D0
        do i=1,60
          call lim_calc(Ugr, Uin, msg_lev)
          call lim_save(i)
          if (i.lt.30) Ugr=Ugr*sqrt(2D0)
          if (i.ge.30) Ugr=Ugr/sqrt(2D0)
        enddo

        Ugr=0D0
        call lim_calc(Ugr, Uin, msg_lev)
        Ugr=0.01D0

        do i=61,90
          call lim_calc(Ugr, Uin, msg_lev)
          call lim_save(i)
          Ugr=Ugr*sqrt(2D0)
        enddo


      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize fields A and L (phases, 0..2pi)
! ltype -- initial value for the L field:
!   0 - constant
!   1 - same as A
!   2 - x gradient
!   3 - vortex
!   4 - ++ and +- vortex pairs
! atype -- initial value for the A field:
!   0 - random phase
!   1 - random
!   2 - random with a hole
      subroutine lim_init(seed, ltype, atype)
        implicit none
        include 'lim.fh'
        integer seed, ltype, atype
        real*8 x,y,v

        call srand(seed)
        do ix=1,nx
          do iy=1,ny
            x = 2D0*dble(ix)/dble(nx)-1D0
            y = 2D0*dble(iy)/dble(ny)-1D0

            ! A field
            if (atype.eq.0) then ! unit length
              v = dpi*rand(0)
              Ax(ix,iy) = dcos(v)
              Ay(ix,iy) = dsin(v)
            else if (atype.eq.1) then
              Ax(ix,iy) = 2D0*rand(0)-1D0
              Ay(ix,iy) = 2D0*rand(0)-1D0
            else if (atype.eq.2) then
              Ax(ix,iy) = 2D0*rand(0)-1D0
              Ay(ix,iy) = 2D0*rand(0)-1D0
              if (x**2 + y**2 < 0.75D0) then
                Ax(ix,iy) = 0D0
                Ay(ix,iy) = 0D0
              endif
            endif

            ! L field
            if (ltype.eq.0) L(ix,iy) = 0D0
            if (ltype.eq.1) L(ix,iy) = datan2(Ay(ix,iy),Ax(ix,iy))
            if (ltype.eq.2) L(ix,iy) = 2D0*dpi*x
            if (ltype.eq.3) L(ix,iy) = datan2(y+0.25D0,x-0.30D0)
            if (ltype.eq.4) L(ix,iy) = datan2(y+0.25D0,x-0.30D0)
     .                               + datan2(y+0.45D0,x-0.20D0)
     .                               + datan2(y+0.45D0,x-0.20D0)
     .                               - datan2(y-0.65D0,x+0.62D0)
     .                               + datan2(y-0.55D0,x+0.70D0)

          enddo
        enddo
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine lim_calc(Ugr, Uin, msg_lev)
        implicit none
        include 'lim.fh'
        real*8 Ugr, Uin

        integer error, lw, msg_lev
        parameter (lw = 14*nx*ny)
        real*8 w(lw)
        external calc_en

        Ug = Ugr
        Ui = Uin
        call tn(error, nx*ny, L, E, dE, w,lw, calc_en, msg_lev)
      end

      subroutine calc_en(n,LL,E,dE)
        implicit none
        include 'lim.fh'
        integer n
        real*8 LL(nx,ny), lx,ly, v, gmod
        E=0D0
        do ix=1,nx
          do iy=1,ny

            ! interaction energy and its derivative d/dLL(ix,iy)
            v = datan2(Ay(ix,iy),Ax(ix,iy))
            if (1.eq.1) then
              E = E - Ui*dcos(v-LL(ix,iy))**2
              dE(ix,iy) = -2D0*Ui*dcos(v-LL(ix,iy))*
     .                            dsin(v-LL(ix,iy))
            else
              E = E - Ui*dcos(v-LL(ix,iy))
              dE(ix,iy) = -2D0*Ui*dsin(v-LL(ix,iy))
            endif

            ! gradient energy
            if (ix.lt.nx.and.iy.lt.ny) then 
              lx = gmod(LL(ix+1,iy)-LL(ix,iy))
              ly = gmod(LL(ix,iy+1)-LL(ix,iy))
              E = E + Ug*(lx**2 + ly**2)
            endif

            ! its derivative d/dLL(ix,iy)
            if (ix.lt.nx) dE(ix,iy) = dE(ix,iy)
     .           - 2D0*Ug*gmod(LL(ix+1,iy)-LL(ix,iy))
            if (ix.gt.1)  dE(ix,iy) = dE(ix,iy)
     .           - 2D0*Ug*gmod(LL(ix-1,iy)-LL(ix,iy))
            if (iy.lt.ny) dE(ix,iy) = dE(ix,iy)
     .           - 2D0*Ug*gmod(LL(ix,iy+1)-LL(ix,iy))
            if (iy.gt.1)  dE(ix,iy) = dE(ix,iy)
     .           - 2D0*Ug*gmod(LL(ix,iy-1)-LL(ix,iy))

          enddo
        enddo
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine lim_save(frame)
        implicit none
        include 'lim.fh'
        character fname*1024
        integer frame

        write (fname, "(A1,I0.4,A4)") "f", frame, ".pnm"
        open(unit=5,file=trim(fname))

        ! convert to colors and print L
        call print_head(nx,ny)
        do iy=1,ny
          do ix=1,nx
            call print_col(L(ix,iy))
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

      subroutine print_col(ph)
        implicit none
        real*8 ph, c, x, r,g,b, mmod
        integer  i
        include 'lim.fh'
        c = 6D0*(ph/dpi - floor(ph/dpi))

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

      ! normalize phase difference to be a -pi..pi saw
      function gmod(x)
        implicit none
        real*8 gmod, x, dpi
        parameter (dpi=6.2831853071795864D0)
        gmod = (x - dpi*floor(x/dpi+0.5D0))
      end



