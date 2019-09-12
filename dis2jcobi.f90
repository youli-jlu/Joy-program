!change distance matrix of H2O-Ar to  Jacobi coordinates
!First, change distance matrix to cartesian.
!Then, translation to H2O mass centre
!Last, calculate jacobi. Here I use rOH instead of Q2,Q3.


        program main
          implicit none
          real*8 dm(6),jacobi(6)
          real*8 dm_ro(6)  !reorder's distance matrix
          integer :: filestate= 0

          open(780,file="rr6.out",action="read")
          open(880,file="result",action="write")

          write(880,*) "         R1(A)         R2(A)        TH(degree)     phi      theta      R(A)    "

          do while(.true.)
            dm=0.0d0
            read(780,*,iostat=filestate) dm
            if ( filestate/=0 ) stop
            dm_ro(1)=dm(2)
            dm_ro(2)=dm(4)
            dm_ro(3)=dm(6)
            dm_ro(4)=dm(1)
            dm_ro(5)=dm(3)
            dm_ro(6)=dm(5)
            call dis2jacobi(dm_ro,jacobi)
            write(880,*) jacobi
          enddo
          close(780)
          close(880)

          end program main


          !change distance matrix to H2O-Ar jacobi coordinates
        Subroutine dis2jacobi(dm,jacobi)
          implicit none
          real*8 dm(6)
          real*8 d(4,4),coord(4,3)
          real*8 r1,r2,theta1,r_ar,theta2,phi
          real*8 jacobi(6)
          integer i,j,k
          real*8,parameter :: pi = 3.14159265359
          real*8 mc(3)
          real*8,external :: dihedral

          !read distance matrix
          !in this program, label [1,2,3,4] refer to [O,H1,H2,Ar]
          !open(778,file="tetra.dat",action="read")
          d=0.0d0
          k=1
          do i=1,3
            do j=i+1,4
              d(i,j)=dm(k)  ! I trust your input and set no input filter 
              d(j,i)=d(i,j)
              k=k+1
            enddo
          enddo
          !close(778)

          !change distance matrix to cartesian coordintes
          !coord(i,j) refers to j axies coordinates of [i]'s atoms
          !specifically,I remain the dihedral phi.
          call dis2cartesian(d,coord,theta1)

          !change cartesian to jacobi coordinates
          !!!!!r1,r2 means OH distances. If you wanto use Q2,Q3, you need to change this
          r1=d(1,2)
          r2=d(1,3)
          

          !calculate H2O 's mass centre and change our coordinates
          open(890,file="checkcoord",action="write",position="append")

          write(890,*) "old coord  "
          do i = 1,4
          write(890,*) coord(i,:)
          enddo

          call mass_centre(coord)
          mc=0.0d0

          write(890,*) "new coord  "
          do i =1,4
          write(890,*) coord(i,:)
          enddo
          write(890,*) "  "

          r_ar = dsqrt(coord(4,1)*coord(4,1)+coord(4,2)*coord(4,2)+coord(4,3)*coord(4,3))
          theta2=( (coord(1,1)*coord(4,1)+coord(1,2)*coord(4,2)+coord(1,3)*coord(4,3) ) &
            & / ( dsqrt(coord(4,1)**2+coord(4,2)**2+coord(4,3)**2) * dsqrt(coord(1,1)**2+coord(1,2)**2+coord(1,3)**2) ) )

          theta2 = nint( theta2 * 1.0d5) * 1d-5 !avoid numeriacal error
          theta2=dacos(theta2)

          phi = dihedral(coord(2,:),coord(1,:),mc,coord(4,:))

          jacobi(1)=r1
          jacobi(2)=r2
          jacobi(3)=theta1*180.0d0/pi
          jacobi(4)=phi*180.0d0/pi
          jacobi(5)=theta2*180.0d0/pi
          jacobi(6)=r_ar

        end subroutine dis2jacobi


        subroutine mass_centre(coord)
          implicit none
          integer,parameter :: atom=4
          real*8 coord(4,3)
          integer i,j,k
          real*8 m(atom),x(atom),y(atom),z(atom)
          real*8 m0,mx,my,mz
          real*8 smx,smy,smz
          character(len=2) ch(atom)

          x=coord(:,1)
          y=coord(:,2)
          z=coord(:,3)
          m(1)=15.9994d0
          m(2)=1.00784d0
          m(3)=1.00784d0

          m0=0.0d0
          mx=0.0d0
          my=0.0d0
          mz=0.0d0

          !!!we only calculate H2O 's mc
          do i=1,atom-1
          m0=m0+m(i)
          ! write(*,*) m(i)
          mx=mx+m(i)*x(i)
          my=my+m(i)*y(i)
          mz=mz+m(i)*z(i)
          enddo

          !!!change all the coordinates
          do i=1,atom
          x(i)=x(i)-mx/m0
          y(i)=y(i)-my/m0
          z(i)=z(i)-mz/m0
          enddo
          coord(:,1)=x
          coord(:,2)=y
          coord(:,3)=z


        end subroutine mass_centre



!!change distance matrix to cartesian coordinates
!!using Cayley–Menger determinant to calculate volume 

        subroutine dis2cartesian(d,coord,theta)
          implicit none
          real*8 d(4,4)
          real*8 phi,theta,theta314,theta214
          integer i,j,k,l,noi,noj
          real*8,parameter :: pi = 3.14159265359
          real*8 coord(4,3)  !!coordinates of O,H1,H2,Ar
          real*8 h_ar

          !!2019.9.12: I made a big mistake. I really don't need to use volume.
          !!As I can just build HOH coord easily, Ican just sovle the distance equation
          !!same way to set O,H1,H2
          coord(1,:) = 0.0d0
          coord(2,1) = d(1,2) !H1 is in (r12,0,0)
          coord(2,2:3) = 0.0d0
          theta= ( d(1,2)**2+d(1,3)**2-d(2,3)**2 ) / ( 2.0d0*d(1,2)*d(1,3) )
          theta= dacos(theta)
          coord(3,1) = d(1,3)*dcos(theta)
          coord(3,2) = d(1,3)*dsin(theta)
          coord(3,3) = 0.0d0


          !set Ar coord by distance equation
          coord(4,1) = (d(1,2)**2+d(1,4)**2-d(2,4)**2) / (2.0d0*d(1,2))
          coord(4,2) = (d(3,4)**2-d(1,4)**2+2.0d0*coord(3,1)*coord(4,1)-coord(3,1)**2-coord(3,2)**2) &
            & / (-2.0d0 * coord(3,2))
          coord(4,3) = (d(3,4)**2-(coord(4,1)-coord(3,1))**2-(coord(4,2)-coord(3,2))**2)
          if (coord(4,3).lt.-0.001d0) then
            write(*,*) "distance matrix wrong"
            stop
          else 
            if(coord(4,3).lt.1.0d-4) coord(4,3)=0.0d0 ! avoid numerical error in calculation. cutoff=10-2
          end if 

          coord(4,3) = dsqrt(coord(4,3))




        end subroutine dis2cartesian




        ! caculate dihedral or plane123 and plane234

        function dihedral(p1,p2,p3,p4)
          implicit none 
        !  real*8,dimension(3) :: p1,p2,p3,p4 !point cartesian coordinates 
        !  real*8,dimension(3) :: e12,e23,e24
        !  real*8,dimension(3) :: n1,n2,n3
          real*8 :: p1(3),p2(3),p3(3),p4(3) !point cartesian coordinates 
          real*8 :: e12(3),e23(3),e24(3)
          real*8 :: n1(3),n2(3),n3(3)
          real*8  dihedral

          e12 = p1-p2
          e23 = p2-p3
          e24 = p2-p4
          call cross(e12,e23,n1)
          call cross(e23,e24,n2)
          dihedral = ( dabs(n1(1)*n2(1)+n1(2)*n2(2)+n1(3)*n2(3)) & 
            & / ( dsqrt( n1(1)**2+n1(2)**2+n1(3)**2) * dsqrt( n2(1)**2+n2(2)**2+n2(3)**2) ) )
          dihedral = nint( dihedral * 1.0d5) * 1d-5 !avoid numeriacal error        
          return 

          end function dihedral

        subroutine cross(e1,e2,n1)
          implicit none
          real*8 :: e1(3)
          real*8 :: e2(3)
          real*8 :: n1(3)

          n1(1) = e1(2) * e2(3) - e1(3) * e2(2)
          n1(2) = e1(3) * e2(1) - e1(1) * e2(3)
          n1(3) = e1(1) * e2(2) - e1(2) * e2(1)


          end subroutine cross



