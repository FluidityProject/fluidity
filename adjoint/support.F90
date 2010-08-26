#include "fdebug.h"

module pod_support

  use ale_module
  use allsorts
  use fldebug
  use spaerr_module
  use snapsvd_module

  ! Beware! Implicit typing used in this module!
  !implicit none

  private

  public :: buildsnapmat, fulltopodproj, fulltopodproj_multivar, &
    & podtofullproj_multivar, setobsv

contains

      subroutine fulltopodproj(n,nsvd,varin,leftsvd,tcoef)
!   projects varin on POD space: tcoef = leftsvd^t*varin

      implicit none

      integer n,nsvd
      REAL varin(n), leftsvd(n,nsvd),tcoef(nsvd)
      integer i,j

      do i=1,nsvd
         tcoef(i) = 0.0d0
         do j=1,n
            tcoef(i) = tcoef(i)+leftsvd(j,i)*varin(j)
         enddo
      enddo

      return
      endsubroutine fulltopodproj

      subroutine fulltopodproj_multivar(nsvd_total,nsvd,nsvd_u,nsvd_phi,nvar,mxpoi,varin,leftsvd,tcoef)
!   projects varin on POD space: tcoef = leftsvd^t*varin

      implicit none

      integer nsvd_total,nsvd_u,nsvd_phi,nvar,mxpoi,nsvd
      REAL varin(nvar*mxpoi), tcoef(nsvd_total)
      double precision :: leftsvd(nvar*mxpoi,nsvd)
      integer i,j,k

      do k=1,nvar
         if(k.ne.nvar) then
            do i=1,nsvd_u
               tcoef((k-1)*nsvd_u+i) = 0.0d0
               do j=1,mxpoi
                  tcoef( (k-1)*nsvd_u+i ) = tcoef((k-1)*nsvd_u+i)+           &
                       leftsvd((k-1)*mxpoi+j,i)*varin((k-1)*mxpoi+j)
               enddo
            enddo
         else
            do i=1,nsvd_phi
               tcoef((k-1)*nsvd_u+i) = 0.0d0
               do j=1,mxpoi
                  tcoef( (k-1)*nsvd_u+i ) = tcoef((k-1)*nsvd_u+i)+           &
                       leftsvd((k-1)*mxpoi+j,i)*varin((k-1)*mxpoi+j)
               enddo
            enddo
         endif
      enddo


      return
      endsubroutine fulltopodproj_multivar

      subroutine podtofullproj(n,nsvd,varin,leftsvd,tcoef)
      
      implicit none
      
      integer n,nsvd
      REAL varin(n), leftsvd(n,nsvd),tcoef(nsvd)
      integer i,j

        do i=1,n
          varin(i) = 0.0d0
        enddo
        
        do i=1,n
           do j=1,nsvd
              varin(i) = varin(i)+tcoef(j)*leftsvd(i,j)
           enddo
       enddo

      return
      end subroutine podtofullproj
                        
      subroutine podtofullproj_multivar(nsvd_total,nsvd,nsvd_u,nsvd_phi,nvar,mxpoi,varin,leftsvd,tcoef)
      
      implicit none
      
      integer nsvd_total,nsvd_u,nsvd_phi,nvar,mxpoi,nsvd
      REAL varin(nvar*mxpoi), tcoef(nsvd*nvar)
      double precision :: leftsvd(nvar*mxpoi,nsvd)
      integer i,j,k

        do i=1,nvar*mxpoi
          varin(i) = 0.0d0
        enddo
        
        do k=1,nvar
           if(k.ne.nvar) then
              do i=1,nsvd_u
                 do j=1,mxpoi
                    varin((k-1)*mxpoi+j) = varin((k-1)*mxpoi+j)+tcoef((k-1)*nsvd_u+i)*leftsvd((k-1)*mxpoi+j,i)
                 enddo
              enddo
           else
              do i=1,nsvd_phi
                 do j=1,mxpoi
                    varin((k-1)*mxpoi+j) = varin((k-1)*mxpoi+j)+tcoef((k-1)*nsvd_u+i)*leftsvd((k-1)*mxpoi+j,i)
                 enddo
              enddo
           endif
        enddo

!        ewrite(3,*) varin(1),tcoef(1),leftsvd(1,1)

      return
      end subroutine podtofullproj_multivar
      
      subroutine buildsnapmat(snapmatrix,smean,iuseobs,ntime,  &
                             obsv,iflagobs,iusemean,Tnsnap)
                             
      implicit none

!********** builds snapshots matrix

      include 'paramnew.h'
      integer iuseobs,iusemean,ntime,i,j,Tnsnap
      REAL smean(istate)
      double precision :: snapmatrix(istate,nrsnapshots) 
      REAL sfactor
      integer iflagobs(0:ntimemax,mxpoi)
      REAL obsv(0:ntimemax,istate)

      REAL phi0(0:ntimemax,mxpoi)
      REAL u0(0:ntimemax,mxpoi)
      REAL v0(0:ntimemax,mxpoi)
      REAL w0(0:ntimemax,mxpoi)
      integer nsnap,itime
      common /fwdtraj/ phi0,u0,v0,w0            

!******************* snapmatrix construction

       
      if((iusemean.eq.1).or.(iuseobs.eq.1)) then
        do i=1,istate
          smean(i) = 0.0d0
        enddo
      endif
 
       if(iuseobs .eq. 1) then
         nsnap=0   
         do itime=0,ntime
!           if (3*int(itime/3) .eq. itime) then
!           if (5*int(itime/5) .eq. itime) then
           if (Tnsnap*int(itime/Tnsnap) .eq. itime) then
              nsnap = nsnap+1
              do i=1,istate
                 snapmatrix(i,nsnap) = obsv(itime,i)                 
              enddo              
           end if         
         enddo

       else 
         nsnap=0   
         do itime=0,ntime
!           if (3*int(itime/3) .eq. itime) then
!           if (5*int(itime/5) .eq. itime) then
           if (Tnsnap*int(itime/Tnsnap) .eq. itime) then
              nsnap = nsnap+1
              do i=1,mxpoi
                 snapmatrix(i,nsnap) = u0(itime,i)
                 snapmatrix(i+mxpoi,nsnap) = v0(itime,i) 
                 snapmatrix(i+2*mxpoi,nsnap) = w0(itime,i) 
                 snapmatrix(i+3*mxpoi,nsnap) = phi0(itime,i)                  
              enddo  


           end if
       
!              if(nsnap.eq.40) then
!                 ewrite(3,*) itime,nsnap
!                 ewrite(3,*) snapmatrix(100,nsnap),u0(itime,100)
!              endif
         enddo
       end if 

!        EWRITE(3,*) 'mx,mxpoi,istate,ntimemax,nrsnapshots',mx,mxpoi,istate,ntimemax,nrsnapshots

       if (nsnap .ne. nrsnapshots) then 
          write(*,*)'error:nsnap .ne. nrsnapshots:',nsnap,nrsnapshots
          stop
       end if  

       if (iusemean .ne. 0) then
          if((iusemean.eq.1).or.(iuseobs.eq.1)) then
           do i=1,istate
              do j=1,nrsnapshots
                smean(i) = smean(i) + snapmatrix(i,j) 
              enddo
           enddo
           do i=1,istate 
               smean(i) = smean(i)/nrsnapshots
           enddo 
          endif 
!************** subtract the mean value from the snapshots
          do i=1,istate 
             do j=1,nrsnapshots
                snapmatrix(i,j) = snapmatrix(i,j)-smean(i)
             enddo 
          enddo       
       end if 
          
       sfactor = dsqrt(dble(nrsnapshots))
       do i=1,istate 
          do j=1,nrsnapshots
                snapmatrix(i,j) = snapmatrix(i,j)/sfactor
          enddo 
       enddo 

       return
       end subroutine buildsnapmat

      subroutine buildsnapmat_dualwei(snapmatrix,smean,iuseobs,ntime,  &
                             obsv,iflagobs,iusemean,Tnsnap,dual_wei)

!********** builds snapshots matrix

      include 'paramnew.h'
      integer iuseobs,iusemean,ntime,i,j,Tnsnap
      REAL smean(istate),snapmatrix(istate,nrsnapshots) 
      REAL sfactor
      integer iflagobs(0:ntimemax,mxpoi)
      REAL obsv(0:ntimemax,istate)

      REAL phi0(0:ntimemax,mxpoi)
      REAL u0(0:ntimemax,mxpoi)
      REAL v0(0:ntimemax,mxpoi)
      REAL w0(0:ntimemax,mxpoi)

      REAL dual_wei(mxpoi)
      common /fwdtraj/ phi0,u0,v0,w0            

!******************* snapmatrix construction

       
      if((iusemean.eq.1).or.(iuseobs.eq.1)) then
        do i=1,istate
          smean(i) = 0.0d0
        enddo
      endif
 
       if(iuseobs .eq. 1) then
         nsnap=0   
         do itime=0,ntime
!           if (3*int(itime/3) .eq. itime) then
!           if (5*int(itime/5) .eq. itime) then
           if (Tnsnap*int(itime/Tnsnap) .eq. itime) then
              nsnap = nsnap+1
              do i=1,istate
                 snapmatrix(i,nsnap) = obsv(itime,i)                 
              enddo              
           end if         
         enddo

       else 
         nsnap=0   
         do itime=0,ntime
!           if (3*int(itime/3) .eq. itime) then
!           if (5*int(itime/5) .eq. itime) then
           if (Tnsnap*int(itime/Tnsnap) .eq. itime) then
              nsnap = nsnap+1
              do i=1,mxpoi
                 snapmatrix(i,nsnap) = u0(itime,i)
                 snapmatrix(i+mxpoi,nsnap) = v0(itime,i) 
                 snapmatrix(i+2*mxpoi,nsnap) = w0(itime,i) 
                 snapmatrix(i+3*mxpoi,nsnap) = phi0(itime,i)                  
              enddo  


           end if
       
!              if(nsnap.eq.40) then
!                 ewrite(3,*) itime,nsnap
!                 ewrite(3,*) snapmatrix(100,nsnap),u0(itime,100)
!              endif
         enddo
       end if 

!        EWRITE(3,*) 'mx,mxpoi,istate,ntimemax,nrsnapshots',mx,mxpoi,istate,ntimemax,nrsnapshots

       if (nsnap .ne. nrsnapshots) then 
          write(*,*)'error:nsnap .ne. nrsnapshots:',nsnap,nrsnapshots
          stop
       end if  

       if (iusemean .ne. 0) then
          if((iusemean.eq.1).or.(iuseobs.eq.1)) then
           do i=1,istate
              do j=1,nrsnapshots
                smean(i) = smean(i) + snapmatrix(i,j) 
              enddo
           enddo
           do j=1,nvar
           do i=1,mxpoi 
!               smean(i) = smean(i)/nrsnapshots
               smean((j-1)*mxpoi+i) = smean(j-1)*mxpoi+(i)*dual_wei(i)
           enddo
           enddo
          endif 
!************** subtract the mean value from the snapshots
          do i=1,istate 
             do j=1,nrsnapshots
                snapmatrix(i,j) = snapmatrix(i,j)-smean(i)
             enddo 
          enddo       
       end if 
          
!       sfactor = dsqrt(dble(nrsnapshots))
     do k=1,nvar
       do i=1,mxpoi 
          do j=1,nrsnapshots
                snapmatrix((k-1)*mxpoi+i,j) = snapmatrix((k-1)*mxpoi+i,j)*sqrt(dual_wei(i))
          enddo 
       enddo 
    enddo

       return
       end subroutine buildsnapmat_dualwei


!---------------------------------------------------------------------*
!This subroutine reads in the input file.
!The info read is:  the number of grid points in x and y (nx,ny),
!                   time step, final time, and time steps per plotting,
!                   p, q, alpha
!                   pstag, qstag.
!where .true. means that it is staggered and .false. means it is unstaggered.
!Written by F.X. Giraldo on 10/95
!---------------------------------------------------------------------*
      subroutine init(phii,ui,vi,node,coord,f,                               &
                     npoin,xmin,xmax,ymin,ymax,comega,nx,ny,dx,dy,dt,       &
                     ntime,rade,iplot,omega,alpha,velmax,cfl,p,q,alf,       &
                     pstag,qstag,idata)
      include 'paramnew.h'
      REAL coord(mxpoi,2)      
      REAL phii(mxpoi), ui(mxpoi), vi(mxpoi), f(mxpoi)
      REAL dummy(imax,jmax), ff(mxpoi)           
      integer node(imax,jmax), p, q, ntime, idata
      logical pstag,qstag
                                         
       
      nx = imax
      ny = jmax
      dt = 200
      time_final = 6      
      p = 4     !1
      q = 1
      alf =  1./3. !0.0
      pstag = .false.
      qstag = .false.
      
                                  !check bounds 
      if (nx.gt.imax.or.ny.gt.jmax) then 
         write(*,'(" Error! - Need to Enlarge IMAX and JMAX")')
         write(*,'(" nx ny imax jmax = ",4(i3,1x))')nx,ny,imax,jmax
         stop
      endif
                                  !Set some constants
      pi=4.0*atan(1.0)
      rade=6.37e+06
      time_final=time_final*3600.0
      ntime=nint(time_final/dt)                           
      if(ntime .gt. ntimemax) then
       write(6,*) 'maximum nr of time steps = ',ntimemax,' exceeded'
       write(6,*) 'ntime =',ntime
       stop
      end if        
      xmin=0.0
      xmax=2.0*pi
      ymin=-pi/2.0
      ymax=pi/2.0
      xl=xmax-xmin
      yl=ymax-ymin
      dx=xl/(nx)
      dy=yl/(ny)
      phi_mean=5.768e4
      omega=20.0
      comega=7.292e-05  
      velmax=-1e5
      alpha_fcor=0.0
      alpha=0.0
                                !set the Initial Conditions

      ip=0
      do j=1,ny
         olat=ymin + real(j-0.5)*dy
         do i=1,nx
            olon=xmin + real(i-0.5)*dx
            ip=ip+1
            node(i,j)=ip
            coord(ip,1)=olon
            coord(ip,2)=olat
            f(ip)=2.0*comega*( -cos(olon)*cos(olat)*sin(alpha_fcor) +        &
                               sin(olat)*cos(alpha_fcor) )
            if(sin(olat) .gt. 0) then
              ff(ip) = 2.0*comega*max(sin(olat),0.5)
            else
              ff(ip) = 2.0*comega*min(sin(olat),-0.5)
            end if        
            ui(ip)=omega*sin(olon)*(sin(olat)**3 -                           &
                                   3*sin(olat)*cos(olat)**2)
            vi(ip)=omega*sin(olat)**2*cos(olon)
            phii(ip)=phi_mean +                                              &
               2*comega*rade*omega*sin(olat)**3*cos(olat)*sin(olon)
         end do   
      end do        

      if (idata .ne. 0) then ! read file data
       open(10, file = 'h5x5.dat')
       do j=1,ny
          read(10,100) (dummy(i,j), i=1,nx)
       end do 
       close(10)
      
       ip=0
      
       do j=1,ny
        do i=1,nx
          ip = ip+1                              
          phii(ip) = dummy(i,j)*g          
        enddo
       enddo 
       
       if (idata .eq. 1) then 
           call geost(ui,vi,phii,nx,ny,ff,coord,dx,dy,rade)
       else
           open(10, file = 'u5x5.dat')
           do j=1,ny
              read(10,100) (dummy(i,j), i=1,nx)
           end do 
           close(10)
      
           ip=0
      
           do j=1,ny
             do i=1,nx
                ip = ip+1                              
                ui(ip) = dummy(i,j)          
             enddo
           enddo 

           open(10, file = 'v5x5.dat')
           do j=1,ny
              read(10,100) (dummy(i,j), i=1,nx)
           end do 
           close(10)
      
           ip=0
      
           do j=1,ny
              do i=1,nx
                 ip = ip+1                              
                 vi(ip) = dummy(i,j)          
              enddo
           enddo 

       end if             
      end if
    
      ip=0
      do j=1,ny  
         do i=1,nx
          ip=ip+1
          vel1=abs(ui(ip)) + abs(vi(ip)) + sqrt(2*phii(ip))
          velmax=max(velmax,vel1) 
         end do   
      end do  
                         
      dl=sqrt(dx**2 + dy**2)
      cfl=dt*velmax/(dl*rade)

      npoin=nx*ny

!      ewrite(3,*) ' dt dx dy velmax = ',dt,dx,dy,velmax
!      ewrite(3,*) ' ** CFL = ',cfl

100   format(500(1X, E16.8))      
101   format(100(1X, I2))  
                 
      return
      end subroutine init

!*---------------------------------------------------------------------*
!*This subroutine reads in the input file.
!*The info read is:  the number of grid points in x and y (nx,ny),
!*                   time step, final time, and time steps per plotting,
!*                   p, q, alpha
!*                   pstag, qstag.
!*where .true. means that it is staggered and .false. means it is unstaggered.
!*Written by F.X. Giraldo on 10/95
!*---------------------------------------------------------------------

!*
      subroutine initpar(node,coord,f,                              &
                     npoin,comega,nx,ny,dx,dy,dt,                  &
                     ntime,rade,omega,alpha,p,q,alf,               &
                     pstag,qstag)
      include 'paramnew.h'
      
      REAL f(mxpoi)
      REAL dt,dx,dy,rade,comega,alpha,alf
      REAL coord(mxpoi,2), xmin,xmax, ymin, ymax
      integer node(imax,jmax),p,q, nx,ny,npoin,ntime      
      logical pstag, qstag

                                        
       
      nx = imax
      ny = jmax
      dt = 200
      time_final = 6      
      p = 4     !1
      q = 1
      alf =  1./3. !0.0
      pstag = .false.
      qstag = .false.
      
                                  !check bounds 
      if (nx.gt.imax.or.ny.gt.jmax) then 
         write(*,'(" Error! - Need to Enlarge IMAX and JMAX")')
         write(*,'(" nx ny imax jmax = ",4(i3,1x))')nx,ny,imax,jmax
         stop
      endif
                                  !Set some constants
      pi=4.0*atan(1.0)
      rade=6.37e+06
      time_final=time_final*3600.0
      ntime=nint(time_final/dt)                           
      if(ntime .gt. ntimemax) then
       write(6,*) 'maximum nr of time steps = ',ntimemax,' exceeded'
       write(6,*) 'ntime =',ntime
       stop
      end if        
      xmin=0.0
      xmax=2.0*pi
      ymin=-pi/2.0
      ymax=pi/2.0
      xl=xmax-xmin
      yl=ymax-ymin
      dx=xl/(nx)
      dy=yl/(ny)
      phi_mean=5.768e4
      omega=20.0
      comega=7.292e-05  
      velmax=-1e5
      alpha_fcor=0.0
      alpha=0.0
                                !set the Initial Conditions

      ip=0
      do j=1,ny
         olat=ymin + real(j-0.5)*dy
         do i=1,nx
            olon=xmin + real(i-0.5)*dx
            ip=ip+1
            node(i,j)=ip
            coord(ip,1)=olon
            coord(ip,2)=olat
            f(ip)=2.0*comega*( -cos(olon)*cos(olat)*sin(alpha_fcor) +        &
                               sin(olat)*cos(alpha_fcor) )                   
         end do   
      end do        
                            
      npoin=nx*ny
                 
      return
      end subroutine initpar






      subroutine  geost(u0,v0,phi0,nx,ny,ff,coord,dx,dy,rade) 
      include 'paramnew.h'
      REAL u0(mxpoi),v0(mxpoi),phi0(mxpoi)
      REAL ff(mxpoi),coord(mxpoi,2)
      integer nx,ny
      REAL dummy(imax,jmax),dphidlong(imax,jmax)
      REAL dphidlat(imax,jmax)     
      
      
      ip = 0
      do j=1,ny
        do i=1,nx
          ip = ip+1                    
          dummy(i,j) = phi0(ip)          
        enddo
      enddo      
      
      ip = 0
      do j=1,ny
        do i=1,nx
          ip = ip+1
            if (i.gt. 1 .and. i .lt. nx) then
            dphidlong(i,j) = (dummy(i+1,j) - dummy(i-1,j))/ (2.*dx)
            elseif (i.eq. 1) then
            dphidlong(i,j) = (dummy(i+1,j) - dummy(i,j))/(dx)             
            else
              dphidlong (i,j)= (dummy(i,j) - dummy(i-1,j))/(dx)               
            end if
            if (j.gt. 1 .and. j .lt. ny) then
             dphidlat(i,j) = (dummy(i,j+1)- dummy(i,j-1) ) /(2.*dy)
            elseif(j .eq. 1) then
             dphidlat(i,j) =  (dummy(i,j+1)- dummy(i,j) ) /(dy)
            else
               dphidlat(i,j) =  (dummy(i,j)- dummy(i,j-1) ) /(dy)
            end if             
            u0(ip) = -dphidlat(i,j)/(rade*ff(ip))            
            v0(ip) = dphidlong(i,j)/(rade*ff(ip)*cos(coord(ip,2)))   
! cos(coord(ip,2)) may be large at Eq., so we will smooth later                        
        enddo
      enddo
                         
      ip = 0
      do j=1,ny
        do i=1,nx
          ip = ip+1 
          if(i .eq. nx/2) then
            v0(ip) = 0.5*(v0(ip-2)+v0(ip+2))
          elseif (i .eq. nx/2 +1) then
             v0(ip) = 0.5*(v0(ip-1)+v0(ip+1)) 
          end if
          if( j.eq. ny/2) then
            v0(ip)  = 0.5*(v0(ip-2*nx)+v0(ip+2*nx))
            u0(ip)  = 0.5*(u0(ip-2*nx)+u0(ip+2*nx))
          elseif( j.eq. ny/2+1) then
            v0(ip)  = 0.5*(v0(ip-2*nx)+v0(ip+2*nx))
            u0(ip)  = 0.5*(u0(ip-2*nx)+u0(ip+2*nx))
          endif   
        enddo
      enddo
      
      return
      end subroutine  geost
      
      
      subroutine perturb(phi0,u0,v0,nx,ny)
      include 'paramnew.h'
      REAL phi0(mxpoi),u0(mxpoi),v0(mxpoi)
      integer nx,ny

      double precision  xmin,xmax, ymin, ymax
      double precision dx,dy,rade,comega,omega,phi_mean

      go to 10

      ip = 0
      do j=1,ny
        do i=1,nx 
          ip = ip+1
          if (j.gt. 1 .and. j.lt. ny ) then
           rpert = rand() 
           phi0(ip) = phi0(ip)  + (0.5-rpert)*phi0(ip)*0.02d0           
           rpert = rand()
           u0(ip) = u0(ip) + (0.5-rpert)* u0(ip)*0.1d0           
           rpert = rand()
           v0(ip) =  v0(ip) + (0.5-rpert)* v0(ip)*0.1d0 
          end if      
        enddo
      enddo

10    continue
      pi=4.0*atan(1.0)
      rade=6.37e+06
      xmin=0.0
      xmax=2.0*pi
      ymin=-pi/2.0
      ymax=pi/2.0
      xl=xmax-xmin
      yl=ymax-ymin
      dx=xl/(nx)
      dy=yl/(ny)
      phi_mean=5.768e4
      omega=20.0
      comega=7.292e-05  

!cccc  set initial guess by shifting the state one grid point
!cccc  in longitudinal direction

      ip=0
      do j=1,ny
         olat=ymin + real(j-0.5)*dy
         do i=1,nx
            olon=xmin + real(i-0.5)*dx
            ip=ip+1
            u0(ip)=omega*sin(olon-dx)*(sin(olat)**3 -                           &
                                   3*sin(olat)*cos(olat)**2)
            v0(ip)=omega*sin(olat)**2*cos(olon-dx)
            phi0(ip)=phi_mean +                                                 &
               2*comega*rade*omega*sin(olat)**3*cos(olat)*sin(olon-dx)
         end do   
      end do        

           
      return
      end subroutine perturb


      subroutine backinclude(backterm,varin,bweights,cost,advarin)

      implicit none

      include 'paramnew.h'

      integer i
      double precision  backterm(istate), bweights(istate)
      double precision varin(istate), cost, advarin(istate)

      do i=1,istate 
         cost = cost + bweights(i)*(varin(i)-backterm(i))**2
      enddo

      do i=1,istate
         advarin(i) = advarin(i) +                                       &
                     2.0d0*bweights(i)*(varin(i)-backterm(i))
      enddo

      return
      end subroutine backinclude

!****************************************************************

      subroutine setobsv(ntime,u0,v0,w0,phi0,obsv,iflagobs,iobs)

      include 'paramnew.h'

      integer ntime, iobs, itime
      REAL u0(0:ntimemax,mxpoi),v0(0:ntimemax,mxpoi),w0(0:ntimemax,mxpoi)
      REAL phi0(0:ntimemax,mxpoi)
      integer iflagobs(0:ntimemax,mxpoi)
      REAL obsv(0:ntimemax,istate)

      do itime=0,ntime
         do ip=1,mxpoi
           obsv(itime,ip) = u0(itime,ip)
           obsv(itime,ip+mxpoi) = v0(itime,ip)
           obsv(itime,ip+2*mxpoi) = w0(itime,ip)
           obsv(itime,ip+3*mxpoi) = phi0(itime,ip)
           if (iobs*int(itime/iobs) .eq. itime) then  
              iflagobs(itime,ip) = 1
           else
              iflagobs(itime,ip) = 0
           end if      
         enddo
      enddo   
 
      return
      end subroutine setobsv


      subroutine pertobsv(ntime,obsv,iflagobs)

      implicit none

      include 'paramnew.h'
      integer ntime,itime,ip
      integer iflagobs(0:ntimemax,mxpoi)
      double precision obsv(0:ntimemax,istate)
      double precision rpert

       open(10, file = 'gaussrandm0std1.dat')
        do itime=0,ntime
         do ip=1,mxpoi
           if (iflagobs(itime,ip) .eq. 1) then
             read(10,50)rpert
             obsv(itime,ip) = obsv(itime,ip)  + 1.d-1*rpert
             read(10,50)rpert
             obsv(itime,ip+mxpoi) = obsv(itime,ip+mxpoi) + 1.d-1*rpert
             read(10,50)rpert
             obsv(itime,ip+2*mxpoi) = obsv(itime,ip+2*mxpoi) + 1.d0*rpert
             if(mxpoi*itime .gt. 80000+ip) rewind(unit = 10)
            end if
          enddo
        enddo
        close(10)
50    format(1x,E16.8)

      return
      end subroutine pertobsv


      subroutine writeout(ntime,nx,ny,u0,v0,phi0,file1,file2,file3)

      implicit none

      include 'paramnew.h'

      double precision u0(0:ntimemax,mxpoi),v0(0:ntimemax,mxpoi) 
      double precision phi0(0:ntimemax,mxpoi)        

      integer ntime, ip,j,nx,ny
      character*12 file1,file2,file3

        open(10, file = file1)
        do j=1,ny
         write(10,100)(u0(0,(j-1)*nx + ip), ip=1,nx)
        enddo
        do j=1,ny
         write(10,100)(u0(ntime,(j-1)*nx + ip), ip=1,nx)
        enddo        
        close(10)
        open(20, file = file2)
        do j=1,ny
         write(20,100)(v0(0,(j-1)*nx + ip), ip=1,nx)
        enddo
        do j=1,ny
         write(20,100)(v0(ntime,(j-1)*nx + ip), ip=1,nx)
        enddo        
        close(20)
        open(30, file = file3)
        do j=1,ny
         write(30,100)(phi0(0,(j-1)*nx + ip), ip=1,nx)
        enddo
        do j=1,ny
         write(30,100)(phi0(ntime,(j-1)*nx + ip), ip=1,nx)
        enddo        
        close(30)                
100   format(100(1X, E16.8))

      return
      end subroutine writeout


      subroutine setweights(npoin,weights,bweights)

      include 'paramnew.h'  

      integer npoin,ip
      double precision weights(istate), bweights(istate)    

      do ip = 1,npoin
          weights(ip) = 0.5d0
          weights(ip+npoin) = 0.5d0 
          weights(ip+2*npoin) = 1.0d0/5.768d4
      enddo 
      
      do ip = 1,istate
          bweights(ip) = weights(ip) !*1.d-2
      enddo

!        open(10, file = 'gaussrandm0std1.dat')
!        do i=1,istate
!           read(10,50)rpert
!           backterm(i) = varin(i) + 1.d-1*rpert ! define background
!           varin(i) = backterm(i)
!           bweights(i) = weights(i)*1.d-2
!        enddo     
!        close(10) 

      return
      end subroutine setweights

      subroutine buildsnapmat_PG(PGNODS,snapmatrix,smean,ntime,  &
                             iflagobs,iusemean,Tnsnap,PG)
!********** builds snapshots matrix

      include 'paramnew.h'
      integer iuseobs,iusemean,ntime,i,j,Tnsnap,PGNODS
      REAL smean(PGNODS),snapmatrix(PGNODS,nrsnapshots) 
      REAL sfactor
      REAL PG(0:ntimemax,PGNODS)

!******************* snapmatrix construction

       
      if(iusemean.eq.1) then
        do i=1,istate
          smean(i) = 0.0d0
        enddo
      endif
 
         nsnap=0   
         do itime=0,ntime
           if (Tnsnap*int(itime/Tnsnap) .eq. itime) then
              nsnap = nsnap+1
              do i=1,PGNODS
                 snapmatrix(i,nsnap) = PG(itime,i)
              enddo  
    
          end if 
         enddo

!        EWRITE(3,*) 'mx,PGNODS,ntimemax,nrsnapshots',mx,PGNODS,ntimemax,nrsnapshots

       if (nsnap .ne. nrsnapshots) then 
          write(*,*)'error:nsnap .ne. nrsnapshots:',nsnap,nrsnapshots
          stop
       end if  

       if (iusemean .ne. 0) then
          if((iusemean.eq.1).or.(iuseobs.eq.1)) then
           do i=1,PGNODS
              do j=1,nrsnapshots
                smean(i) = smean(i) + snapmatrix(i,j) 
              enddo
           enddo
           do i=1,PGNODS
               smean(i) = smean(i)/nrsnapshots
           enddo 
          endif 
!************** subtract the mean value from the snapshots
          do i=1,PGNODS
             do j=1,nrsnapshots
                snapmatrix(i,j) = snapmatrix(i,j)-smean(i)
             enddo 
          enddo       
       end if 
          
       sfactor = dsqrt(dble(nrsnapshots))
       do i=1,PGNODS
          do j=1,nrsnapshots
                snapmatrix(i,j) = snapmatrix(i,j)/sfactor
          enddo 
       enddo 

       return
       end subroutine buildsnapmat_PG

SUBROUTINE INTERPL_MESH(NONODS,TOTELE,NLOC,MLOC,NDGLNO,PNDGLN,FREDOP,NPRESS,X,Y,Z,U,V,W,P,&
                     LTIME,ACCTIM,TIMMES,DT)

  implicit none

  include 'paramnew.h'
  INTEGER NONODS,TOTELE,NLOC,MLOC,FREDOP,NPRESS
  INTEGER NDGLNO(TOTELE*NLOC)
  INTEGER PNDGLN(TOTELE*MLOC)
  REAL LTIME,ACCTIM,TIMMES,DT
  REAL X(NONODS),Y(NONODS),Z(NONODS)
  REAL U(NONODS),V(NONODS),W(NONODS),P(FREDOP*NPRESS)
  REAL phi0(0:ntimemax,mxpoi)
  REAL u0(0:ntimemax,mxpoi)
  REAL v0(0:ntimemax,mxpoi)
  REAL w0(0:ntimemax,mxpoi)
  common /fwdtraj/ phi0,u0,v0,w0

!LOCAL VARIABLES
  INTEGER NONODS1,TOTELE1,NLOC1,MLOC1,FREDOP1,NPRESS1
  INTEGER NFIELDS
  REAL,DIMENSION(:),ALLOCATABLE:: RLOCAL,R
  REAL,DIMENSION(:),ALLOCATABLE:: X1,Y1,Z1
  INTEGER,DIMENSION(:),ALLOCATABLE:: NDGLNO1,PNDGLN1
  INTEGER FIELDS2(4),FIELDS1(4)

  INTEGER   ::IERROR
  INTEGER II,JJ,KK,ISTATE_AD


  IF(mxpoi.NE.NONODS) THEN
!    ewrite(3,*) 'mxpoi.ne.nonods'
!    stop 356
  ENDIF

IF( TIMMES.GT.LTIME) THEN
   DO II=1,NONODS
       U0(NINT((ACCTIM-T0)/DT),II)=U(II)
       V0(NINT((ACCTIM-T0)/DT),II)=V(II)
       W0(NINT((ACCTIM-T0)/DT),II)=W(II)
       PHI0(NINT((ACCTIM-T0)/DT),II)=P(II)
   ENDDO

ELSE
   IF(NINT( (ACCTIM-T0)/DT).EQ.0) THEN
     OPEN(1,file='interpol.dat')
       WRITE(1,*) TOTELE,NONODS,NLOC,MLOC,FREDOP,NPRESS
       WRITE(1,*) (X(II),II=1,NONODS)
       WRITE(1,*) (Y(II),II=1,NONODS)
       WRITE(1,*) (Z(II),II=1,NONODS)
       WRITE(1,*) (NDGLNO(II),II=1,TOTELE*NLOC)
       WRITE(1,*) (PNDGLN(II),II=1,TOTELE*MLOC)
     CLOSE(1)
   ENDIF

  
  OPEN(1,file='interpol.dat')
  
  READ(1,*) TOTELE1,NONODS1,NLOC1,MLOC1,FREDOP1,NPRESS1
!  ewrite(3,*)  TOTELE1,NONODS1,NLOC1,MLOC1,FREDOP1,NPRESS1

  ALLOCATE(NDGLNO1(TOTELE1*NLOC1))
  ALLOCATE(PNDGLN1(TOTELE1*MLOC1))
  ALLOCATE(X1(NONODS1))
  ALLOCATE(Y1(NONODS1))
  ALLOCATE(Z1(NONODS1))

  READ(1,*) (X1(JJ),JJ=1,NONODS1)
  READ(1,*) (Y1(JJ),JJ=1,NONODS1)
  READ(1,*) (Z1(JJ),JJ=1,NONODS1)
!  ewrite(3,*) 'finish reading locations, x1,y1,z1'
  READ(1,*) (NDGLNO1(JJ),JJ=1,TOTELE1*NLOC1)
  READ(1,*) (PNDGLN1(JJ),JJ=1,TOTELE1*MLOC1)
  CLOSE(1)

  ALLOCATE(R(3*NONODS+FREDOP*NPRESS))
  ALLOCATE(RLOCAL(3*NONODS1+FREDOP1*NPRESS1))
  R=0.0
  RLOCAL=0.0

     ! set up the reference mesh0 system
      !........................................................

       FIELDS1(1) = 1
       FIELDS1(2) = 1*NONODS1+1
       FIELDS1(3) = 2*NONODS1+1
       FIELDS1(4) = 3*NONODS1+1

       RLOCAL=0.0
       !set up the adaptive mesh system
       FIELDS2(1) = 1
       FIELDS2(2) = 1*NONODS+1
       FIELDS2(3) = 2*NONODS+1
       FIELDS2(4) = 3*NONODS+1

       R(1:NONODS)=U(1:NONODS)
       R(1*NONODS+1:2*NONODS)=V(1:NONODS)
       R(2*NONODS+1:3*NONODS)=W(1:NONODS)

     ! interpolete U,V,W,P into the refernce mesh0
     !...................................................................

       NFIELDS = 3
       CALL FLTetra4toTetra4(NONODS,TOTELE,X,Y,Z,NDGLNO,R,FIELDS2,          &
          NFIELDS,NONODS1,TOTELE1,X1,Y1,Z1,NDGLNO1,RLOCAL,FIELDS1,IERROR)
          U0(NINT((ACCTIM-T0)/DT),1:NONODS1)= RLOCAL(1:NONODS1)
          V0(NINT((ACCTIM-T0)/DT),1:NONODS1)= RLOCAL(1*NONODS1+1:2*NONODS1)
          W0(NINT((ACCTIM-T0)/DT),1:NONODS1)= RLOCAL(2*NONODS1+1:3*NONODS1)

       NFIELDS = 1
       FIELDS1(1) = 1
       FIELDS2(1) = 1
       R(1:FREDOP*NPRESS)=P(1:FREDOP*NPRESS)
       CALL FLTetra4toTetra4(FREDOP*NPRESS,TOTELE,X,Y,Z,PNDGLN,R,FIELDS2,          &
          NFIELDS,FREDOP1*NPRESS1,TOTELE1,X1,Y1,Z1,PNDGLN1,RLOCAL,FIELDS1,IERROR)

          PHI0(NINT((ACCTIM-T0)/DT),1:NONODS1)= RLOCAL(1:FREDOP1*NPRESS1)

  DEALLOCATE(R)
  DEALLOCATE(RLOCAL)

  DEALLOCATE(NDGLNO1)
  DEALLOCATE(PNDGLN1)
  DEALLOCATE(X1)
  DEALLOCATE(Y1)
  DEALLOCATE(Z1)

ENDIF

END SUBROUTINE INTERPL_MESH

SUBROUTINE CalHessian(NONODS,XNONOD,NLOC,NGI,TOTELE,NDIM,T,HESSIAN,         &
                        X,Y,Z,XONDGL,NDGLNO,N,NLX,NLY,NLZ,DCYL,D3,weight)

implicit none

INTEGER, INTENT(IN)        :: NONODS,XNONOD,NLOC,NGI,TOTELE,NDIM
LOGICAL, INTENT(IN)        :: DCYL,D3
REAL, INTENT(IN),DIMENSION(NONODS)        :: T
REAL, INTENT(IN),DIMENSION(NONODS)        :: X,Y,Z
INTEGER, INTENT(IN),DIMENSION(TOTELE*NLOC) :: NDGLNO,XONDGL
REAL , INTENT(IN),DIMENSION(NLOC,NGI)        ::N,NLX,NLY,NLZ

REAL, INTENT(INOUT),DIMENSION(NONODS*NDIM*NDIM)        :: HESSIAN
                                
! local
REAL        :: DETWEI(NGI),WEIGHT(NGI),VOLUME
REAL,DIMENSION(NLOC,NGI)        :: NX,NY,NZ
REAL,DIMENSION(NONODS)        :: TX,TY,TZ,TXX,TXY,TXZ,TYY,TYZ,TZZ
REAL,DIMENSION(NONODS)        :: ML
INTEGER                        ::ELE,ILOC,JLOC,GLOBI,GLOBJ,II,JJ,KK,GI

       ewrite(3,*) 'in Cal_Hessian'
       write(3,*) 'dd'
       TX=0.0
       TY=0.0
       TZ=0.0
       ML=0.0
       DO 340 ELE=1,TOTELE
!
! Calculate DETWEI,RA,NX,NY,NZ for element ELE
           CALL DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI,      &    
                N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL,               &
                NX,NY,NZ) 
!
        DO ILOC=1,NLOC
                GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
                DO JLOC=1,NLOC
                GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
                   DO GI=1,NGI
                            TX(GLOBI)=TX(GLOBI)+N(ILOC,GI)*NX(JLOC,GI)*DETWEI(GI)*T(GLOBJ)
                         TY(GLOBI)=TY(GLOBI)+N(ILOC,GI)*NY(JLOC,GI)*DETWEI(GI)*T(GLOBJ)
                         TZ(GLOBI)=TZ(GLOBI)+N(ILOC,GI)*NZ(JLOC,GI)*DETWEI(GI)*T(GLOBJ)
                         ML(GLOBI)=ML(GLOBI)+N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
!                         write(3,*) GLOBI,GLOBJ,TX(GLOBI),N(ILOC,GI),NX(JLOC,GI),DETWEI(GI),T(GLOBJ)
!                         stop 14
                   END DO
                END DO
        END DO
340   CONTINUE
        
           DO GLOBI=1,NONODS
                 TX(GLOBI)=TX(GLOBI)/ML(GLOBI)
                 TY(GLOBI)=TY(GLOBI)/ML(GLOBI)
                 TZ(GLOBI)=TZ(GLOBI)/ML(GLOBI)
           END DO
                
        
       TXX=0.0
       TYY=0.0
       TZZ=0.0
       TXY=0.0
       TXZ=0.0
       TYZ=0.0
       DO 360 ELE=1,TOTELE
!
        DO ILOC=1,NLOC
                GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
                DO JLOC=1,NLOC
                GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
                   DO GI=1,NGI
                         TXX(GLOBI)=TXX(GLOBI)+N(ILOC,GI)*NX(JLOC,GI)*DETWEI(GI)*TX(GLOBJ)
                         TYY(GLOBI)=TYY(GLOBI)+N(ILOC,GI)*NY(JLOC,GI)*DETWEI(GI)*TY(GLOBJ)
                         TZZ(GLOBI)=TZZ(GLOBI)+N(ILOC,GI)*NZ(JLOC,GI)*DETWEI(GI)*TZ(GLOBJ)
                         
                         TXY(GLOBI)=TXY(GLOBI)+N(ILOC,GI)*NX(JLOC,GI)*DETWEI(GI)*TY(GLOBJ)
                         TXZ(GLOBI)=TXZ(GLOBI)+N(ILOC,GI)*NX(JLOC,GI)*DETWEI(GI)*TZ(GLOBJ)
                         TYZ(GLOBI)=TYZ(GLOBI)+N(ILOC,GI)*NY(JLOC,GI)*DETWEI(GI)*TZ(GLOBJ)
                    END DO
                END DO
        END DO
                        
360   CONTINUE
        
   DO GLOBI=1,NONODS
         TXX(GLOBI)=TXX(GLOBI)/ML(GLOBI)
         TYY(GLOBI)=TYY(GLOBI)/ML(GLOBI)
         TZZ(GLOBI)=TZZ(GLOBI)/ML(GLOBI)
         TXY(GLOBI)=TXY(GLOBI)/ML(GLOBI)
         TXZ(GLOBI)=TXZ(GLOBI)/ML(GLOBI)
         TYZ(GLOBI)=TYZ(GLOBI)/ML(GLOBI)
   END DO
open(1,file='Txyz1a.dat')
  write(1,*) (TXX(II),II=1,NONODS)
  write(1,*) (Tyy(II),II=1,NONODS)
  write(1,*) (Tzz(II),II=1,NONODS)
  write(1,*) (TXY(II),II=1,NONODS)
  write(1,*) (TXZ(II),II=1,NONODS)
  write(1,*) (Tyz(II),II=1,NONODS)
close(1)

        DO GLOBI=1,NONODS
                HESSIAN( (GLOBI-1)*NDIM*NDIM+1) = TXX(GLOBI)
                HESSIAN( (GLOBI-1)*NDIM*NDIM+2) = TXY(GLOBI)
                HESSIAN( (GLOBI-1)*NDIM*NDIM+3) = TXZ(GLOBI)
                HESSIAN( (GLOBI-1)*NDIM*NDIM+4) = TXY(GLOBI)
                HESSIAN( (GLOBI-1)*NDIM*NDIM+5) = TYY(GLOBI)
                HESSIAN( (GLOBI-1)*NDIM*NDIM+6) = TYZ(GLOBI)
                HESSIAN( (GLOBI-1)*NDIM*NDIM+7) = TXZ(GLOBI)
                HESSIAN( (GLOBI-1)*NDIM*NDIM+8) = TYZ(GLOBI)
                HESSIAN( (GLOBI-1)*NDIM*NDIM+9) = TZZ(GLOBI)
        ENDDO
                
END SUBROUTINE CalHessian



SUBROUTINE READ_SHAPE_XYZ_MUPTXX(TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS,        & 
                        M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,WEIGHT,FINDRM,CENTRM,                        &
                        DT,LTIME,ACCTIM,TIMMES,D3,DCYL,                                                &
                        U,V,W,P,RMEM,IMEM,NRMEM,NIMEM,                                                &
                        GEOBAL,SCFACTH0,                &                        
                        MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,X,Y,Z,NDGLNO,XONDGL,PNDGLN,    &
                        stotel,snloc,sndgln)

  implicit none

  include 'paramnew.h'
  INTEGER  TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS
  INTEGER  stotel,snloc,sndgln
  LOGICAL  D3,DCYL
  INTEGER  M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,WEIGHT,FINDRM,CENTRM
  REAL           DT,LTIME,ACCTIM,TIMMES
  INTEGER  NRMEM,NIMEM
  REAL           RMEM(NRMEM)
  INTEGER  IMEM(NIMEM)
  INTEGER  U,V,W,P
  INTEGER  GEOBAL
  REAL           SCFACTH0
  INTEGER  MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,X,Y,Z,NDGLNO,XONDGL,PNDGLN
  INTEGER  II,JJ,KK,I,J,K

  IF(NINT((ACCTIM-T0)/DT).EQ.0) THEN
        OPEN(1,file='shape.dat')
           
           WRITE(1,*) TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS
           WRITE(1,*) D3,DCYL
           
           DO II=1,NGI
              WRITE(1,*) (RMEM(M+(II-1)*NLOC+I-1),I=1,MLOC)
              WRITE(1,*) (RMEM(MLX+(II-1)*NLOC+I-1),I=1,MLOC)
              WRITE(1,*) (RMEM(MLY+(II-1)*NLOC+I-1),I=1,MLOC)
              WRITE(1,*) (RMEM(MLZ+(II-1)*NLOC+I-1),I=1,MLOC)
           ENDDO
           
           DO II=1,NGI
              WRITE(1,*) (RMEM(N+(II-1)*NLOC+I-1),I=1,NLOC)
              WRITE(1,*) (RMEM(NLX+(II-1)*NLOC+I-1),I=1,NLOC)
              WRITE(1,*) (RMEM(NLY+(II-1)*NLOC+I-1),I=1,NLOC)
              WRITE(1,*) (RMEM(NLZ+(II-1)*NLOC+I-1),I=1,NLOC)
           ENDDO
           
           WRITE(1,*) (RMEM(WEIGHT+I-1),I=1,NGI)

           WRITE(1,*) (IMEM(FINDRM-1+I),I=1,NONODS+1)
           WRITE(1,*) (IMEM(CENTRM-1+I),I=1,NONODS)

        CLOSE(1)
           
        OPEN(1,file='XYZ_MUPTXX.dat')
           WRITE(1,*) DT,LTIME
           !WRITE(1,*) OPTOME,GEOBAL
           FLAbort("Broken - ask Stephan")
           !WRITE(1,*) OMEGA,OMEGA1,OMEGA2,OMEGA3,OMEGA4,SCFACTH0
           WRITE(1,*) (RMEM(MUPTXX+I-1),I=1,NONODS)
           WRITE(1,*) (RMEM(MUPTXY+I-1),I=1,NONODS)
           WRITE(1,*) (RMEM(MUPTXZ+I-1),I=1,NONODS)
           WRITE(1,*) (RMEM(MUPTYY+I-1),I=1,NONODS)
           WRITE(1,*) (RMEM(MUPTYZ+I-1),I=1,NONODS)
           WRITE(1,*) (RMEM(MUPTZZ+I-1),I=1,NONODS)
           WRITE(1,*) (RMEM(X+I-1),I=1,NONODS)
           WRITE(1,*) (RMEM(Y+I-1),I=1,NONODS)
           WRITE(1,*) (RMEM(Z+I-1),I=1,NONODS)
           WRITE(1,*) (IMEM(NDGLNO+I-1),I=1,TOTELE*NLOC)
           WRITE(1,*) (IMEM(XONDGL+I-1),I=1,TOTELE*NLOC)
           WRITE(1,*) (IMEM(PNDGLN+I-1),I=1,TOTELE*MLOC)
        CLOSE(1)
  ENDIF
        OPEN(1,file='surface-shape.dat')
           WRITE(1,*) stotel,snloc
           WRITE(1,*) (IMEM(sndgln+I-1),I=1,stotel*snloc)
        CLOSE(1)
           ewrite(3,*)'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'
           ewrite(3,*)'NINT( ACCTIM/DT)',NINT(ACCTIM/DT),ACCTIM,LTIME
!!call read_N(NLOC,MLOC,NGI,RMEM(N),RMEM(NLX),RMEM(NLY),RMEM(NLZ),RMEM(M),RMEM(MLX),RMEM(MLY),RMEM(MLZ))

  IF( NINT( (ACCTIM-T0)/DT).LE.ntimemax) THEN
        CALL INTERPL_MESH(NONODS,TOTELE,NLOC,MLOC,IMEM(NDGLNO),IMEM(PNDGLN),FREDOP,NPRESS, &
              RMEM(X),RMEM(Y),RMEM(Z),RMEM(U),RMEM(V),RMEM(W),RMEM(P),                           &
              LTIME,ACCTIM,TIMMES,DT)
  ENDIF

END SUBROUTINE READ_SHAPE_XYZ_MUPTXX

SUBROUTINE READ_SHAPE_XYZ_MUPTXX2(TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS,        & 
                        M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,WEIGHT,FINDRM,CENTRM,                        &
                        D3,DCYL,DT,                                                &
                        GEOBAL,SCFACTH0,                &                        
                        MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,X,Y,Z,NDGLNO,XONDGL,PNDGLN,    &
                        stotel,snloc,sndgln)

  implicit none

  include 'paramnew.h'
  INTEGER  TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS
  INTEGER  stotel,snloc
  REAL sndgln(stotel*snloc)
  REAL DT,LTIME
  LOGICAL  D3,DCYL
  REAL,DIMENSION(MLOC,NGI)::   M,MLX,MLY,MLZ
  REAL,DIMENSION(NLOC,NGI):: N,NLX,NLY,NLZ
  REAL WEIGHT(NGI)
  INTEGER FINDRM(NONODS+1),CENTRM(NONODS)
  REAL P(FREDOP*NPRESS)
  INTEGER  GEOBAL
  REAL           SCFACTH0
  REAL,DIMENSION(NONODS):: MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,X,Y,Z
  INTEGER,DIMENSION(TOTELE*NLOC):: NDGLNO,XONDGL,PNDGLN
  INTEGER  II,JJ,KK,I,J,K

        LTIME=100.0
        OPEN(1,file='shape.dat')
           
           WRITE(1,*) TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS
           WRITE(1,*) D3,DCYL
           
           DO II=1,NGI
              WRITE(1,*) (M(I,II),I=1,MLOC)
              WRITE(1,*) (MLX(I,II),I=1,MLOC)
              WRITE(1,*) (MLY(I,II),I=1,MLOC)
              WRITE(1,*) (MLZ(I,II),I=1,MLOC)
           ENDDO
           
           DO II=1,NGI
              WRITE(1,*) (N(I,II),I=1,NLOC)
              WRITE(1,*) (NLX(I,II),I=1,NLOC)
              WRITE(1,*) (NLY(I,II),I=1,NLOC)
              WRITE(1,*) (NLZ(I,II),I=1,NLOC)
           ENDDO
           
           WRITE(1,*) (WEIGHT(I),I=1,NGI)

           WRITE(1,*) (FINDRM(I),I=1,NONODS+1)
           WRITE(1,*) (CENTRM(I),I=1,NONODS)

        CLOSE(1)
           
        OPEN(1,file='XYZ_MUPTXX.dat')
           WRITE(1,*) DT,LTIME
           !WRITE(1,*) OPTOME,GEOBAL
           FLAbort("Broken - ask Stephan")
           !WRITE(1,*) OMEGA,OMEGA1,OMEGA2,OMEGA3,OMEGA4,SCFACTH0
           WRITE(1,*) (MUPTXX(I),I=1,NONODS)
           WRITE(1,*) (MUPTXY(I),I=1,NONODS)
           WRITE(1,*) (MUPTXZ(I),I=1,NONODS)
           WRITE(1,*) (MUPTYY(I),I=1,NONODS)
           WRITE(1,*) (MUPTYZ(I),I=1,NONODS)
           WRITE(1,*) (MUPTZZ(I),I=1,NONODS)
           WRITE(1,*) (X(I),I=1,NONODS)
           WRITE(1,*) (Y(I),I=1,NONODS)
           WRITE(1,*) (Z(I),I=1,NONODS)
           WRITE(1,*) (NDGLNO(I),I=1,TOTELE*NLOC)
           WRITE(1,*) (XONDGL(I),I=1,TOTELE*NLOC)
           WRITE(1,*) (PNDGLN(I),I=1,TOTELE*MLOC)
        CLOSE(1)

        OPEN(1,file='surface-shape.dat')
           WRITE(1,*) stotel,snloc
           WRITE(1,*) (sndgln(I),I=1,stotel*snloc)
        CLOSE(1)


END SUBROUTINE READ_SHAPE_XYZ_MUPTXX2

SUBROUTINE read_N(NLOC,MLOC,NGI,N,NLX,NLY,NLZ,M,MLX,MLY,MLZ)
       INTEGER NLOC,MLOC,NGI
       REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
       REAL M(NLOC,NGI),MLX(NLOC,NGI),MLY(NLOC,NGI),MLZ(NLOC,NGI)
       integer ii,jj,kk
open(1,file='NN.dat')
  DO JJ=1,NGI
     write(1,*) (M(II,JJ),II=1,MLOC)
     write(1,*) (MLX(II,JJ),II=1,MLOC)
     write(1,*) (MLY(II,JJ),II=1,MLOC)
     write(1,*) (MLZ(II,JJ),II=1,MLOC)
  ENDDO
  DO JJ=1,NGI
     write(1,*) (N(II,JJ),II=1,MLOC)
     write(1,*) (NLX(II,JJ),II=1,MLOC)
     write(1,*) (NLY(II,JJ),II=1,MLOC)
     write(1,*) (NLZ(II,JJ),II=1,MLOC)
  ENDDO
close(1)
END SUBROUTINE read_N

SUBROUTINE READ_OPTIMISED_INITIAL(NONODS,NOBCU,NOBCV,NOBCW,ACCTIM,        &
        U,V,W,P,BCU1,BCV1,BCW1,BCU2,BCV2,BCW2)

  implicit none

  include 'paramnew.h'

  INTEGER NONODS,NOBCU,NOBCV,NOBCW
  REAL ACCTIM
  REAL U(NONODS),V(NONODS),W(NONODS),P(NONODS)
  REAL BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
  INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
  REAL phi0(0:ntimemax,mxpoi)
  REAL u0(0:ntimemax,mxpoi)
  REAL v0(0:ntimemax,mxpoi)
  REAL w0(0:ntimemax,mxpoi)
  common /fwdtraj/ phi0,u0,v0,w0
  INTEGER I,J,K,II,JJ,KK
        IF(.false.) THEN
        IF(ABS(ACCTIM-T0).LT.1.0E-8) THEN
           DO II=1,NONODS
              U(II) = U0(0,II)
              V(II) = V0(0,II)
              W(II) = W0(0,II)
              P(II) = PHI0(0,II)
           ENDDO
           if(.true.) then
           DO II=1,NOBCU
              I=BCU2(II)
              U(I)=BCU1(II)
              U0(0,I)=U(I)
           ENDDO
           DO II=1,NOBCV
              I=BCV2(II)
              V(I)=BCV1(II)
              V0(0,I)=V(I)
           ENDDO
           DO II=1,NOBCW
              I=BCW2(II)
              W(I)=BCW1(II)
              W0(0,I)=W(I)
           ENDDO
           endif
        ENDIF
        ENDIF
END SUBROUTINE READ_OPTIMISED_INITIAL

SUBROUTINE H1_SPACE(my_epsilon,my_epsilon2,snapmatrix,smean,istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots,LEFTSVD,SVDVAL)
  
  implicit none
  
  INTEGER :: istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots
  real :: smean(istate)
  double precision :: snapmatrix(istate,nrsnapshots)
  double precision :: LEFTSVD(istate,nsvd), svdval(nsvd_total)
  REAL    :: my_epsilon,my_epsilon2
  
  
  INTEGER :: TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS
  LOGICAL :: D3,DCYL
  REAL,DIMENSION(:),ALLOCATABLE     ::    MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ
  REAL,DIMENSION(:),ALLOCATABLE     ::    WEIGHT,DETWEI
  INTEGER,DIMENSION(:),ALLOCATABLE  ::    NDGLNO,XONDGL,PNDGLN,SNDGLN0
  REAL,DIMENSION(:),ALLOCATABLE     ::    X,Y,Z,NU,NV,NW,UG,VG,WG
  REAL,DIMENSION(:,:),ALLOCATABLE   ::    N,NLX,NLY,NLZ,M,MLX,MLY,MLZ
  REAL,DIMENSION(:,:),ALLOCATABLE   ::    NX,NY,NZ,MX,MY,MZ
  
  ! Local variabels
  INTEGER :: II,JJ,KK,I,J,K,GI,ILOC,JLOC
  REAL DT,LTIME
  INTEGER GEOBAL
  REAL SCFACTH0
  
  
  OPEN(1,file='shape.dat')
  
  READ(1,*) TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS
  READ(1,*) D3,DCYL
  ewrite(3,*) 'TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS',  &
       TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC,FREDOP,NPRESS
  ewrite(3,*) D3,DCYL
  ALLOCATE(N(NLOC,NGI))
  ALLOCATE(NLX(NLOC,NGI))
  ALLOCATE(NLY(NLOC,NGI))
  ALLOCATE(NLZ(NLOC,NGI))
  
  ALLOCATE(M(MLOC,NGI))
  ALLOCATE(MLX(MLOC,NGI))
  ALLOCATE(MLY(MLOC,NGI))
  ALLOCATE(MLZ(MLOC,NGI))
  ALLOCATE(NX(NLOC,NGI))
  ALLOCATE(NY(NLOC,NGI))
  ALLOCATE(NZ(NLOC,NGI))
  ALLOCATE(MX(MLOC,NGI))
  ALLOCATE(MY(MLOC,NGI))
  ALLOCATE(MZ(MLOC,NGI))
  
  II=TOTELE*NLOC
  ALLOCATE(NDGLNO(II))
  ALLOCATE(XONDGL(II))
  II=TOTELE*MLOC
  ALLOCATE(PNDGLN(II))
  
  ewrite(3,*) 'after locating M,MLX..'
  ALLOCATE(MUPTXX(MXPOI))
  ALLOCATE(MUPTXY(MXPOI))
  ALLOCATE(MUPTXZ(MXPOI))
  ALLOCATE(MUPTYY(MXPOI))
  ALLOCATE(MUPTYZ(MXPOI))
  ALLOCATE(MUPTZZ(MXPOI))
  ewrite(3,*) 'after locating the viscosities'
  ALLOCATE(WEIGHT(NGI))
  ALLOCATE(X(MXPOI))
  ALLOCATE(Y(MXPOI))
  ALLOCATE(Z(MXPOI))
  
  N=0.0
  NLX=0.0
  NLY=0.0
  NLZ=0.0
  M=0.0
  MLX=0.0
  MLY=0.0
  MLZ=0.0
  NX=0.0
  NY=0.0
  NZ=0.0
  MX=0.0
  MY=0.0
  MZ=0.0

  MUPTXX=0.0
  MUPTXY=0.0
  MUPTXZ=0.0
  MUPTYY=0.0
  MUPTYZ=0.0
  MUPTZZ=0.0
  WEIGHT=0.0
  X=0.0
  Y=0.0
  Z=0.0
  NDGLNO=0
  XONDGL=0
  PNDGLN=0
  
  
  DO JJ=1,NGI
     READ(1,*) (M(II,JJ),II=1,MLOC)
     READ(1,*) (MLX(II,JJ),II=1,MLOC)
     READ(1,*) (MLY(II,JJ),II=1,MLOC)
     READ(1,*) (MLZ(II,JJ),II=1,MLOC)
  ENDDO
  
  ewrite(3,*)'finish reading M,MLX...'
  
  DO JJ=1,NGI
     READ(1,*) (N(II,JJ),II=1,NLOC)
     READ(1,*) (NLX(II,JJ),II=1,NLOC)
     READ(1,*) (NLY(II,JJ),II=1,NLOC)
     READ(1,*) (NLZ(II,JJ),II=1,NLOC)
  ENDDO
  READ(1,*) (WEIGHT(JJ),JJ=1,NGI)

  CLOSE(1)
  
  
  OPEN(1,file='XYZ_MUPTXX.dat')
  READ(1,*) DT,LTIME
  FLAbort("Broken - Stephan")
  !READ(1,*) OPTOME,GEOBAL
  !READ(1,*) OMEGA,OMEGA1,OMEGA2,OMEGA3,OMEGA4,SCFACTH0
  ewrite(3,*) 'DT,LTIME',DT,LTIME
  
  READ(1,*) (MUPTXX(JJ),JJ=1,NONODS)
  READ(1,*) (MUPTXY(JJ),JJ=1,NONODS)
  READ(1,*) (MUPTXZ(JJ),JJ=1,NONODS)
  READ(1,*) (MUPTYY(JJ),JJ=1,NONODS)
  READ(1,*) (MUPTYZ(JJ),JJ=1,NONODS)
  READ(1,*) (MUPTZZ(JJ),JJ=1,NONODS)
  ewrite(3,*) 'finishing reading viscosity'
  READ(1,*) (X(JJ),JJ=1,NONODS)
  READ(1,*) (Y(JJ),JJ=1,NONODS)
  READ(1,*) (Z(JJ),JJ=1,NONODS)
  ewrite(3,*) 'finish reading locations, x,y,z'
  READ(1,*) (NDGLNO(JJ),JJ=1,TOTELE*NLOC)
  READ(1,*) (XONDGL(JJ),JJ=1,TOTELE*NLOC)
  READ(1,*) (PNDGLN(JJ),JJ=1,TOTELE*MLOC)

  CLOSE(1)


  CALL DERIVATIVE_LEFTSVD(LEFTSVD,SVDVAL,snapmatrix, &
     istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots,X,Y,Z,      &
     M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,     &
     MX,MY,MZ,NX,NY,NZ,D3,DCYL,my_epsilon,my_epsilon2, &
     NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,NONODS,XNONOD)


  DEALLOCATE(N)
  DEALLOCATE(NLX)
  DEALLOCATE(NLY)
  DEALLOCATE(NLZ)
  ewrite(3,*) 'after locating N,NLX..'
  DEALLOCATE(M)
  DEALLOCATE(MLX)
  DEALLOCATE(MLY)
  DEALLOCATE(MLZ)
  DEALLOCATE(NX)
  DEALLOCATE(NY)
  DEALLOCATE(NZ)
  DEALLOCATE(MX)
  DEALLOCATE(MY)
  DEALLOCATE(MZ)
  II=TOTELE*NLOC
  DEALLOCATE(NDGLNO)
  DEALLOCATE(XONDGL)
  II=TOTELE*MLOC
  DEALLOCATE(PNDGLN)

  ewrite(3,*) 'after locating M,MLX..'
  DEALLOCATE(MUPTXX)
  DEALLOCATE(MUPTXY)
  DEALLOCATE(MUPTXZ)
  DEALLOCATE(MUPTYY)
  DEALLOCATE(MUPTYZ)
  DEALLOCATE(MUPTZZ)
  ewrite(3,*) 'after locating the viscosities'
  DEALLOCATE(WEIGHT)
  DEALLOCATE(X)
  DEALLOCATE(Y)
  DEALLOCATE(Z)

END SUBROUTINE H1_SPACE

SUBROUTINE DERIVATIVE_LEFTSVD(LEFTSVD,SVDVAL,snapmatrix, &
     istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots,X,Y,Z,      &
     M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,     &
     MX,MY,MZ,NX,NY,NZ,D3,DCYL,my_epsilon,my_epsilon2,     &
     NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,NONODS,XNONOD)

  implicit none

  INTEGER,INTENT(IN) :: istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots
  double precision, INTENT(IN),DIMENSION(istate,nrsnapshots) :: snapmatrix
  double precision :: LEFTSVD(istate,nsvd), svdval(nsvd_total)
  
  INTEGER,INTENT(IN) :: TOTELE,NLOC,NGI,MLOC,NONODS,XNONOD
  REAL,INTENT(IN),DIMENSION(NONODS) :: X,Y,Z
  INTEGER, INTENT(IN),DIMENSION(TOTELE*NLOC):: NDGLNO,XONDGL
  INTEGER, INTENT(IN),DIMENSION(TOTELE*MLOC):: PNDGLN
  REAL,INTENT(IN),DIMENSION(NLOC,NGI) :: N,NLX,NLY,NLZ
  REAL,INTENT(IN),DIMENSION(MLOC,NGI) :: M,MLX,MLY,MLZ
  REAL,INTENT(INOUT),DIMENSION(NLOC,NGI) :: NX,NY,NZ
  REAL,INTENT(INOUT),DIMENSION(MLOC,NGI) :: MX,MY,MZ
  LOGICAL :: D3,DCYL
  LOGICAL :: IF_X,IF_Y,IF_Z
  REAL :: WEIGHT(NGI),DETWEI(NGI)
  REAL :: my_epsilon,my_epsilon2
  
  double precision, dimension(:, :), allocatable :: snapmatrix_dx, snapmatrix_dy, snapmatrix_dz
  double precision, dimension(:, :), allocatable :: snapmatrix_dxx, snapmatrix_dyy,snapmatrix_dzz
  double precision, dimension(:, :), allocatable :: leftsvd_dx
  REAL,DIMENSION(:),ALLOCATABLE :: svdval_dx
  
  ! Local variabels
  INTEGER :: II,JJ,KK,I,J,K,GI,ILOC,JLOC,ELE
  
  ALLOCATE(snapmatrix_dx(istate,nrsnapshots))
  ALLOCATE(snapmatrix_dy(istate,nrsnapshots))
  ALLOCATE(snapmatrix_dz(istate,nrsnapshots))
  snapmatrix_dx=0.0
  snapmatrix_dy=0.0
  snapmatrix_dz=0.0
  
! print*,'my_epsilon2=',my_epsilon2
 
  CALL DERIVATIVE_SNAPMATRIX(.true.,.false.,.false.,snapmatrix, snapmatrix_dx,       &
     istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots,X,Y,Z,  &
     M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,     &
     MX,MY,MZ,NX,NY,NZ,D3,DCYL,     &
     NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,NONODS,XNONOD)
  CALL DERIVATIVE_SNAPMATRIX(.false.,.true.,.false.,snapmatrix, snapmatrix_dy,       &
     istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots,X,Y,Z,  &
     M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,     &
     MX,MY,MZ,NX,NY,NZ,D3,DCYL,     &
     NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,NONODS,XNONOD)
IF(D3) THEN
  CALL DERIVATIVE_SNAPMATRIX(.false.,.false.,.true.,snapmatrix, snapmatrix_dz,       &
     istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots,X,Y,Z,  &
     M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,     &
     MX,MY,MZ,NX,NY,NZ,D3,DCYL,     &
     NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,NONODS,XNONOD)
ENDIF

IF(abs(my_epsilon2).GT.1.0E-12) then
  ALLOCATE(snapmatrix_dxx(istate,nrsnapshots))
  ALLOCATE(snapmatrix_dyy(istate,nrsnapshots))
  ALLOCATE(snapmatrix_dzz(istate,nrsnapshots))
  snapmatrix_dxx=0.0
  snapmatrix_dyy=0.0
  snapmatrix_dzz=0.0
  CALL DERIVATIVE_SNAPMATRIX(.true.,.false.,.false.,snapmatrix_dx, snapmatrix_dxx,       &
       istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots,X,Y,Z,  &
       M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,     &
       MX,MY,MZ,NX,NY,NZ,D3,DCYL,     &
       NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,NONODS,XNONOD)
  CALL DERIVATIVE_SNAPMATRIX(.false.,.true.,.false.,snapmatrix_dy, snapmatrix_dyy,       &
       istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots,X,Y,Z,  &
       M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,     &
       MX,MY,MZ,NX,NY,NZ,D3,DCYL,     &
       NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,NONODS,XNONOD)
  IF(D3) THEN
     CALL DERIVATIVE_SNAPMATRIX(.false.,.false.,.true.,snapmatrix_dz, snapmatrix_dzz,       &
          istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots,X,Y,Z,  &
          M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,     &
          MX,MY,MZ,NX,NY,NZ,D3,DCYL,     &
          NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,NONODS,XNONOD)
  ENDIF
ENDIF
IF(abs(my_epsilon2).GT.1.0E-12) then
   OPEN(10,FILE='snapmatrix_dxx.dat')
   do ii=1,istate
      WRITE(10, 2000) (snapmatrix_dxx(ii,jj), jj=1,nrsnapshots)
   enddo
   close(10)
ELSE
   OPEN(10,FILE='snapmatrix_dx.dat')
   do ii=1,istate
      WRITE(10, 2000) (snapmatrix_dx(ii,jj), jj=1,nrsnapshots)
   enddo
   close(10)
ENDIF
OPEN(10,FILE='snapmatrix_dy.dat')
do ii=1,istate
   WRITE(10, 2000) (snapmatrix_dy(ii,jj), jj=1,nrsnapshots)
enddo
close(10)
OPEN(10,FILE='snapmatrix_dz.dat')
do ii=1,istate
   WRITE(10, 2000) (snapmatrix_dz(ii,jj), jj=1,nrsnapshots)
enddo
close(10)
!stop 1
2000 format (20(1X,E24.16)) 

IF(abs(my_epsilon2).GT.1.0E-12) then
   ewrite(3,*) 'calculation of POD vectors for U in H1 space'
   call snapsvd_H1(mxpoi,nrsnapshots,my_epsilon2,snapmatrix(1:mxpoi,:),  &
        snapmatrix_dxx(1:mxpoi,:),snapmatrix_dyy(1:mxpoi,:),snapmatrix_dzz(1:mxpoi,:), &
        nsvd,nsvd,leftsvd(1:mxpoi,:),svdval(1:nsvd),D3)
   
   ewrite(3,*) 'calculation of POD vectors for V in H1 space'
   call snapsvd_H1(mxpoi,nrsnapshots,my_epsilon2,snapmatrix(1*mxpoi+1:2*mxpoi,:), &
        snapmatrix_dxx(1*mxpoi+1:2*mxpoi,:),snapmatrix_dyy(1*mxpoi+1:2*mxpoi,:), &
        snapmatrix_dzz(1*mxpoi+1:2*mxpoi,:),nsvd,nsvd,leftsvd(1*mxpoi+1:2*mxpoi,:),svdval(nsvd+1:2*nsvd),D3)
   
   IF(D3) THEN
      ewrite(3,*) 'calculation of POD vectors for W in H1 space'
      call snapsvd_H1(mxpoi,nrsnapshots,my_epsilon2,snapmatrix(2*mxpoi+1:3*mxpoi,:), &
           snapmatrix_dxx(2*mxpoi+1:3*mxpoi,:),snapmatrix_dyy(2*mxpoi+1:3*mxpoi,:),           &
           snapmatrix_dzz(2*mxpoi+1:3*mxpoi,:),nsvd,nsvd,leftsvd(2*mxpoi+1:3*mxpoi,:),svdval(2*nsvd+1:3*nsvd),D3)
   ENDIF
   
   ewrite(3,*) 'calculation of POD vectors for PHI in H1 space'
   call snapsvd_H1(mxpoi,nrsnapshots,my_epsilon2,snapmatrix(3*mxpoi+1:4*mxpoi,:), &
        snapmatrix_dxx(3*mxpoi+1:4*mxpoi,:),snapmatrix_dyy(3*mxpoi+1:4*mxpoi,:),&
        snapmatrix_dzz(3*mxpoi+1:4*mxpoi,:),nsvd,nsvd,leftsvd(3*mxpoi+1:4*mxpoi,:),svdval(3*nsvd+1:4*nsvd),D3)
   
ELSE
   ewrite(3,*) 'calculation of POD vectors for U in H1 space'
   call snapsvd_H1(mxpoi,nrsnapshots,my_epsilon,snapmatrix(1:mxpoi,:),  &
        snapmatrix_dx(1:mxpoi,:),snapmatrix_dy(1:mxpoi,:),snapmatrix_dz(1:mxpoi,:), &
        nsvd,nsvd,leftsvd(1:mxpoi,:),svdval(1:nsvd),D3)
   
   ewrite(3,*) 'calculation of POD vectors for V in H1 space'
   call snapsvd_H1(mxpoi,nrsnapshots,my_epsilon,snapmatrix(1*mxpoi+1:2*mxpoi,:), &
        snapmatrix_dx(1*mxpoi+1:2*mxpoi,:),snapmatrix_dy(1*mxpoi+1:2*mxpoi,:), &
        snapmatrix_dz(1*mxpoi+1:2*mxpoi,:),nsvd,nsvd,leftsvd(1*mxpoi+1:2*mxpoi,:),svdval(nsvd+1:2*nsvd),D3)
   
   IF(D3) THEN
      ewrite(3,*) 'calculation of POD vectors for W in H1 space'
      call snapsvd_H1(mxpoi,nrsnapshots,my_epsilon,snapmatrix(2*mxpoi+1:3*mxpoi,:), &
           snapmatrix_dx(2*mxpoi+1:3*mxpoi,:),snapmatrix_dy(2*mxpoi+1:3*mxpoi,:),           &
           snapmatrix_dz(2*mxpoi+1:3*mxpoi,:),nsvd,nsvd,leftsvd(2*mxpoi+1:3*mxpoi,:),svdval(2*nsvd+1:3*nsvd),D3)
   ENDIF
   
   ewrite(3,*) 'calculation of POD vectors for PHI in H1 space'
   call snapsvd_H1(mxpoi,nrsnapshots,my_epsilon,snapmatrix(3*mxpoi+1:4*mxpoi,:), &
        snapmatrix_dx(3*mxpoi+1:4*mxpoi,:),snapmatrix_dy(3*mxpoi+1:4*mxpoi,:),&
        snapmatrix_dz(3*mxpoi+1:4*mxpoi,:),nsvd,nsvd,leftsvd(3*mxpoi+1:4*mxpoi,:),svdval(3*nsvd+1:4*nsvd),D3)
ENDIF

deallocate(snapmatrix_dx)
deallocate(snapmatrix_dy)
deallocate(snapmatrix_dz)
IF(abs(my_epsilon2).GT.1.0E-12) then
   DEALLOCATE(snapmatrix_dxx)
   DEALLOCATE(snapmatrix_dyy)
   DEALLOCATE(snapmatrix_dzz)
ENDIF
END SUBROUTINE DERIVATIVE_LEFTSVD

SUBROUTINE DERIVATIVE_SNAPMATRIX(IF_X,IF_Y,IF_Z,snapmatrix,snapmatrix_dx, &
     istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots,X,Y,Z,  &
     M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,     &
     MX,MY,MZ,NX,NY,NZ,D3,DCYL,     &
     NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,NONODS,XNONOD)

  implicit none

  INTEGER,INTENT(IN) :: istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots
  double precision, INTENT(IN),DIMENSION(istate,nrsnapshots) :: snapmatrix
  double precision, INTENT(INOUT),DIMENSION(istate,nrsnapshots) :: snapmatrix_dx

  INTEGER,INTENT(IN) :: TOTELE,NLOC,NGI,MLOC,NONODS,XNONOD
  REAL,INTENT(IN),DIMENSION(NONODS) :: X,Y,Z
  INTEGER, INTENT(IN),DIMENSION(TOTELE*NLOC):: NDGLNO,XONDGL
  INTEGER, INTENT(IN),DIMENSION(TOTELE*MLOC):: PNDGLN
  REAL,INTENT(IN),DIMENSION(NLOC,NGI) :: N,NLX,NLY,NLZ
  REAL,INTENT(IN),DIMENSION(MLOC,NGI) :: M,MLX,MLY,MLZ
  REAL,INTENT(INOUT),DIMENSION(NLOC,NGI) :: NX,NY,NZ
  REAL,INTENT(INOUT),DIMENSION(MLOC,NGI) :: MX,MY,MZ
  LOGICAL :: D3,DCYL
  LOGICAL :: IF_X,IF_Y,IF_Z
  REAL :: WEIGHT(NGI),DETWEI(NGI)
! Local variabels
  INTEGER :: II,JJ,KK,I,J,K,GI,ILOC,JLOC,ELE,IGL,JGL,IGLP,JGLP
    REAL VOLUME
  REAL,DIMENSION(:),ALLOCATABLE :: ml

  ALLOCATE(ml(nonods))
    snapmatrix_dx = 0.0
    ml=0.0
     DO ELE=1,TOTELE  !ELE LOOP
         CALL DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,MLOC,NGI,      &    
                M,MLX,MLY,MLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL,               &
                MX,MY,MZ) 
         CALL DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI,      &    
                N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL,               &
                NX,NY,NZ) 

       DO ILOC = 1,NLOC !ILOC LOOP
           IGL=NDGLNO((ELE-1)*NLOC+ILOC)
         DO JLOC = 1,NLOC
           JGL=NDGLNO((ELE-1)*NLOC+JLOC)
           DO GI=1,NGI
             IF(IF_X) THEN
               snapmatrix_dx(IGL,:)=snapmatrix_dx(IGL,:)+ N(ILOC,GI)*DETWEI(GI)*snapmatrix(JGL,:)*NX(JLOC,GI)
               snapmatrix_dx(mxpoi+IGL,:)=snapmatrix_dx(mxpoi+IGL,:)+ N(ILOC,GI)*DETWEI(GI)*snapmatrix(mxpoi+JGL,:)*NX(JLOC,GI)
               snapmatrix_dx(2*mxpoi+IGL,:)=snapmatrix_dx(2*mxpoi+IGL,:)+ N(ILOC,GI)*DETWEI(GI)*snapmatrix(2*mxpoi+JGL,:)*NX(JLOC,GI)
             ELSE IF(IF_Y) THEN
               snapmatrix_dx(IGL,:)=snapmatrix_dx(IGL,:)+ N(ILOC,GI)*DETWEI(GI)*snapmatrix(JGL,:)*NY(JLOC,GI)
               snapmatrix_dx(mxpoi+IGL,:)=snapmatrix_dx(mxpoi+IGL,:)+ N(ILOC,GI)*DETWEI(GI)*snapmatrix(mxpoi+JGL,:)*NY(JLOC,GI)
               snapmatrix_dx(2*mxpoi+IGL,:)=snapmatrix_dx(2*mxpoi+IGL,:)+ N(ILOC,GI)*DETWEI(GI)*snapmatrix(2*mxpoi+JGL,:)*NY(JLOC,GI)
             ELSE IF(IF_Z) THEN
               snapmatrix_dx(IGL,:)=snapmatrix_dx(IGL,:)+ N(ILOC,GI)*DETWEI(GI)*snapmatrix(JGL,:)*NZ(JLOC,GI)
               snapmatrix_dx(mxpoi+IGL,:)=snapmatrix_dx(mxpoi+IGL,:)+ N(ILOC,GI)*DETWEI(GI)*snapmatrix(mxpoi+JGL,:)*NZ(JLOC,GI)
               snapmatrix_dx(2*mxpoi+IGL,:)=snapmatrix_dx(2*mxpoi+IGL,:)+ N(ILOC,GI)*DETWEI(GI)*snapmatrix(2*mxpoi+JGL,:)*NZ(JLOC,GI)
             ENDIF
              ml(igl)=ml(igl)+ N(ILOC,GI)*DETWEI(GI)*N(JLOC,GI)
           ENDDO
         ENDDO ! JLOC loop
       ENDDO !ILOC LOOP

       DO ILOC = 1,MLOC !ILOC LOOP
          IGLP=PNDGLN((ELE-1)*MLOC+ILOC)
         DO JLOC = 1,MLOC
          JGLP=PNDGLN((ELE-1)*MLOC+JLOC)
           DO GI=1,NGI
             IF(IF_X) THEN
              snapmatrix_dx(3*mxpoi+IGLP,:)=snapmatrix_dx(3*mxpoi+IGLP,:)+ M(ILOC,GI)*DETWEI(GI)*snapmatrix(3*mxpoi+JGLP,:)*MX(JLOC,GI)
             ELSE IF(IF_Y) THEN
              snapmatrix_dx(3*mxpoi+IGLP,:)=snapmatrix_dx(3*mxpoi+IGLP,:)+ M(ILOC,GI)*DETWEI(GI)*snapmatrix(3*mxpoi+JGLP,:)*MY(JLOC,GI)
             ELSE IF(IF_Z) THEN
              snapmatrix_dx(3*mxpoi+IGLP,:)=snapmatrix_dx(3*mxpoi+IGLP,:)+ M(ILOC,GI)*DETWEI(GI)*snapmatrix(3*mxpoi+JGLP,:)*MZ(JLOC,GI)
             ENDIF
           ENDDO
         ENDDO ! JLOC loop
        ENDDO !ILOC LOOP
     ENDDO !END ELE LOOP


           do igl=1,nonods
               snapmatrix_dx(IGL,:)=snapmatrix_dx(IGL,:)/ml(igl)
               snapmatrix_dx(mxpoi+IGL,:)=snapmatrix_dx(mxpoi+IGL,:)/ml(igl)
               snapmatrix_dx(2*mxpoi+IGL,:)=snapmatrix_dx(2*mxpoi+IGL,:)/ml(igl)
               snapmatrix_dx(3*mxpoi+IGL,:)=snapmatrix_dx(3*mxpoi+IGL,:)/ml(igl)
           end do
     
  DEALLOCATE(ml)

END SUBROUTINE DERIVATIVE_SNAPMATRIX

SUBROUTINE DERIVATIVE_SMEAN(IF_X,IF_Y,IF_Z,smean, &
     istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots,X,Y,Z,  &
     M,MLX,MLY,MLZ,N,NLX,NLY,NLZ,     &
     MX,MY,MZ,NX,NY,NZ,D3,DCYL,my_epsilon,    &
     NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,NONODS,XNONOD)

  implicit none

  INTEGER,INTENT(IN) :: istate,NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXPOI,nrsnapshots
  REAL, INTENT(INOUT),DIMENSION(istate) :: smean
  REAL my_epsilon

  INTEGER,INTENT(IN) :: TOTELE,NLOC,NGI,MLOC,NONODS,XNONOD
  REAL,INTENT(IN),DIMENSION(NONODS) :: X,Y,Z
  INTEGER, INTENT(IN),DIMENSION(TOTELE*NLOC):: NDGLNO,XONDGL
  INTEGER, INTENT(IN),DIMENSION(TOTELE*MLOC):: PNDGLN
  REAL,INTENT(IN),DIMENSION(NLOC,NGI) :: N,NLX,NLY,NLZ
  REAL,INTENT(IN),DIMENSION(MLOC,NGI) :: M,MLX,MLY,MLZ
  REAL,INTENT(IN),DIMENSION(NLOC,NGI) :: NX,NY,NZ
  REAL,INTENT(IN),DIMENSION(MLOC,NGI) :: MX,MY,MZ
  LOGICAL :: D3,DCYL
  LOGICAL :: IF_X,IF_Y,IF_Z
  REAL :: WEIGHT(NGI),DETWEI(NGI)
  REAL,DIMENSION(:),ALLOCATABLE :: smean_dx

! Local variabels
  INTEGER :: II,JJ,KK,I,J,K,GI,ILOC,JLOC,ELE,IGL,JGL,IGLP,JGLP
    REAL VOLUME

  ewrite(3,*) 'in DERIVATIVE_SMEAN'

  ALLOCATE(SMEAN_DX(ISTATE))
  SMEAN_DX=0.0

     DO ELE=1,TOTELE  !ELE LOOP

       DO ILOC = 1,NLOC !ILOC LOOP
           IGL=NDGLNO((ELE-1)*NLOC+ILOC)
           DO GI=1,NGI
             IF(IF_X) THEN
              smean_dx(IGL)=smean_dx(IGL) + smean(IGL)*NX(ILOC,GI)
               smean_dx(mxpoi+IGL)=smean_dx(mxpoi+IGL) + smean(mxpoi+IGL)*NX(ILOC,GI)
               smean_dx(2*mxpoi+IGL)=smean_dx(2*mxpoi+IGL)+smean(2*mxpoi+IGL)*NX(ILOC,GI)
             ELSE IF(IF_Y) THEN
               smean_dx(IGL)=smean_dx(IGL) + smean(IGL)*NY(ILOC,GI)
               smean_dx(mxpoi+IGL)=smean_dx(mxpoi+IGL) + smean(mxpoi+IGL)*NY(ILOC,GI)
               smean_dx(2*mxpoi+IGL)=smean_dx(2*mxpoi+IGL)+ smean(2*mxpoi+IGL)*NY(ILOC,GI)
             ELSE IF(IF_Z) THEN
               smean_dx(IGL)=smean_dx(IGL) + smean(IGL)*NZ(ILOC,GI)
               smean_dx(mxpoi+IGL)=smean_dx(mxpoi+IGL) + smean(mxpoi+IGL)*NZ(ILOC,GI)
               smean_dx(2*mxpoi+IGL)=smean_dx(2*mxpoi+IGL)+ smean(2*mxpoi+IGL)*NZ(ILOC,GI)
             ENDIF
           ENDDO
       ENDDO !ILOC LOOP
       DO ILOC = 1,MLOC !ILOC LOOP
           IGLP=PNDGLN((ELE-1)*MLOC+ILOC)
           DO GI=1,NGI
              smean_dx(3*mxpoi+IGLP)=smean_dx(3*mxpoi+IGLP)+ smean(3*mxpoi+IGL)*MX(ILOC,GI)
           ENDDO
       ENDDO !ILOC LOOP
     ENDDO !END ELE LOOP

     SMEAN = SMEAN+my_epsilon*SMEAN_DX
        OPEN(10,FILE='smean_dx.dat')
           WRITE(10, 2000) (smean_dx(ii), ii=1,istate)
        close(10)

  DEALLOCATE(SMEAN_DX)


2000 format (20(1X,E24.16)) 

END SUBROUTINE DERIVATIVE_SMEAN

          SUBROUTINE JACIMP1(HESIAN,VECERR,NNODP,NONODS,NDIM, &
! Working arrays...
           V1,V2,A,VV,   &
           DD,D1,D2)
! Combine the metrics in HESIAN,VECERR and put the result in VECERR
          
          implicit none
          
          INTEGER NNODP,NONODS,NDIM
          REAL HESIAN(NONODS*NDIM*NDIM),VECERR(NONODS*NDIM*NDIM)
!           REAL X(NONODS),Y(NONODS),Z(NONODS)
          REAL V1(NDIM,NDIM),V2(NDIM,NDIM),A(NDIM,NDIM)
          REAL VV(NDIM,NDIM)
          REAL DD(NDIM),D1(NDIM),D2(NDIM)
! Local variables...
          INTEGER NOD,NDIM2
          NDIM2=NDIM**2
          DO NOD=1,NNODP
            CALL JACIM1(HESIAN((NOD-1)*NDIM2+1), &
            VECERR((NOD-1)*NDIM2+1),VECERR((NOD-1)*NDIM2+1),NDIM,  &
! Working arrays...
            V1,V2,A,VV,DD,D1,D2,.true.)
          END DO
          RETURN

  end subroutine jacimp1

SUBROUTINE DUALWEIGHT_SPACE(dual_wei_space,Metric,NDIM,D3,X,Y,Z,   &
                  NONODS,TOTELE,NLOC,NCOLM,NDGLNO,FINDRM,COLM)

   implicit none

   INTEGER, INTENT(IN)  :: NDIM, NONODS,TOTELE, NLOC,NCOLM
   REAL, INTENT(IN), DIMENSION(NONODS) ::  X,Y,Z,Metric
   INTEGER, INTENT(IN), DIMENSION(TOTELE*NLOC) :: NDGLNO
   INTEGER, INTENT(IN), DIMENSION(NONODS+1) :: FINDRM
   INTEGER, INTENT(IN), DIMENSION(NCOLM) :: COLM
   REAL,INTENT(OUT), DIMENSION(NONODS) :: dual_wei_space
   LOGICAL :: D3

! Local variable

  INTEGER :: ELE,NOD,II,JJ,KK,I,J,K,NOD1,NOD2,NDIM2
  INTEGER :: ILOC,JLOC,KLOC,GI,COUNT,LOWER,UPPER,NL
  REAL LENGTH1,LENGTH,LENGTH_ave
  REAL VV(3)

   VV=0.0
   NDIM2 = NDIM*NDIM
   DO ELE=1,TOTELE
     DO NOD1=1,NONODS
       LOWER=FINDRM(NOD1)
       UPPER=FINDRM(NOD1+1)-1
       NL=0
       LENGTH =0.0       
       DO COUNT=LOWER, UPPER
          NOD2=COLM(COUNT)
          NL=NL+1
          VV(1)=X(NOD2)-X(NOD1)
          VV(2)=Y(NOD2)-Y(NOD1)
          IF(D3) VV(3)=Z(NOD2)-Z(NOD1)
          CALL JACED1(LENGTH1,VV,Metric((NOD1-1)*NDIM2+1),Metric((NOD2-1)*NDIM2+1),NDIM)        
          LENGTH = LENGTH +LENGTH1
       ENDDO
       dual_wei_space(NOD1) = LENGTH/NL
    ENDDO

   ENDDO


END SUBROUTINE DUALWEIGHT_SPACE

          SUBROUTINE JACED1(LENGTH,VV,M1,M2,NDIM) 
          
          implicit none
          
          INTEGER NDIM
          REAL LENGTH
          REAL VV(NDIM),M1(NDIM,NDIM),M2(NDIM,NDIM)
! Calculate LENGTH=(VV^T M VV)
          INTEGER I,J
          REAL VDUMI
!
          LENGTH=0.0
          DO I=1,NDIM
            VDUMI=0.0 
            DO J=1,NDIM
              VDUMI=VDUMI+0.5*(M1(I,J)+M2(I,J))*VV(J)
            END DO
            LENGTH=LENGTH+VV(I)*VDUMI
          END DO

          RETURN

  end subroutine jaced1

  SUBROUTINE EXASOL_UVW_WRITE(NONODS,U,V,W,X,Y,Z,ACCTIM,DT,D3)
    
    implicit none
    
    !----------------------------------------------------
    ! This subroutine is used to get the exact solutions
    !----------------------------------------------------
    
    INTEGER,PARAMETER ::NOSTIM=251,NOSNOD=2
    !  INTEGER,PARAMETER ::NOSTIM=21,NOSNOD=10
    INTEGER ::ITSTIME
    REAL    ::ACCTIM
    REAL               ::UEXAC(NOSNOD),VEXAC(NOSNOD),           &
         WEXAC(NOSNOD),HEXAC(NOSNOD)
    REAL    ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
    REAL    ::STIME(NOSTIM)
    INTEGER ::NONODS
    REAL    ::U(NONODS),V(NONODS),W(NONODS)
    REAL    ::X(NONODS),Y(NONODS),Z(NONODS)
    REAL    ::DT,LTIME
    LOGICAL ::D3
    CHARACTER(40) ::UNAME,VNAME,WNAME,SXNAME,SYNAME,SZNAME,STIMENAME
    REAL,PARAMETER ::PIE=3.141592654
    ! Local variables.......
    INTEGER ::II,KK,I,NOD
    REAL    ::DIST,A
    
    UEXAC=0.0
    VEXAC=0.0
    WEXAC=0.0
    SX=0.0
    SY=0.0
    SZ=0.0
    STIME=0.0
    
    
    DO KK =2,NOSTIM
       STIME(KK) = STIME(KK-1)+DT
    ENDDO
    
    
    SX(1) = 3.0
    SY(1) = 0.0
    SX(2) = 6.0
    SY(2) = 0.0
    
    DO II=1,NOSNOD
       DO I=1,NONODS
          A = (X(I)-SX(II))**2+(Y(I)-SY(II))**2
          IF( I.EQ.1) THEN
             DIST = A
             NOD=I
          ELSE
             IF(A.LT.DIST) THEN
                DIST = A
                NOD=I
             ENDIF
          ENDIF
       ENDDO
       UEXAC(II) = U(NOD)
       VEXAC(II) = V(NOD)
       IF(D3) WEXAC(II) = W(NOD)
    ENDDO
    
    
    
   IF( ABS(ACCTIM-0.0).LT.1.0E-6 ) THEN
      OPEN(1,FILE='AdjSourceData.dat')
      WRITE(1,*) 'SXSXSXSXSXSXSXSXSXSXSXSXSXSX'
      WRITE(1,*) (SX(II),II=1,NOSNOD)
      WRITE(1,*) 'SYSYSYSYSYSYSYSYSYSYSYSYSYSY'
      WRITE(1,*) (SY(II),II=1,NOSNOD)
      IF(D3) THEN
         WRITE(1,*) 'SZSZSZSZSZSZSZSZSZSZSZSZSZ'
         WRITE(1,*) (SZ(II),II=1,NOSNOD)
      ENDIF
      
      WRITE(1,*) 'STIMESTIMESTIMESTIMESTIMESTIME'
      WRITE(1,*) (STIME(II),II=1,NOSTIM)
   ELSE
      OPEN(1,FILE='AdjSourceData.dat',POSITION='APPEND')
   ENDIF
   
   WRITE(1,*)  ACCTIM
   WRITE(1,*) (UEXAC(II),II=1,NOSNOD)
   WRITE(1,*) (VEXAC(II),II=1,NOSNOD)
   IF(D3) THEN
      WRITE(1,*) (WEXAC(II),II=1,NOSNOD)
   ENDIF
   
   CLOSE(1)
   
   
 END SUBROUTINE EXASOL_UVW_WRITE


 SUBROUTINE CalViscosity(NONODS,XNONOD,NLOC,NGI,TOTELE,U,V,MUPTXX,MUPTYY,viscosity,           &
      X,Y,Z,XONDGL,NDGLNO,N,NLX,NLY,NLZ,DCYL,D3,weight)
   
   implicit none
   
   INTEGER, INTENT(IN)           :: NONODS,XNONOD,NLOC,NGI,TOTELE
   LOGICAL, INTENT(IN)           :: DCYL,D3
   REAL, INTENT(IN),DIMENSION(NONODS)           :: U,V,MUPTXX,MUPTYY
   REAL, INTENT(IN),DIMENSION(NONODS)           :: X,Y,Z
   INTEGER, INTENT(IN),DIMENSION(TOTELE*NLOC) :: NDGLNO,XONDGL
   REAL , INTENT(IN),DIMENSION(NLOC,NGI)           ::N,NLX,NLY,NLZ
   
   ! local
   REAL           :: DETWEI(NGI),WEIGHT(NGI),VOLUME
   REAL,DIMENSION(NLOC,NGI)           :: NX,NY,NZ
   REAL,DIMENSION(NONODS)           :: UX,UY,UZ,UXX,UYY,UZZ
   REAL,DIMENSION(NONODS)           :: VX,VY,VZ,VXX,VYY,VZZ,VISCOSITY
   REAL,DIMENSION(NONODS)           :: ML
   INTEGER                                 ::ELE,ILOC,JLOC,GLOBI,GLOBJ,II,JJ,KK,GI
   
   ewrite(3,*) 'in Cal_Hessian'
   print*,'dd'
   UX=0.0
   UY=0.0
   UZ=0.0
   VX=0.0
   VY=0.0
   VZ=0.0
   ML=0.0
   DO 340 ELE=1,TOTELE
      !
      
      ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
      CALL DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI,      &    
           N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL,               &
           NX,NY,NZ) 
      !
      DO ILOC=1,NLOC
         GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
         DO JLOC=1,NLOC
            GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
            DO GI=1,NGI
               UX(GLOBI)=UX(GLOBI)+N(ILOC,GI)*NX(JLOC,GI)*DETWEI(GI)*U(GLOBJ)
               UY(GLOBI)=UY(GLOBI)+N(ILOC,GI)*NY(JLOC,GI)*DETWEI(GI)*U(GLOBJ)
               VX(GLOBI)=VX(GLOBI)+N(ILOC,GI)*NX(JLOC,GI)*DETWEI(GI)*V(GLOBJ)
               VY(GLOBI)=VY(GLOBI)+N(ILOC,GI)*NY(JLOC,GI)*DETWEI(GI)*V(GLOBJ)
               ML(GLOBI)=ML(GLOBI)+N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
            END DO
         END DO
      END DO
340   CONTINUE
      
      DO GLOBI=1,NONODS
         UX(GLOBI)=UX(GLOBI)/ML(GLOBI)
         UY(GLOBI)=UY(GLOBI)/ML(GLOBI)
         VX(GLOBI)=VX(GLOBI)/ML(GLOBI)
         VY(GLOBI)=VY(GLOBI)/ML(GLOBI)
      END DO
      
      
      UXX=0.0
      UYY=0.0
      UZZ=0.0
      VXX=0.0
      VYY=0.0
      VZZ=0.0
      DO 360 ELE=1,TOTELE
         !
         DO ILOC=1,NLOC
            GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
            DO JLOC=1,NLOC
               GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
               DO GI=1,NGI
                  UXX(GLOBI)=UXX(GLOBI)+N(ILOC,GI)*NX(JLOC,GI)*DETWEI(GI)*UX(GLOBJ)
                  UYY(GLOBI)=UYY(GLOBI)+N(ILOC,GI)*NY(JLOC,GI)*DETWEI(GI)*UY(GLOBJ)
                  VXX(GLOBI)=VXX(GLOBI)+N(ILOC,GI)*NX(JLOC,GI)*DETWEI(GI)*VX(GLOBJ)
                  VYY(GLOBI)=VYY(GLOBI)+N(ILOC,GI)*NY(JLOC,GI)*DETWEI(GI)*VY(GLOBJ)
                  
               END DO
            END DO
         END DO
         
360      CONTINUE
         DO GLOBI=1,NONODS
            UXx(GLOBI)=UxX(GLOBI)/ML(GLOBI)
            UYy(GLOBI)=UYy(GLOBI)/ML(GLOBI)
            VxX(GLOBI)=VXx(GLOBI)/ML(GLOBI)
            VYy(GLOBI)=VYy(GLOBI)/ML(GLOBI)
         END DO
         
         VISCOSITY=0.0           
         DO GLOBI=1,NONODS
            VISCOSITY(GLOBI)=  SQRT((MUPTXX(GLOBI)*(UXX(GLOBI)+UYY(GLOBI)))**2+(MUPTYY(GLOBI)*(VXX(GLOBI)+VYY(GLOBI)))**2)
         END DO
         
         
       END SUBROUTINE CalViscosity

end module pod_support
