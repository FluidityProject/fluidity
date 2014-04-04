!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA


#include "fdebug.h"

   module ann_bp
     use spud
     implicit none
     ! sample_num: number of samples. here equals to timesteps .In : NO. of input variables of each sample. 
     !Out :NO. of output variables of each sample. 
     ! Neuron: NO.of neuron，Train  
     ! integer, save :: sample_num = 200  !Data
     integer, save :: in_num=2
     integer, save :: out_num=1
     integer, save :: Neuron =45  
     integer, save :: Train_num=200000
     real, save :: A=0.2
     real, save :: B=0.4
     real, save :: aa=0.2
     real, save :: bb=0.3
     integer :: timestep,nsvd,sample_num, coef_num
     ! real, dimension(:,:),allocatable ::d_in, d_out,w,v,dv,dw
     ! real, dimension(:),  allocatable :: OutputData,o !Maxin,Minin,Maxout,Minout
     ! double precision e;
     ! real Maxin(2),Minin(2),Maxout(1),Minout(1);
     real, dimension(:,:),allocatable :: d_in_allcoef ! store coef_pod_all
     real d_in(1600,2),d_out(1600,1);
     real w(45,2),o(45),v(1,45);
     real Maxin(2),Minin(2),Maxout(1),Minout(1);
     real OutputData(1);
     real dv(1,45),dw(45,2);
     real e;
      real, dimension(:,:), allocatable ::pod_coef_all_diff, pod_coef_all_ann,pod_coef_all_diff_obj
  
    contains 
     
     subroutine ann_bp_main(total_timestep,timestep,output,l)
      integer, intent(in) ::  total_timestep,timestep,l
      real, intent(inout) ::output
      integer :: i, k,m,n
      real, dimension(:),allocatable ::s
      real ::ss
      real :: coef, coefdiff
      integer :: ctrl_timestep,  ctrl_timestep1
       call get_option(&
            '/reduced_model/pod_basis_formation/pod_basis_count', nsvd)
       coef_num=3*nsvd        
       allocate(s(3*nsvd))
      ! sample_num=820!timestep
        sample_num=timestep-3

     
       allocate(pod_coef_all_diff(total_timestep,coef_num))
       allocate(pod_coef_all_diff_obj(total_timestep,coef_num))
       allocate(pod_coef_all_ann(total_timestep,coef_num))
       pod_coef_all_ann=0
       pod_coef_all_diff=0
       pod_coef_all_diff_obj=0
       print *, 'timestep',timestep
       open(1,file='coef_pod_all_obv')
         read(1,*)((pod_coef_all_ann(m,n),n=1,coef_num),m=1,timestep) !initialisation
       close(1)
       pod_coef_all_diff=pod_coef_all_ann
       do m=1, timestep-1
         pod_coef_all_diff(m,:)= pod_coef_all_ann(m+1,:)-pod_coef_all_ann(m,:)
       enddo
       do m=1, timestep-2
         pod_coef_all_diff_obj(m,:)= pod_coef_all_ann(m+2,:)-pod_coef_all_ann(m+1,:)
       enddo
  
       open(21,file='coef_pod_all_diff')
          write(21,*)((pod_coef_all_diff(m,n),n=1,coef_num),m=1,timestep) 
       close(21)
       open(21,file='coef_pod_all_diff_obj')
          write(21,*)((pod_coef_all_diff_obj(m,n),n=1,coef_num),m=1,timestep) 
       close(21)
      ! stop 93
       if (have_option("/reduced_model/training")) then
          ctrl_timestep=1500
          ctrl_timestep1=1500
         else 
          ctrl_timestep=555555
          ctrl_timestep1=19
       endif
          if(timestep .eq. ctrl_timestep) then  ! if run then eq.5555
	    ! call readData(l)
             call readtestData2(l)
           ! call readtestData() 
             call initBPNetwork() 
             call trainNetwork() 
             call writeweight(l)    
               coef= pod_coef_all_ann(timestep-2,l)
                
             !  coefdiff= pod_coef_all_diff(timestep-2,1)
             coefdiff= pod_coef_all_ann(timestep-1,l)
              call bpresultok(coef, coefdiff,ss)
              output=ss
             ! output=ss+pod_coef_all_ann(timestep-1,1)
              print*, 'coef+coefdiff=', output
             elseif(timestep >ctrl_timestep1) then 
             !  print *, 'coefbefore',pod_coef_all_diff(timestep-2,1)
             !  coef= pod_coef_all_ann(timestep-2,1)
             !  coefdiff= pod_coef_all_diff(timestep-2,1)
               call readweight(l)
               coef= pod_coef_all_ann(timestep-2,l)
               coefdiff= pod_coef_all_ann(timestep-1,l)
               print *, 'coef,coefdiff',coef,coefdiff
               call bpresultok(coef, coefdiff,ss)
               print *, 'coef+coefdiff=', ss
            ! output=ss+pod_coef_all_ann(timestep-1,1)
              output=ss
             ! print*, 'coef+coefdiff=', output
            !  stop 93			   
            endif
      
     
       !  open(40,file='pod_coef')
        ! write(40,*)(s(i),i=1,3*nsvd) 
       !  close(40)

      
      deallocate(pod_coef_all_diff)
      deallocate(pod_coef_all_ann)
      deallocate(pod_coef_all_diff_obj)
     end subroutine ann_bp_main

     subroutine readtestData() 
       integer :: k,j
       integer :: i, timestep
        
       open(1,file='in')
       ! do i=1, 820
         ! do j=1,2
            read(1,*)((d_in(i,j),j=1,2),i=1,820) 
        !  enddo
       ! enddo
       close(1)
         
  
        open(1,file='out')
         read(1,*)(d_out(i,1),i=1,820) 
        close(1)
          do i=1,820
            print *, 'out' ,d_in(i,1), d_in(i,2), d_out(i,1) 
          enddo
        !stop 153
     end subroutine readtestData

 
     subroutine readData(l)
       integer :: k,j
       !integer, intent(in) :: i, timestep
       double precision :: d_in_tmp(200), d_in_diff(200)
       integer, intent(in)  ::l
       DO k=1, sample_num
          d_in(k,1)=pod_coef_all_ann(k,l)
       ENDDO 
         
        DO k=1, sample_num
         d_in(k,2)=pod_coef_all_diff(k,l)
         d_out(k,1)=pod_coef_all_diff_obj(k,l) !test output
        ENDDO
         

     end subroutine readData
    
     subroutine readtestData2(l)
       integer :: k,j
       integer, intent(in)  ::l
       DO k=1, sample_num
          d_in(k,1)=pod_coef_all_ann(k,l)
       ENDDO 
         
        DO k=1, sample_num
         d_in(k,2)=pod_coef_all_ann(k+1,l)
         d_out(k,1)=pod_coef_all_ann(k+2,l) !test output
        ENDDO
         
    ! print *,'d_ind_out', d_in,d_out
     end subroutine readtestData2


     
     subroutine initBPNetwork() ! normalisation, between [0,1]
       integer :: i,j
       real::x
        Do i=1,in_num    ! cannot use minval because some of the vector is null.          
          Minin(i)=d_in(1,i)
          Maxin(i)=d_in(1,i)      
         do j=1,sample_num
           if(Minin(i)>d_in(j,i)) then
              Minin(i)=d_in(j,i)          
           endif
           if(Maxin(i)<d_in(j,i)) then
              Maxin(i)=d_in(j,i)
           endif
         enddo
        ENDDO
 
       Do i=1,out_num
            Minout(i)=d_out(1,i)
            Maxout(i)=d_out(1,i)
         do j=1,sample_num
            if(Minout(i)>d_out(j,i)) then
              Minout(i)=d_out(j,i)
           endif
           if(Maxout(i)<d_out(j,i)) then
              Maxout(i)=d_out(j,i)
           endif
       enddo 
         ! print *, 'Minin(1),Maxin(1),Minin(2),Maxin(2)',Minin(1),Maxin(1),Minin(2),Maxin(2)     
      ENDDO

      open(1,file='maxmin')
          write(1,*) Minin(1),Maxin(1),Minin(2),Maxin(2),Maxout(1),Minout(1)
      close(1)
       
      DO i=1, in_num
        do j=1,sample_num
          d_in(j,i)=(d_in(j,i)-Minin(i)+1)/(Maxin(i)-Minin(i)+1)
      enddo
      ENDDO

      DO i=1,out_num
        do j=1,sample_num
          d_out(j,i)=(d_out(j,i)-Minout(i)+1)/(Maxout(i)-Minout(i)+1);
      enddo
      ENDDO
   
          do i=1,820
          !   print *, 'out' ,d_in(i,1), d_in(i,2), d_out(i,1) 
           enddo
     ! stop 208
     call random_seed ()
     !  print *, 'Neuroninnit', Neuron
       DO i=1,Neuron
        do j=1,in_num             
          call random_number (x)          
          w(i,j)=x*2-1 !rand()*2.0/RAND_MAX-1;
	  dw(i,j)=0;
      enddo
      ENDDO
   call random_seed () 
     DO i=1,Neuron
        do j=1,out_num
             
           call random_number (x)
           v(j,i)=x*2-1!rand()*2.0/RAND_MAX-1;
	   dv(j,i)=0;
      enddo
      ENDDO
    ! print *, 'w,v', w,v
     end subroutine initBPNetwork

     subroutine computO(var)
      integer :: var,i,j
      real ::  s,y  !s :sum
       DO i=1,Neuron
         s=0
        do j=1,in_num
          s=s+w(i,j)*d_in(var,j)
        enddo
	  o(i)=1./(1.+exp(-s))
        ENDDO     

      DO i=1,out_num
         s=0
        do j=1,Neuron
           s=s+v(i,j)*o(j)
         enddo
	   OutputData(i)=s
        ENDDO
     end subroutine computO


   subroutine backUpdate(var)
     integer :: var,i,j
     real :: t

      DO i=1,Neuron
         t=0
        do j=1,out_num
         t=t+(OutputData(j)-d_out(var,j))*v(j,i)
         dv(j,i)=A*dv(j,i)+B*(OutputData(j)-d_out(var,j))*o(i)
	 v(j,i)=v(j,i)-dv(j,i)  
        enddo

        do j=1,in_num
         dw(i,j)=aa*dw(i,j)+bb*t*o(i)*(1-o(i))*d_in(var,j)
	 w(i,j)=w(i,j)-dw(i,j)
        enddo
       
      ENDDO   

   end subroutine backUpdate

  subroutine  trainNetwork()
     integer :: var,i,j,training_num
     real :: e
     training_num=0
 	
 
 
 120 e=0
     DO i=1, sample_num
          call computO(i)
          do j=1,out_num
           e=e+abs((OutputData(j)-d_out(i,j))/d_out(i,j))
          enddo
          call backUpdate(i) 
      ENDDO
     ! print *, 'sample_num , whole_error and traning number' ,sample_num, e/sample_num  ,training_num   
     training_num=training_num+1
 
    if (e/sample_num>0.001 .and.training_num<Train_num) goto 120
    print *,'training is over the training number is ', training_num
  end subroutine

  function bpresult(var1,var2) result (s)
   real :: var1,var2,s,y
   integer :: i,j
   var1=(var1-Minin(1)+1)/(Maxin(1)-Minin(1)+1);
   var2=(var2-Minin(2)+1)/(Maxin(2)-Minin(2)+1);
   do i=1,Neuron 
      s=0
      s=w(i,1)*var1+w(i,2)*var2 
      o(i)=1/(1+exp(-s)) 
   enddo
   s=0
    do j=1,Neuron
       s=s+v(1,j)*o(j)
    enddo 
     s=s*(Maxout(1)-Minout(1)+1)+Minout(1)-1 
  end function bpresult

 subroutine bpresultok(var1,var2,s)
   real::var1,var2
   real, intent(inout) :: s 
   integer :: i,j 
  ! call initBPNetwork()
  ! Neuron=45
  ! print *, 'var1,var2',var1,var2
   var1=(var1-Minin(1)+1)/(Maxin(1)-Minin(1)+1);
   var2=(var2-Minin(2)+1)/(Maxin(2)-Minin(2)+1);
    print *, 'var1,var2,bpresult',var1,var2,Neuron
     print *, 'Minin(1),Maxin(1),Minin(2),Maxin(2)',Minin(1),Maxin(1),Minin(2),Maxin(2)
    !stop 286
    ! open(1,file='neuron.txt')
       !    read(1,*)((w(i,j),i=1,45),j=1,2) 
     !      read(1,*)(v(1,i),i=1,45) 
   !  close(1)

   do i=1,Neuron 
      s=0
      s=w(i,1)*var1+w(i,2)*var2 
     !   print *, 'we',w(i,1),w(i,2) 
      o(i)=1/(1+exp(-s)) 
    !    print *, 'out',o(i)
     !stop 326
   enddo
   s=0
    do j=1,Neuron
       s=s+v(1,j)*o(j)
    enddo 
     s=s*(Maxout(1)-Minout(1)+1)+Minout(1)-1 
     print *, 's' , s
  end subroutine bpresultok

   subroutine writeweight(l)
     integer, intent(in) :: l
     integer :: i,j
     if(l.eq. 1) then 
      open(1,file='weight_w')
          write(1,*)((w(i,j),j=1,2),i=1,45) 
         ! write(1,*)(o(i),i=1,45) 
       close(1)
      else 
         open(1,file='weight_w', position='append',ACTION='WRITE')
          write(1,*)((w(i,j),j=1,2),i=1,45) 
         ! write(1,*)(o(i),i=1,45) 
       close(1)
     endif 
      if(l.eq. 1) then 
      open(2,file='weight_o')
          !write(1,*)((w(i,j),j=1,2),i=1,45) 
          write(2,*)(o(i),i=1,45) 
       close(2)
      else 
         open(2,file='weight_o', position='append',ACTION='WRITE')
          !write(1,*)((w(i,j),j=1,2),i=1,45) 
          write(2,*)(o(i),i=1,45) 
       close(2)
     endif 
     
       if(l.eq. 1) then 
      open(1,file='weight_v')
          write(1,*)((v(i,j),j=1,45),i=1,1) 
       close(1)
      else 
         open(1,file='weight_v', position='append',ACTION='WRITE')
          write(1,*)((v(i,j),j=1,45),i=1,1) 
       close(1)
     endif 
   end subroutine writeweight

   subroutine readweight(l)
     integer, intent(in) :: l
     integer :: i,j,k
     real :: w_all(360,45,2),o_all(360,45),v_all(360,1,45)
      open(1,file='maxmin')
          read(1,*) Minin(1),Maxin(1),Minin(2),Maxin(2),Maxout(1),Minout(1)
      close(1)
     
      do k=1, coef_num
      open(1,file='weight_w')
          read(1,*)((w_all(k,i,j),j=1,2),i=1,45)           
      close(1)       
      enddo
 
      do k=1, coef_num
      open(2,file='weight_o')
          read(2,*)(o_all(k,i),i=1,45)           
       close(2)       
      enddo
      
        do k=1, coef_num
      open(1,file='weight_v')
          read(1,*)((v_all(k,i,j),j=1,45),i=1,1)           
      close(1)       
      enddo
      w(:,:)=w_all(l,:,:)
      o(:)=o_all(l,:)
      v(:,:)=v_all(l,:,:)
   end subroutine readweight

   end module ann_bp
