module ggint
  use precision, only : prec
  use file_OUT, only : LOUT
  implicit none
  integer, parameter :: newhn=60
  integer, parameter :: ggn=400,hn=200,ghn=100,tho=1000,RECB=36,RECS=20!,RECB=32,RECS=16!
  integer, parameter :: ROIsum=50
  real(prec) :: hermiteh_gh(0:hn,ggn),ghabsciss(ghn),ghweights(ghn)!,ggabsciss(:),ggweights(:),hermiteh_gg(4,0:hn,ggn)
  real(prec) :: dzejmu(4,-5:tho),norm(0:tho),aleph=1._prec
  real(prec) :: factorials(0:tho),binomials(0:tho,0:tho),coo(0:tho,0:tho)
  real(prec), parameter :: sqrtpi&
       &=1.77245385090551602729816748334114518_prec
  real(prec), parameter :: sqrttwo&
       &=1.41421356237309504880168872420969807857_prec

contains

  subroutine density_matrix(mscf,n,eorb,orbvec,occdens,occdensE)
    implicit none
    integer    :: i,j,mscf,n,al,ii
    real(prec) :: eorb(:),orbvec(:,:),occdens(:,:),occdensE(:,:)
    real(prec) :: val,valE,t1,t0
    occdens = 0._prec
    occdensE = 0._prec
    do j = 1, mscf
       do i = 1, j
          val  = sum(orbvec(i,1:n)*orbvec(j,1:n))
          valE = sum(eorb(1:n)*orbvec(i,1:n)*orbvec(j,1:n))
          occdens(i,j) = val
          occdens(j,i) = val
          occdensE(i,j) = valE
          occdensE(j,i) = valE
       end do
    end do
  end subroutine density_matrix

  subroutine make_vecP(m,mscf,n,vecp,int6f12table,orbvec,occdens)
    implicit none
    integer    :: m,mscf,n,i,j,B,A,al,ij,k,l,p,r
    real(prec) :: temp,temp2,val1,val2
    real(prec) :: vecp(:,:,:),orbvec(:,:),occdens(:,:)
    real(prec) :: int6f12table(0:,0:,0:,0:,0:,0:)
    vecp = 0._prec
    do j = 0, m-1
       do i = 0, j
          ij = j*(j+1)/2+i+1
          do B = 1, n
             do A = 1, B
                temp2 = 0._prec
                do l = 0, mscf-1
                   do k = 0, mscf-1
                      temp = 0._prec
                      do r = 0, mscf-1, 2
                         do p = 0, mscf-1, 2
                            val1 = occdens(r+1,p+1)*(int6(r,k,l,j,p,i,int6f12table)+int6(r,k,l,i,p,j,int6f12table))
                            temp = temp + val1
                         end do
                      end do
                      do r = 1, mscf-1, 2
                         do p = 1, mscf-1, 2
                            val1 = occdens(r+1,p+1)*(int6(r,k,l,j,p,i,int6f12table)+int6(r,k,l,i,p,j,int6f12table))
                            temp = temp + val1
                         end do
                      end do
                      temp2 = temp2 + orbvec(k+1,A)*orbvec(l+1,B)*temp
                   end do
                end do
                vecp(ij,A,B) = vecp(ij,A,B) + temp2
                if (A.NE.B) vecp(ij,B,A) = vecp(ij,B,A) + temp2
             end do
          end do
       end do
    end do
  end subroutine make_vecP

  real(prec) function int6(i,j,k,l,a,b,inttable)
    implicit none
    integer    :: i,j,k,l,ii,jj,kk,ll,a,b,aa,bb
    real(prec) :: inttable(0:,0:,0:,0:,0:,0:)
    call sort4(i,j,k,l,ii,jj,kk,ll)
    call sort2(a,b,aa,bb)
    int6 = inttable(ii,jj,kk,ll,aa,bb)
  end function int6

  subroutine read_f12(intf12)
    implicit none
    integer    :: ii,i,j,k,l,number
    real(prec) :: temp
    real(prec) :: intf12(0:,0:,0:,0:)
    open(11,file='dat/file_f12.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    intf12 = 0._prec
    do ii = 1, 779968 !11029041 !for 80
       read(11,rec=ii) number,temp
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       k=iand(ishft(number,-16),255)
       l=ishft(number,-24)
       intf12(i,j,k,l) = temp
       intf12(k,j,i,l) = temp
       intf12(i,l,k,j) = temp
       intf12(k,l,i,j) = temp
       ! print*, i,j,k,l,temp
    end do
    ! stop
    close(11)
  end subroutine read_f12

  subroutine create_ROI(m,mscf,n,matps,matpt,matpes,matpet,intf12,occdens,occdensE)
    implicit none
    integer    :: i,j,k,l,counter,m,n,ij,kl,p,q,r,rs,mscf
    real(prec) :: intf12(0:,0:,0:,0:)
    real(prec) :: temp1,temp2,temp3,temp4,ttemp1,ttemp2,ttemp3,ttemp4
    real(prec) :: temp1E,temp2E,temp3E,temp4E
    real(prec) :: val,valE,val1,val2,val3,val4,occdens(:,:),occdensE(:,:)
    real(prec) :: matps(:,:),matpt(:,:),matpes(:,:),matpet(:,:)
    matps = 0._prec
    matpt = 0._prec
    matpes = 0._prec
    matpet = 0._prec
    do j = 0, m-1
       do i = 0, j
          ij = j*(j+1)/2+i+1
          do l = 0, m-1
             do k = 0, l
                kl = l*(l+1)/2+k+1
                ttemp1 = 0._prec
                ttemp2 = 0._prec
                ttemp3 = 0._prec
                ttemp4 = 0._prec
                temp1E = 0._prec
                temp2E = 0._prec
                temp3E = 0._prec
                temp4E = 0._prec
                do p = 0, mscf-1
                   do r = 0, mscf-1
                      temp1 = 0._prec
                      temp2 = 0._prec
                      temp3 = 0._prec
                      temp4 = 0._prec
                      do q = 0, ROIsum
                         val1 = intf12(q,r,j,i)
                         val2 = intf12(q,p,l,k)
                         val3 = intf12(q,r,i,j)
                         val4 = intf12(q,p,k,l)
                         temp1 = temp1 + val1*val2
                         temp2 = temp2 + val3*val4
                         temp3 = temp3 + val1*val4
                         temp4 = temp4 + val2*val3
                      end do
                      val = occdens(p+1,r+1)
                      valE = occdensE(p+1,r+1)
                      ttemp1 = ttemp1 + val*temp1
                      ttemp2 = ttemp2 + val*temp2
                      ttemp3 = ttemp3 + val*temp3
                      ttemp4 = ttemp4 + val*temp4
                      temp1E = temp1E + valE*temp1
                      temp2E = temp2E + valE*temp2
                      temp3E = temp3E + valE*temp3
                      temp4E = temp4E + valE*temp4
                   end do
                end do
                matps(ij,kl) = ttemp1+ttemp2+ttemp3+ttemp4
                matpt(ij,kl) = ttemp1+ttemp2-ttemp3-ttemp4
                matpes(ij,kl) = temp1E+temp2E+temp3E+temp4E
                matpet(ij,kl) = temp1E+temp2E-temp3E-temp4E
             end do
          end do
       end do
    end do
  end subroutine create_ROI

  ! subroutine make_p1(m,n,matps,matpt,eorb,matpes,matpet,intf12)
  !   implicit none
  !   integer    :: m,n
  !   integer    :: i,j,k,l,al,p,ij,kl
  !   real(prec) :: temp,tempe,beta
  !   real(prec) :: tttemp
  !   real(prec) :: matps(:,:), matpt(:,:), matpes(:,:), matpet(:,:)
  !   real(prec) :: eorb(n),intf12(0:,0:,0:,0:)
  !   beta = 2*aleph+1._prec
  !   matps  = 0._prec
  !   matpt  = 0._prec
  !   matpes = 0._prec
  !   matpet = 0._prec
  !   do l = 0, m-1
  !      do k = 0, l
  !         kl=l*(l+1)/2+k+1
  !         do j = 0, m-1
  !            do i = 0, j
  !               ij=j*(j+1)/2+i+1
  !               temp = 0._prec
  !               tempe = 0._prec
  !               do p = 0, ROIsum
  !                  do al = 0, n-1
  !                     tttemp = intf12(i,j,al,p)*intf12(al,p,k,l)
  !                     temp = temp + tttemp
  !                     tempe = tempe + tttemp*eorb(al+1)
  !                  end do
  !               end do
  !               matps(ij,kl) = matps(ij,kl) + temp
  !               matpt(ij,kl) = matpt(ij,kl) + temp
  !               matpes(ij,kl) = matpes(ij,kl) + tempe
  !               matpet(ij,kl) = matpet(ij,kl) + tempe
  !               temp = 0._prec
  !               tempe = 0._prec
  !               do p = 0, ROIsum
  !                  do al = 0, n-1
  !                     tttemp = intf12(j,i,al,p)*intf12(al,p,l,k)
  !                     temp = temp + tttemp
  !                     tempe = tempe + tttemp*eorb(al+1)
  !                  end do
  !               end do
  !               matps(ij,kl) = matps(ij,kl) + temp
  !               matpt(ij,kl) = matpt(ij,kl) + temp
  !               matpes(ij,kl) = matpes(ij,kl) + tempe
  !               matpet(ij,kl) = matpet(ij,kl) + tempe
  !               temp = 0._prec
  !               tempe = 0._prec
  !               do p = 0, ROIsum
  !                  do al = 0, n-1
  !                     tttemp = intf12(j,i,al,p)*intf12(al,p,k,l)
  !                     temp = temp + tttemp
  !                     tempe = tempe + tttemp*eorb(al+1)
  !                  end do
  !               end do
  !               matps(ij,kl) = matps(ij,kl) + temp
  !               matpt(ij,kl) = matpt(ij,kl) - temp
  !               matpes(ij,kl) = matpes(ij,kl) + tempe
  !               matpet(ij,kl) = matpet(ij,kl) - tempe
  !               temp = 0._prec
  !               tempe = 0._prec
  !               do p = 0, ROIsum
  !                  do al = 0, n-1
  !                     tttemp = intf12(i,j,al,p)*intf12(al,p,l,k)
  !                     temp = temp + tttemp
  !                     tempe = tempe + tttemp*eorb(al+1)
  !                  end do
  !               end do
  !               matps(ij,kl) = matps(ij,kl) + temp
  !               matpt(ij,kl) = matpt(ij,kl) - temp
  !               matpes(ij,kl) = matpes(ij,kl) + tempe
  !               matpet(ij,kl) = matpet(ij,kl) - tempe
  !            end do
  !         end do
  !      end do
  !   end do
  ! end subroutine make_p1

  subroutine make_mats(m,mscf,n,int6f122table,matjs,matjt,matms,matmt,orbvec)
    implicit none
    integer :: m,n,counter,mscf
    integer :: i,j,k,l,p,r,A,B,ij,kl,ii
    real(prec) :: temp,temp1,temp2,temp3,temp4
    real(prec) :: val, val1, val2, val3, val4
    real(prec) :: int6f122table(0:,0:,0:,0:,0:,0:)
    real(prec) :: orbvec(:,:)
    real(prec) :: matjs(:,:,:,:), matjt(:,:,:,:)
    real(prec) :: matms(:,:,:,:), matmt(:,:,:,:)
    matjs=0._prec
    matjt=0._prec
    matms=0._prec
    matmt=0._prec
    counter = 0
    do B = 1, n
       do A = 1, B
          do l = 0, m-1
             do k = 0, l
                kl=l*(l+1)/2+k+1
                do j = 0, m-1
                   do i = 0, j
                      ij=j*(j+1)/2+i+1
                      temp1 = 0._prec
                      temp2 = 0._prec
                      temp3 = 0._prec
                      temp4 = 0._prec
                      do p = 0, mscf-1
                         do r = 0, p
                            if (modulo(i+j+k+l+p+r,2).EQ.0) then
                               temp=orbvec(r+1,A)*orbvec(p+1,B)
                               val1=int6(r,p,i,k,j,l,int6f122table)
                               val2=int6(r,p,j,l,i,k,int6f122table)
                               val3=int6(r,p,i,l,j,k,int6f122table)
                               val4=int6(r,p,j,k,i,l,int6f122table)
                               temp1=temp1+temp*val1
                               temp2=temp2+temp*val2
                               temp3=temp3+temp*val3
                               temp4=temp4+temp*val4
                               if (r.NE.p) then
                                  temp=orbvec(p+1,A)*orbvec(r+1,B)
                                  temp1=temp1+temp*val1
                                  temp2=temp2+temp*val2
                                  temp3=temp3+temp*val3
                                  temp4=temp4+temp*val4
                               end if
                            end if
                         end do
                      end do
                      val=temp1+temp2+temp3+temp4
                      matjs(ij,kl,A,B)=val
                      matjs(ij,kl,B,A)=val
                      val=temp1+temp2-temp3-temp4
                      matjt(ij,kl,A,B)=val
                      matjt(ij,kl,B,A)=val
                      val=temp1-temp2+temp3-temp4
                      matms(ij,kl,A,B)=val
                      matms(ij,kl,B,A)=val
                      val=temp1-temp2-temp3+temp4
                      matmt(ij,kl,A,B)=val
                      matmt(ij,kl,B,A)=val
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine make_mats

  subroutine read_int6f12_file(g,int6f12table,fun)
    implicit none
    integer :: fun
    integer :: number,numberab,ii,i,j,k,l,a,b
    real(prec) :: lolzix,g
    real(prec) :: int6f12table(0:,0:,0:,0:,0:,0:)
    open(11,file='bas/file_int6f12.F',status='unknown',form='unformatted',access='direct',RECL=RECB)
    int6f12table = 0._prec
    do ii = 1, fun
       read(11,rec=ii) number,numberab,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       k=iand(ishft(number,-16),255)
       l=ishft(number,-24)
       a=iand(numberab,255)
       b=iand(ishft(numberab,-8),255)
       lolzix = lolzix*g
       int6f12table(i,j,k,l,a,b) = lolzix
    end do
    close(11)
  end subroutine read_int6f12_file

  subroutine create_int6f12_file(m,intinterm12,caunta)
    implicit none
    integer :: a,b,i,j,k,l
    integer, intent(in) :: m
    integer, intent(out) :: caunta
    real(prec) :: integral
    real(prec) :: intinterm12(0:,0:,0:,0:)
    open(11,file='bas/file_int6f12.F',status='unknown',form='unformatted',access='direct',RECL=RECB)
    caunta = 0
    do l = 0, m-1
       do k = 0, l
          do j = 0, k
             do i = 0, j
                do b = 0, m-1
                   do a = 0, b
                      if ((modulo(i+j+k+l+a+b,2).EQ.0)) then
                         integral=inttildef12(i,j,k,l,a,b,intinterm12)
                         caunta = caunta + 1
                         write(11,rec=caunta) i+ishft(j,8)+ishft(k,16)+ishft(l,24),a+ishft(b,8),integral
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do
    close(11)
  end subroutine create_int6f12_file

  subroutine read_int6f122_file(g,int6f122table,fun)
    implicit none
    integer :: fun
    integer :: number,numberab,ii,i,j,k,l,a,b
    real(prec) :: lolzix,g
    real(prec) :: int6f122table(0:,0:,0:,0:,0:,0:)
    open(11,file='bas/file_int6f122.F',status='unknown',form='unformatted',access='direct',RECL=RECB)
    int6f122table=0._prec
    do ii = 1, fun
       read(11,rec=ii) number,numberab,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       k=iand(ishft(number,-16),255)
       l=ishft(number,-24)
       a=iand(numberab,255)
       b=iand(ishft(numberab,-8),255)
       lolzix = lolzix*g
       int6f122table(i,j,k,l,a,b) = lolzix
    end do
    close(11)
  end subroutine read_int6f122_file

  subroutine create_int6f122_file(m,intinterm22,caunta)
    implicit none
    integer :: a,b,i,j,k,l
    integer, intent(in) :: m
    integer, intent(out) :: caunta
    real(prec) :: integral,beta
    real(prec) :: intinterm22(0:,0:,0:,0:)
    beta = 3._prec*aleph+1._prec
    open(11,file='bas/file_int6f122.F',status='unknown',form='unformatted',access='direct',RECL=RECB)
    caunta = 0
    do l = 0, m-1
       do k = 0, l
          do j = 0, k
             do i = 0, j
                do b = 0, m-1
                   do a = 0, b
                      if ((modulo(i+j+a+b+k+l,2).EQ.0)) then
                         integral=inttildef122(i,j,k,l,a,b,intinterm22)
                         caunta = caunta + 1
                         write(11,rec=caunta) i+ishft(j,8)+ishft(k,16)+ishft(l,24),a+ishft(b,8),integral
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do
    close(11)
  end subroutine create_int6f122_file

  subroutine read_vec0_file(vec0,counter)
    implicit none
    integer :: number,ii,i,j,A,B,jj,counter
    real(prec) :: lolzix
    real(prec) :: vec0(:,:,:)
    open(11,file='bas/file_int2_vec.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    vec0 = 0._prec
    do ii = 1, counter
       read(11,rec=ii) number,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       A=iand(ishft(number,-16),255)
       B=ishft(number,-24)
       jj=j*(j+1)/2+i+1
       vec0(jj,A,B) = lolzix
       vec0(jj,B,A) = lolzix
    end do
    close(11)
  end subroutine read_vec0_file

  subroutine create_vec0_file(m,mscf,n,counter,intf12,orbvec)
    implicit none
    integer :: i,j,k,l,A,B,counter,ii
    integer, intent(in) :: m,n,mscf
    real(prec) :: integral,beta,val
    real(prec) :: intf12(0:,0:,0:,0:),orbvec(:,:)
    beta = 2._prec*aleph+1._prec
    open(11,file='bas/file_int2_vec.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    counter = 0
    do B = 1, n
       do A = 1, B
          do j = 0, m-1
             do i = 0, j
                integral = 0._prec
                do l = 0, mscf-1
                   do k = 0, l
                      if (modulo(i+j+k+l,2).EQ.0) then
                         val = intf12(k,i,l,j)
                         integral = integral + orbvec(k+1,A)*orbvec(l+1,B)*val
                         if (l.NE.k) integral = integral + orbvec(l+1,A)*orbvec(k+1,B)*val
                      end if
                   end do
                end do
                counter = counter + 1
                write(11,rec=counter) i+ishft(j,8)+ishft(A,16)+ishft(B,24),integral
             end do
          end do
       end do
    end do
    close(11)
  end subroutine create_vec0_file

  subroutine read_matH_file(maths,matht,counta)
    implicit none
    integer :: counta
    integer :: number,ii,i,j,k,l,jj,kk
    real(prec) :: lolzix,lolzix2
    real(prec) :: maths(:,:),matht(:,:)
    open(11,file='bas/fileH.F',status='unknown',form='unformatted',access='direct',RECL=RECB)
    maths=0._prec
    matht=0._prec
    do ii = 1, counta
       read(11,rec=ii) number,lolzix,lolzix2
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       k=iand(ishft(number,-16),255)
       l=ishft(number,-24)
       jj=j*(j+1)/2+i+1
       kk=l*(l+1)/2+k+1
       maths(jj,kk) = maths(jj,kk) + lolzix
       matht(jj,kk) = matht(jj,kk) + lolzix2
    end do
    close(11)
  end subroutine read_matH_file

  subroutine create_matH_file(m,counta)
    implicit none
    integer :: i,j,k,l,counta
    integer, intent(in) :: m
    real(prec) :: integral,integral2,beta1,beta2,maths,matht
    open(11,file='bas/fileH.F',status='unknown',form='unformatted',access='direct',RECL=RECB)
    counta = 0
    beta1 = 4._prec*aleph+1._prec
    beta2 = 3._prec*aleph+1._prec
    do l = 0, m-1
       do k = 0, l
          do j = 0, m-1
             do i = 0, j
                if ((modulo(i+j+k+l,2).EQ.0)) then
                   integral=0.5_prec*(intdf122(i,j,k,l,beta1)+(i+j+k+l+2)*intf122(i,j,k,l,beta2))
                   integral2=0.5_prec*(intdf122(i,j,l,k,beta1)+(i+j+l+k+2)*intf122(i,j,l,k,beta2))
                   maths=integral+integral2
                   matht=integral-integral2
                   if (abs(maths).LT.100000*epsilon(0._prec)) then
                      maths = 0._prec
                   endif
                   if (abs(matht).LT.100000*epsilon(0._prec)) then
                      matht = 0._prec
                   endif
                   counta = counta + 1
                   write(11,rec=counta) i+ishft(j,8)+ishft(k,16)+ishft(l,24),maths,matht
                end if
             end do
          end do
       end do
    end do
    ! do l = n, m-1
    !    do k = n, l
    !       do j = n, m-1
    !          do i = n, j
    ! do l = 0, m-1
    !    do k = 0, l
    !       do j = 0, m-1
    !          do i = 0, j
    !             integral = 0._prec
    !             integral2 = 0._prec
    !             if ((i.EQ.k).AND.(j.EQ.l)) then
    !                integral = real(k+l+1,prec)
    !             end if
    !             ! if ((i.EQ.l).AND.(j.EQ.k)) then
    !             !    integral2 = real(k+l+1,prec)
    !             ! end if
    !             maths=integral+integral2
    !             matht=integral-integral2
    !             if ((maths.NE.0._prec).OR.(matht.NE.0._prec)) then
    !                counta = counta + 1
    !                write(11,rec=counta) i+ishft(j,8)+ishft(k,16)+ishft(l,24),maths,matht
    !             end if
    !          end do
    !       end do
    !    end do
    ! end do
    ! do l = 0, n-1
    !    do k = 0, l
    !       do j = n, m-1
    !          do i = n, j
    !             if ((modulo(i+j+k+l,2).EQ.0)) then
    !                integral=(i+j+1)*intf12(i,j,k,l)
    !                integral2=(i+j+1)*intf12(i,j,l,k)
    !                maths=integral+integral2
    !                matht=integral-integral2
    !                if (abs(maths).LT.100000*epsilon(0._prec)) then
    !                   maths = 0._prec
    !                endif
    !                if (abs(matht).LT.100000*epsilon(0._prec)) then
    !                   matht = 0._prec
    !                endif
    !                counta = counta + 1
    !                write(11,rec=counta) i+ishft(j,8)+ishft(k,16)+ishft(l,24),maths,matht
    !             end if
    !          end do
    !       end do
    !    end do
    ! end do
    ! do l = n, m-1
    !    do k = n, l
    !       do j = 0, n-1
    !          do i = 0, j
    !             if ((modulo(i+j+k+l,2).EQ.0)) then
    !                integral=(k+l+1)*intf12(i,j,k,l)
    !                integral2=(k+l+1)*intf12(i,j,l,k)
    !                maths=integral+integral2
    !                matht=integral-integral2
    !                if (abs(maths).LT.100000*epsilon(0._prec)) then
    !                   maths = 0._prec
    !                endif
    !                if (abs(matht).LT.100000*epsilon(0._prec)) then
    !                   matht = 0._prec
    !                endif
    !                counta = counta + 1
    !                write(11,rec=counta) i+ishft(j,8)+ishft(k,16)+ishft(l,24),maths,matht
    !             end if
    !          end do
    !       end do
    !    end do
    ! end do
    close(11)
  end subroutine create_matH_file

  subroutine read_matS_file(matss,matst,counta)
    implicit none
    integer :: counta
    integer :: number,ii,i,j,k,l,jj,kk
    real(prec) :: lolzix,lolzix2
    real(prec) :: matss(:,:),matst(:,:)
    open(11,file='bas/fileS.F',status='unknown',form='unformatted',access='direct',RECL=RECB)
    matss=0._prec
    matst=0._prec
    do ii = 1, counta
       read(11,rec=ii) number,lolzix,lolzix2
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       k=iand(ishft(number,-16),255)
       l=ishft(number,-24)
       jj=j*(j+1)/2+i+1
       kk=l*(l+1)/2+k+1
       matss(jj,kk) = matss(jj,kk) + lolzix
       matst(jj,kk) = matst(jj,kk) + lolzix2
    end do
    close(11)
  end subroutine read_matS_file

  subroutine create_matS_file(m,counta)
    implicit none
    integer :: i,j,k,l,counta
    integer, intent(in) :: m
    real(prec) :: integral,integral2,beta,matss,matst
    open(11,file='bas/fileS.F',status='unknown',form='unformatted',access='direct',RECL=RECB)
    counta = 0
    beta = 3._prec*aleph+1._prec
    do l = 0, m-1
       do k = 0, l
          do j = 0, m-1
             do i = 0, j
                if ((modulo(i+j+k+l,2).EQ.0)) then
                   integral=intf122(i,j,k,l,beta)
                   integral2=intf122(i,j,l,k,beta)
                   matss=integral+integral2
                   matst=integral-integral2
                   if (abs(matss).LT.100000*epsilon(0._prec)) then
                      matss = 0._prec
                   endif
                   if (abs(matst).LT.100000*epsilon(0._prec)) then
                      matst = 0._prec
                   endif
                   counta = counta + 1
                   write(11,rec=counta) i+ishft(j,8)+ishft(k,16)+ishft(l,24),matss,matst
                end if
             end do
          end do
       end do
    end do
    ! do l = n, m-1
    !    do k = n, l
    !       do j = n, m-1
    !          do i = n, j
    !             integral = 0._prec
    !             integral2 = 0._prec
    !             if ((i.EQ.k).AND.(j.EQ.l)) then
    !                integral = 1._prec
    !             end if
    !             ! if ((i.EQ.l).AND.(j.EQ.k)) then
    !             !    integral2 = 1._prec
    !             ! end if
    !             matss=integral+integral2
    !             matst=integral-integral2
    !             if ((matss.NE.0._prec).OR.(matst.NE.0._prec)) then
    !                counta = counta + 1
    !                write(11,rec=counta) i+ishft(j,8)+ishft(k,16)+ishft(l,24),matss,matst
    !             end if
    !          end do
    !       end do
    !    end do
    ! end do
    ! beta = 2._prec*aleph+1._prec
    ! do l = 0, n-1
    !    do k = 0, l
    !       do j = n, m-1
    !          do i = n, j
    !             if ((modulo(i+j+k+l,2).EQ.0)) then
    !                integral=intf12(i,j,k,l)
    !                integral2=intf12(i,j,l,k)
    !                matss=integral+integral2
    !                matst=integral-integral2
    !                if (abs(matss).LT.100000*epsilon(0._prec)) then
    !                   matss = 0._prec
    !                endif
    !                if (abs(matst).LT.100000*epsilon(0._prec)) then
    !                   matst = 0._prec
    !                endif
    !                counta = counta + 1
    !                write(11,rec=counta) i+ishft(j,8)+ishft(k,16)+ishft(l,24),matss,matst
    !             end if
    !          end do
    !       end do
    !    end do
    ! end do
    ! do l = n, m-1
    !    do k = n, l
    !       do j = 0, n-1
    !          do i = 0, j
    !             if ((modulo(i+j+k+l,2).EQ.0)) then
    !                integral=intf12(i,j,k,l)
    !                integral2=intf12(j,i,k,l)
    !                matss=integral+integral2
    !                matst=integral-integral2
    !                if (abs(matss).LT.100000*epsilon(0._prec)) then
    !                   matss = 0._prec
    !                endif
    !                if (abs(matst).LT.100000*epsilon(0._prec)) then
    !                   matst = 0._prec
    !                endif
    !                counta = counta + 1
    !                write(11,rec=counta) i+ishft(j,8)+ishft(k,16)+ishft(l,24),matss,matst
    !             end if
    !          end do
    !       end do
    !    end do
    ! end do
    close(11)
  end subroutine create_matS_file

  ! subroutine read_matL1_file(m,matl1s,matl1t,counta)
  !   implicit none
  !   integer :: m,fun,counta
  !   integer :: number,ii,i,j,k,l,jj,kk
  !   real(prec) :: lolzix
  !   real(prec) :: matl1s(:,:),matl1t(:,:)
  !   open(11,file='bas/fileL1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
  !   matl1s=0._prec
  !   matl1t=0._prec
  !   do ii = 1, counta
  !      read(11,rec=ii) number,lolzix
  !      i=iand(number,255)
  !      j=iand(ishft(number,-8),255)
  !      k=iand(ishft(number,-16),255)
  !      l=ishft(number,-24)
  !      jj=j*(j+1)/2+i+1
  !      kk=l*(l+1)/2+k+1
  !      matl1s(jj,kk) = matl1s(jj,kk) + lolzix
  !   end do
  !   close(11)
  ! end subroutine read_matL1_file

  ! subroutine create_matL1_file(m,n,g,counta)
  !   implicit none
  !   integer :: i,j,k,l,counta
  !   integer, intent(in) :: m,n
  !   real(prec), intent(in) :: g
  !   real(prec) :: matl1s
  !   open(11,file='bas/fileL1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
  !   counta = 0
  !   do l = n, m-1
  !      do k = n, l
  !         do j = n, m-1
  !            do i = n, j
  !               matl1s = 2._prec*g*intgh(i,j,k,l,2._prec)
  !               counta = counta + 1
  !               write(11,rec=counta) i+ishft(j,8)+ishft(k,16)+ishft(l,24),matl1s
  !            end do
  !         end do
  !      end do
  !   end do
  !   close(11)
  ! end subroutine create_matL1_file

  subroutine read_SCF_file(m,twe,gie)
    implicit none
    integer :: m,fun
    integer :: number,ii,i,j,k,l,jj,kk
    real(prec) :: lolzix,gie
    real(prec) :: twe(:,:)
    open(10,file='bas/file2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    fun=int(1._prec/96._prec*(81._prec+15._prec*(-1._prec)**(m-1._prec)+2._prec*(65._prec+3._prec*&
         (-1._prec)**(m-1._prec))*(m-1._prec)+76._prec*(m-1._prec)**2._prec+20._prec*(m-1._prec)**3._prec&
         +2._prec*(m-1._prec)**4._prec))
    do ii = 1, fun
       read(10,rec=ii) number,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       k=iand(ishft(number,-16),255)
       l=ishft(number,-24)
       lolzix = lolzix * gie
       jj=i*m+j+1
       kk=k*m+l+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=j*m+i+1
       kk=k*m+l+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=i*m+j+1
       kk=l*m+k+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=j*m+i+1
       kk=l*m+k+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=i*m+k+1
       kk=j*m+l+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=k*m+i+1
       kk=j*m+l+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=i*m+k+1
       kk=l*m+j+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=k*m+i+1
       kk=l*m+j+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=i*m+l+1
       kk=k*m+j+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=i*m+l+1
       kk=j*m+k+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=l*m+i+1
       kk=k*m+j+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
       jj=l*m+i+1
       kk=j*m+k+1
       twe(jj,kk) = lolzix
       twe(kk,jj) = lolzix
    end do
    close(10)
  end subroutine read_SCF_file

  subroutine create_SCF_file(m)
    implicit none
    integer :: i,j,k,l,counter
    integer, intent(in) :: m
    open(10,file='bas/file2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    counter = 0
    do l = 0, m-1
       do k = 0, l
          do j = 0, k
             do i = 0, j
                if ((modulo(modulo(i,2)+modulo(j,2)+modulo(k,2)+modulo(l,2),2).EQ.0)) then
                   counter = counter + 1
                   write(10,rec=counter) i+ishft(j,8)+ishft(k,16)+ishft(l,24),intgh(i,j,k,l,2._prec)
                end if
             end do
          end do
       end do
    end do
    close(10)
  end subroutine create_SCF_file

  real(prec) function coeff(s,i,k)
    implicit none
    integer :: i,k,s
    coeff = factorials(s)*2._prec**real(s,prec)*binomials(i,s)*binomials(k,s)
  end function coeff

  ! subroutine import_gg
  !   implicit none
  !   integer            :: i
  !   open(10,file='dat/ggweights.dat')
  !   ggweights=0._prec
  !   do i = 1, ggn
  !      read(10,*) ggweights(i)
  !   end do
  !   close(10)

  !   open(10,file='dat/ggabsciss.dat')
  !   ggabsciss=0._prec
  !   do i = 1, ggn
  !      read(10,*) ggabsciss(i)
  !   end do
  !   close(10)
  ! end subroutine import_gg

  subroutine import_gh
    implicit none
    integer            :: i
    open(10,file='dat/ghabsciss.dat')
    ghabsciss=0._prec
    do i = 1, ghn
       read(10,*) ghabsciss(i)
    end do
    close(10)

    open(10,file='dat/ghweights.dat')
    ghweights=0._prec
    do i = 1, ghn
       read(10,*) ghweights(i)
    end do
    close(10)
  end subroutine import_gh

  ! subroutine gener_hermiteh_gg(beta)
  !   implicit none
  !   integer         :: i,j,b
  !   real(prec)      :: beta
  !   b = bet(beta)
  !   hermiteh_gg(b,:,:)=0._prec
  !   do i = 1, ggn
  !      hermiteh_gg(b,0,i)=1._prec
  !      hermiteh_gg(b,1,i)=2._prec*ggabsciss(i)/sqrt(beta)
  !   end do
  !   do j = 2, hn
  !      do i = 1, ggn
  !         hermiteh_gg(b,j,i) = 2._prec*ggabsciss(i)/sqrt(beta)*hermiteh_gg(b,j-1,i)-2._prec*(j-1)*hermiteh_gg(b,j-2,i)
  !      end do
  !   end do
  ! end subroutine gener_hermiteh_gg

  subroutine gener_hermiteh_gh(beta)
    implicit none
    integer         :: i,j
    real(prec)      :: beta
    hermiteh_gh=0._prec
    do i = 1, ghn
       hermiteh_gh(0,i)=1._prec
       hermiteh_gh(1,i)=2._prec*ghabsciss(i)/sqrt(beta)
    end do
    do j=2,hn
       do i = 1, ghn
          hermiteh_gh(j,i) = 2._prec*ghabsciss(i)/sqrt(beta)*hermiteh_gh(j-1,i)-2._prec*(j-1)*hermiteh_gh(j-2,i)
       end do
    end do
  end subroutine gener_hermiteh_gh

  function sortt(vec)
    implicit none
    integer               :: n,j,temp,bubble
    integer, dimension(:) :: vec
    integer, dimension(4) :: sortt
    n = 4
    do while (n.GT.1)
       bubble = 0
       do j = 1, (n-1)
          if (vec(j)>vec(j+1)) then
             temp = vec(j)
             vec(j) = vec(j+1)
             vec(j+1) = temp
             bubble = j
          end if
       end do
       n = bubble
    end do
    sortt = vec
  end function sortt

  subroutine sort2(a,b,aa,bb)
    implicit none
    integer :: a,b,aa,bb
    if (a.GT.b) then
       aa = b
       bb = a
    else
       aa = a
       bb = b
    end if
  end subroutine sort2

  subroutine sort4(a,b,c,d,aa,bb,cc,dd)
    implicit none
    integer, intent(inout) :: a,b,c,d,aa,bb,cc,dd
    integer                :: n,j,temp,bubble
    integer                :: vec(4)
    n = 4
    vec(1) = a
    vec(2) = b
    vec(3) = c
    vec(4) = d
    do while (n.GT.1)
       bubble = 0
       do j = 1, (n-1)
          if (vec(j)>vec(j+1)) then
             temp = vec(j)
             vec(j) = vec(j+1)
             vec(j+1) = temp
             bubble = j
          end if
       end do
       n = bubble
    end do
    aa = vec(1)
    bb = vec(2)
    cc = vec(3)
    dd = vec(4)
  end subroutine sort4

  subroutine make_norm()
    implicit none
    integer :: m
    ! norm(:) = 0._prec
    do m = 0, tho
       norm(m) = fnorm(m)
    end do
  end subroutine make_norm

  real(prec) function fnorm(m)
    implicit none
    integer :: m
    fnorm = 1._prec/sqrt(2._prec**m)*1._prec/sqrt(factorials(m)*sqrtpi)
  end function fnorm

  real(prec) function abn(m,n,pm)
    implicit none
    integer    :: m,n,i
    ! real(prec) :: abn
    logical    :: pm
    abn = 1._prec
    if (n.LE.m) then
       if (pm) then
          do i = 1, n
             abn = abn * real(m+i,prec)
          end do
          abn = sqrt(abn) * sqrt(2._prec**real(n,prec))
       else
          do i = 1, n
             abn = abn / real(m-i+1,prec)
          end do
          abn = sqrt(abn) / sqrt(2._prec**real(n,prec))
       end if
    end if
  end function abn

  subroutine concatenate(s_in,f_in,s_out)
    implicit none
    character(*) :: s_in,s_out
    real(prec)   :: f_in
    character(4) :: stmp
    write(stmp,'(f4.2)') f_in
    s_out = trim(s_in)//stmp
    print*, s_out
  end subroutine concatenate

  subroutine import_coo
    implicit none
    integer :: iter,i,n
    open(10,file='dat/coo.dat')
    coo=0._prec
    do iter = 1, 22951
       read(10,*) i,n,coo(i,n)
    end do
    close(10)
  end subroutine import_coo

  subroutine import_factorials
    implicit none
    integer            :: i
    open(10,file='dat/factorials.dat')
    factorials=0._prec
    do i = 0, tho
       read(10,*) factorials(i)
    end do
    close(10)
  end subroutine import_factorials

  subroutine import_binomials
    implicit none
    integer :: i,j
    open(10,file='dat/binomials.dat')
    binomials=0._prec
    do i = 0, tho
       do j = 0, i
          read(10,*) binomials(i,j)
       end do
    end do
    close(10)
  end subroutine import_binomials

  real(prec) function intgh(mm,oo,pp,rr,beta)
    implicit none
    integer    :: m,o,p,r,i
    integer    :: mm,oo,pp,rr,vec(4),vect(4)
    real(prec) :: beta
    intgh = 0._prec
    if (modulo(modulo(mm,2)+modulo(oo,2)+modulo(pp,2)+modulo(rr,2),2).EQ.0) then
       vect(1) = mm
       vect(2) = oo
       vect(3) = pp
       vect(4) = rr
       vec = sortt(vect)
       m = vec(1)
       o = vec(2)
       p = vec(3)
       r = vec(4)
       do i = 1, ghn
          intgh = intgh + hermiteh_gh(m,i)*hermiteh_gh(o,i)*hermiteh_gh(p,i)*hermiteh_gh(r,i) * ghweights(i)
       end do
       intgh = intgh*norm(m)*norm(o)*norm(p)*norm(r)/sqrt(beta)
    end if
  end function intgh

  real(prec) function preintf12(mm,oo,pp,rr,beta)
    implicit none
    integer    :: m,o,p,r,s,t
    integer    :: mm,oo,pp,rr
    real(prec) :: tempf12,beta
    preintf12 = 0._prec
    if (modulo(m+o+p+r,2).EQ.0) then
       if (pp.lt.mm) then
          m=pp
          p=mm
       else
          m=mm
          p=pp
       end if
       if (rr.lt.oo) then
          r=oo
          o=rr
       else
          r=rr
          o=oo
       end if
       do s = 0, m
          tempf12 = 0._prec
          do t = 0, o
             tempf12 = tempf12 + coeff(t,o,r) * imunuf12norm(m+p-2*s,o+r-2*t,beta)/norm(o+r-2*t)
          end do
          preintf12 = preintf12 + coeff(s,m,p) * tempf12/norm(m+p-2*s)
       end do
       preintf12 = norm(mm)*norm(oo)*norm(pp)*norm(rr)*preintf12
    end if
  end function preintf12

  real(prec) function intf12chem(mm,oo,pp,rr,beta)
    implicit none
    integer, intent(in)    :: mm,oo,pp,rr
    real(prec)             :: beta
    intf12chem = preintf12(mm,pp,oo,rr,beta)
  end function intf12chem

  real(prec) function intf122(mm,oo,pp,rr,beta)
    implicit none
    integer    :: m,o,p,r,s,t
    integer    :: mm,oo,pp,rr
    real(prec) :: tempf12,beta
    intf122 = 0._prec
    if (modulo(modulo(mm,2)+modulo(oo,2)+modulo(pp,2)+modulo(rr,2),2).EQ.0) then
       if (pp.lt.mm) then
          m=pp
          p=mm
       else
          m=mm
          p=pp
       end if
       if (rr.lt.oo) then
          r=oo
          o=rr
       else
          r=rr
          o=oo
       end if
       do s = 0, m
          tempf12 = 0._prec
          do t = 0, o
             tempf12 = tempf12 + coeff(t,o,r) * imunuf122norm(m+p-2*s,o+r-2*t,beta)/norm(o+r-2*t)
          end do
          intf122 = intf122 + coeff(s,m,p) * tempf12/norm(m+p-2*s)
       end do
       intf122 = norm(mm)*norm(oo)*norm(pp)*norm(rr)*intf122
    end if
  end function intf122

  real(prec) function intf122chem(mm,oo,pp,rr,beta)
    implicit none
    integer, intent(in)    :: mm,oo,pp,rr
    real(prec)             :: beta
    intf122chem = intf122(mm,pp,oo,rr,beta)
  end function intf122chem

  real(prec) function intdf122(mm,oo,pp,rr,beta) !check this mofo -- 1,1,2,2 works, but 1,2,1,2 also???
    implicit none
    integer    :: m,o,p,r,s,t
    integer    :: mm,oo,pp,rr
    real(prec) :: tempf12,beta
    intdf122 = 0._prec
    if (modulo(mm+oo+pp+rr,2).EQ.0) then
       if (pp.lt.mm) then
          m=pp
          p=mm
       else
          m=mm
          p=pp
       end if
       if (rr.lt.oo) then
          r=oo
          o=rr
       else
          r=rr
          o=oo
       end if
       do s = 0, m
          tempf12 = 0._prec
          do t = 0, o
             tempf12 = tempf12 + coeff(t,o,r) * imunudf122norm(m+p-2*s,o+r-2*t,beta)/norm(o+r-2*t)
          end do
          intdf122 = intdf122 + coeff(s,m,p) * tempf12/norm(m+p-2*s)
       end do
       intdf122 = norm(mm)*norm(oo)*norm(pp)*norm(rr)*intdf122
    end if
  end function intdf122

  real(prec) function intdf122chem(mm,oo,pp,rr,beta)
    implicit none
    integer, intent(in)    :: mm,oo,pp,rr
    real(prec), intent(in) :: beta
    ! real(prec)             :: intdf122chem
    intdf122chem = intdf122(mm,pp,oo,rr,beta)
  end function intdf122chem

  real(prec) function oldimunudf122(m,n,beta)
    integer    :: m,n,b
    real(prec) :: beta
    oldimunudf122 = 0._prec
    if (mod(m+n,2).NE.1) then
       b = bet(beta)
       oldimunudf122 = 2._prec*sqrtpi*(-1._prec)**n*2._prec**(-(m+n)/2._prec)*(&
            &dzejmu(b,m+n)&
            &-(real(b,prec)-1._prec)/2._prec*( dzejmu(b,m+n+2)+&
            &real(4._prec*(m+n)+2,prec)*dzejmu(b,m+n)+&
            &real(4._prec*(m+n)*(m+n-1),prec)*dzejmu(b,m+n-2) )&
            &+((real(b,prec)-1._prec)/4._prec)**2*&
            &( dzejmu(b,m+n+4)+&
            & real(8._prec*(m+n)+12._prec,prec)*dzejmu(b,m+n+2)+&
            & real(24._prec*(m+n)**2+24._prec*(m+n)+12._prec,prec)*dzejmu(b,m+n)+&
            & real(32._prec*(m+n)**3-48._prec*(m+n)**2+16._prec*(m+n),prec)*dzejmu(b,m+n-2)+&
            & real(16._prec*(m+n)*(m+n-1)*(m+n-2)*(m+n-3),prec)*dzejmu(b,m+n-4)&
            &))
    end if
  end function oldimunudf122

  real(prec) function imunudf122(m,n,beta)
    integer    :: m,n,b
    real(prec) :: beta
    imunudf122 = 0._prec
    if (mod(m+n,2).NE.1) then
       b = bet(beta)
       imunudf122 = 2._prec*sqrtpi*(-1._prec)**n*2._prec**(-(m+n)/2._prec)*(dzejmu(b,m+n)-&
            & 2._prec*(beta-1._prec)/4._prec*(dzejmu(b,m+n+2)+&
            & real(4._prec*(m+n)+2,prec)*dzejmu(b,m+n)+&
            & real(4._prec*(m+n)**2-4._prec*(m+n),prec)*dzejmu(b,m+n-2))+&
            & ((beta-1._prec)/4._prec)**2*(dzejmu(b,m+n+4)+&
            & real(8._prec*(m+n)+12._prec,prec)*dzejmu(b,m+n+2)+&
            & real(24._prec*(m+n)**2+24._prec*(m+n)+12._prec,prec)*dzejmu(b,m+n)+&
            & real(32._prec*(m+n)**3-48._prec*(m+n)**2+16._prec*(m+n),prec)*dzejmu(b,m+n-2)+&
            & real(16._prec*(m+n)*(m+n-1)*(m+n-2)*(m+n-3),prec)*dzejmu(b,m+n-4)))
    end if
  end function imunudf122

  real(prec) function imunudf122norm(m,n,beta)
    integer    :: m,n,b
    real(prec) :: beta,mn
    !some super bogus code...
    if (m.EQ.0.and.n.eq.0) then
       imunudf122norm = 0.5390125124744174978149239120563448922870409268698824918081870965
    else if (m.EQ.2.and.n.eq.0.or.m.EQ.0.and.n.eq.2) then
       imunudf122norm = -0.1345197891935502887834707784128512214525037706050587983590004814
    else if (m.EQ.1.and.n.eq.1) then
       imunudf122norm = 0.1902397102850885286405613807257687855130732683070173500499483870
    else
       !/bogus
       mn = real(m+n,prec)
       imunudf122norm = 0._prec
       if (mod(m+n,2).NE.1) then
          b = bet(beta)
          imunudf122norm = norm(n)*norm(m)/norm(m+n)*&
               & 2._prec*sqrtpi*(-1._prec)**n*2._prec**(-(m+n)/2._prec)*&
               &(dzejmu(b,m+n)-&
               & 2._prec*(beta-1._prec)/4._prec*(dzejmu(b,m+n+2)*abn(m+n,2,.TRUE.)+&
               & (4._prec*(m+n)+2._prec)*dzejmu(b,m+n)+&
               & (4._prec*mn**2-4._prec*mn)*dzejmu(b,m+n-2)*abn(m+n,2,.FALSE.))+&
               & ((beta-1._prec)/4._prec)**2*(dzejmu(b,m+n+4)*abn(m+n,4,.TRUE.)+&
               & (8._prec*mn+12._prec)*dzejmu(b,m+n+2)*abn(m+n,2,.TRUE.)+&
               & (24._prec*mn**2+24._prec*mn+12._prec)*dzejmu(b,m+n)+&
               & (32._prec*mn**3-48._prec*mn**2+16._prec*mn)*dzejmu(b,m+n-2)*abn(m+n,2,.FALSE.)+&
               & (16._prec*mn*(mn-1._prec)*(mn-2._prec)*(mn-3._prec))*dzejmu(b,m+n-4)*abn(m+n,4,.FALSE.)))
       end if
    end if
  end function imunudf122norm

  real(prec) function imunuf122norm(mu,nu,beta)
    implicit none
    integer    :: mu,nu,b
    real(prec) :: beta
    imunuf122norm = 0._prec
    if (mod(mu+nu,2).NE.1) then
       b = bet(beta)
       imunuf122norm = norm(nu)*norm(mu)/norm(mu+nu)*4._prec*sqrtpi*(-1._prec)**nu*2._prec**(-(mu+nu)/2._prec)&
            &*(0.25_prec*dzejmu(b,mu+nu+2)*sqrt(4._prec*real(mu+nu+2)*real(mu+nu+1))&
            &+(real(mu+nu,prec)+0.5_prec)*dzejmu(b,mu+nu)&
            &+sqrt(real(mu+nu,prec)*real(mu+nu-1,prec))*dzejmu(b,mu+nu-2)/2._prec)
    end if
  end function imunuf122norm

  real(prec) function imunuf12norm(mu,nu,beta)
    implicit none
    integer    :: mu,nu,b
    real(prec) :: beta
    imunuf12norm = 0._prec
    if (mod(mu+nu,2).NE.1) then
       b = bet(beta)
       imunuf12norm = norm(nu)*norm(mu)/norm(mu+nu)*sqrttwo*sqrtpi*(-1._prec)**nu*2._prec**(-real(mu+nu,prec)/2._prec)&
            &*(sqrt(2._prec*real(mu+nu+1,prec))*dzejmu(b,mu+nu+1)+sqrt(2._prec*(real(mu+nu,prec)))*dzejmu(b,mu+nu-1))
    end if
  end function imunuf12norm

  real(prec) function imunutildef12norm(mu,nu,beta)
    implicit none
    integer    :: mu,nu,b
    real(prec) :: beta
    imunutildef12norm = 0._prec
    if (mod(mu+nu,2).NE.1) then
       b = bet(beta)
       imunutildef12norm = norm(nu)*norm(mu)/norm(mu+nu)*&
            &sqrtpi*2._prec**((real(nu,prec)-1._prec)/2._prec)*3._prec**(-(real(mu+nu,prec)-1._prec)/2._prec)*&
            &(-1._prec)**(real(nu,prec))*&
            &(sqrt(real(mu+nu+1,prec))*dzejmu(b,mu+nu+1)+sqrt(real(mu+nu,prec))*dzejmu(b,mu+nu-1))
    end if
  end function imunutildef12norm

  real(prec) function preintildef12(mu,nu,beta)
    implicit none
    integer    :: mu,nu,i
    real(prec) :: beta
    preintildef12 = 0._prec
    do i = 0, mu/2
       preintildef12 = preintildef12 + coo(i,mu)*imunutildef12norm(mu-2*i,nu,beta)/norm(mu-2*i)
    end do
    preintildef12 = preintildef12*norm(mu)
  end function preintildef12

  subroutine create_intildef12_file()
    implicit none
    integer    :: counter,i,j,nm,m
    real(prec) :: beta
    m=newhn
    open(10,file='bas/file_t1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    counter = 0
    beta = 1.5_prec*aleph+1._prec
    do j = 0, nm
       do i = 0, nm
          if ((modulo(i+j,2).EQ.0)) then
             counter = counter + 1
             write(10,rec=counter) i+ishft(j,8),preintildef12(i,j,beta)
          end if
       end do
    end do
    close(10)
  end subroutine create_intildef12_file

  subroutine read_intildef12_file(intildef12)
    implicit none
    integer             :: i,j,ii,fun,number,m
    real(prec)          :: mm,lolzix
    real(prec)          :: intildef12(0:,0:)
    m=newhn
    intildef12 = 0._prec
    open(11,file='bas/file_t1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    fun = int((mm**2)/2)
    do ii = 1, fun
       read(11,rec=ii) number,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       intildef12(i,j) = lolzix
       ! print*, i,j,lolzix,preintildef12(i,j,2.5_prec)
    end do
    close(11)
  end subroutine read_intildef12_file

  subroutine create_interm1_file(intildef12)
    implicit none
    integer    :: counter,nm,m,s,t,u,other
    real(prec) :: temp
    real(prec) :: intildef12(0:,0:)
    m=newhn/2
    open(10,file='bas/file_t1i1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    counter = 0
    do other = 0, (m*2)-1
       do t = 0, nm
          do s = 0, t
             temp = 0._prec
             do u = 0, s
                if ((modulo(other+s+t-2*u,2).EQ.0)) then
                   ! write(10,rec=counter) i+ishft(j,8),preintildef12(i,j,2.5_prec)
                   temp = temp + coeff(u,s,t)*intildef12(other,s+t-2*u)/norm(s+t-2*u)
                end if
             end do
             counter = counter + 1
             write(10,rec=counter) other+ishft(s,8)+ishft(t,16),temp
          end do
       end do
    end do
    close(10)
  end subroutine create_interm1_file

  subroutine read_interm1_file(intinterm1)
    implicit none
    integer             :: other,s,t,ii,fun,number,m
    real(prec)          :: mm,lolzix
    real(prec)          :: intinterm1(0:,0:,0:)
    m=newhn/2
    intinterm1 = 0._prec
    open(11,file='bas/file_t1i1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    fun = int((mm**2+mm**3))
    do ii = 1, fun
       read(11,rec=ii) number,lolzix
       other=iand(number,255)
       s=iand(ishft(number,-8),255)
       t=iand(ishft(number,-16),255)
       intinterm1(other,s,t) = lolzix
    end do
    close(11)
  end subroutine read_interm1_file

  subroutine create_interm12_file(intinterm1)
    implicit none
    integer    :: counter,nm,m,w,i,j,s,t,summing
    real(prec) :: temp
    real(prec) :: intinterm1(0:,0:,0:)
    m=newhn/2
    open(10,file='bas/file_t1i2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    counter = 0
    do t = 0, nm
       do s = 0, t
          do j = 0, nm
             do i = 0, nm
                temp = 0._prec
                if (i.LT.j) then
                   summing = i
                else
                   summing = j
                endif
                do w = 0, summing
                   if (modulo(i+j-2*w+s+t,2).EQ.0) then
                      temp = temp + coeff(w,i,j)*intinterm1(i+j-2*w,s,t)/norm(i+j-2*w)
                   end if
                end do
                counter = counter + 1
                write(10,rec=counter) i+ishft(j,8)+ishft(s,16)+ishft(t,24),temp
             end do
          end do
       end do
    end do
    close(10)
  end subroutine create_interm12_file

  subroutine read_interm12_file(intinterm12)
    implicit none
    integer             :: i,j,s,t,ii,fun,number,m
    real(prec)          :: mm,lolzix
    real(prec)          :: intinterm12(0:,0:,0:,0:)
    m=newhn/2
    intinterm12 = 0._prec
    open(11,file='bas/file_t1i2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    fun = int((mm**4+mm**3)/2._prec)
    do ii = 1, fun
       read(11,rec=ii) number,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       s=iand(ishft(number,-16),255)
       t=ishft(number,-24)
       intinterm12(i,j,s,t) = lolzix
    end do
    close(11)
  end subroutine read_interm12_file

  real(prec) function inttildef12(mm,nn,oo,pp,ss,tt,intinterm12)
    implicit none
    integer    :: m,n,o,p,s,t,x,y
    integer    :: mm,nn,oo,pp,ss,tt
    integer    :: v1(4),v2(4)
    real(prec) :: tempf12,ttempf12
    real(prec) :: intinterm12(0:,0:,0:,0:)
    inttildef12 = 0._prec
    if (modulo(mm+nn+oo+pp+ss+tt,2).EQ.0) then
       if (ss.lt.tt) then
          s=ss
          t=tt
       else
          s=tt
          t=ss
       end if
       v1(1) = mm
       v1(2) = nn
       v1(3) = oo
       v1(4) = pp
       v2 = sortt(v1)
       m = v2(1)
       n = v2(2)
       o = v2(3)
       p = v2(4)
       tempf12 = 0._prec
       do y = 0, o
          ttempf12 = 0._prec
          do x = 0, m
             ttempf12 = ttempf12 + coeff(x,m,n)*intinterm12(m+n-2*x,o+p-2*y,s,t)
          end do
          tempf12 = tempf12 + coeff(y,o,p)*ttempf12
       end do
       inttildef12 = norm(mm)*norm(nn)*norm(oo)*norm(pp)*norm(ss)*norm(tt)*tempf12
    end if
  end function inttildef12

  real(prec) function imunutildef122norm(mu,nu,beta)
    implicit none
    integer    :: mu,nu,b
    real(prec) :: beta
    imunutildef122norm = 0._prec
    if (mod(mu+nu,2).NE.1) then
       b = bet(beta)
       imunutildef122norm = norm(nu)*norm(mu)/norm(mu+nu)*&
            &(sqrtpi*2._prec**((real(nu,prec)-1._prec)/2._prec)*3._prec**(-real(mu+nu,prec)/2._prec+1._prec)*&
            &(-1._prec)**(real(nu,prec)))*&
            &(2._prec*sqrt(real(mu+nu+2,prec)*real(mu+nu+1,prec))*&
            &0.25_prec*dzejmu(b,mu+nu+2)+(real(mu+nu,prec)+0.5_prec)*dzejmu(b,mu+nu)+&
            &0.5_prec*sqrt(real(mu+nu,prec)*(real(mu+nu,prec)-1._prec))*dzejmu(b,mu+nu-2))
    end if
  end function imunutildef122norm

  real(prec) function preintildef122(mu,nu,beta)
    implicit none
    integer    :: mu,nu,i
    real(prec) :: beta
    preintildef122 = 0._prec
    do i = 0, mu/2
       preintildef122 = preintildef122 + coo(i,mu)*imunutildef122norm(mu-2*i,nu,beta)/norm(mu-2*i)
    end do
    preintildef122 = preintildef122*norm(mu)
  end function preintildef122

  subroutine create_intildef122_file()
    implicit none
    integer    :: counter,i,j,nm,m
    m=newhn
    open(10,file='bas/file_t2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    counter = 0
    do j = 0, nm
       do i = 0, nm
          if ((modulo(i+j,2).EQ.0)) then
             counter = counter + 1
             write(10,rec=counter) i+ishft(j,8),preintildef122(i,j,4._prec)
          end if
       end do
    end do
    close(10)
  end subroutine create_intildef122_file

  subroutine read_intildef122_file(intildef122)
    implicit none
    integer             :: i,j,ii,fun,number,m
    real(prec)          :: mm,lolzix
    real(prec)          :: intildef122(0:,0:)
    m=newhn
    intildef122 = 0._prec
    open(11,file='bas/file_t2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    fun = int((mm**2)/2)
    do ii = 1, fun
       read(11,rec=ii) number,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       intildef122(i,j) = lolzix
       ! print*, i,j,lolzix,preintildef12(i,j,2.5_prec)
    end do
    close(11)
  end subroutine read_intildef122_file

  subroutine create_interm2_file(intildef122)
    implicit none
    integer    :: counter,nm,m,s,t,u,other
    real(prec) :: temp
    real(prec) :: intildef122(0:,0:)
    m=newhn/2
    open(10,file='bas/file_t2i1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    counter = 0
    do other = 0, (m*2)-1
       do t = 0, nm
          do s = 0, t
             temp = 0._prec
             do u = 0, s
                if ((modulo(other+s+t-2*u,2).EQ.0)) then
                   ! write(10,rec=counter) i+ishft(j,8),preintildef12(i,j,2.5_prec)
                   temp = temp + coeff(u,s,t)*intildef122(other,s+t-2*u)/norm(s+t-2*u)
                end if
             end do
             counter = counter + 1
             write(10,rec=counter) other+ishft(s,8)+ishft(t,16),temp
          end do
       end do
    end do
    close(10)
  end subroutine create_interm2_file

  subroutine read_interm2_file(intinterm2)
    implicit none
    integer             :: other,s,t,ii,fun,number,m
    real(prec)          :: mm,lolzix
    real(prec)          :: intinterm2(0:,0:,0:)
    m=newhn/2
    intinterm2 = 0._prec
    open(11,file='bas/file_t2i1.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    fun = int((mm**2+mm**3))
    do ii = 1, fun
       read(11,rec=ii) number,lolzix
       other=iand(number,255)
       s=iand(ishft(number,-8),255)
       t=iand(ishft(number,-16),255)
       intinterm2(other,s,t) = lolzix
    end do
    close(11)
  end subroutine read_interm2_file

  subroutine create_interm22_file(intinterm2)
    implicit none
    integer    :: counter,nm,m,w,i,j,s,t,summing
    real(prec) :: temp
    real(prec) :: intinterm2(0:,0:,0:)
    m=newhn/2
    open(10,file='bas/file_t2i2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    nm = m-1
    counter = 0
    do t = 0, nm
       do s = 0, t
          do j = 0, nm
             do i = 0, nm
                temp = 0._prec
                if (i.LT.j) then
                   summing = i
                else
                   summing = j
                endif
                do w = 0, summing
                   if ((modulo(i+j-2*w+s+t,2).EQ.0)) then
                      temp = temp + coeff(w,i,j)*intinterm2(i+j-2*w,s,t)/norm(i+j-2*w)
                   end if
                end do
                counter = counter + 1
                write(10,rec=counter) i+ishft(j,8)+ishft(s,16)+ishft(t,24),temp
             end do
          end do
       end do
    end do
    close(10)
  end subroutine create_interm22_file

  subroutine read_interm22_file(intinterm22)
    implicit none
    integer             :: i,j,s,t,ii,fun,number,m
    real(prec)          :: mm,lolzix
    real(prec)          :: intinterm22(0:,0:,0:,0:)
    m=newhn/2
    intinterm22 = 0._prec
    open(10,file='bas/file_t2i2.F',status='unknown',form='unformatted',access='direct',RECL=RECS)
    mm=real(m,prec)
    fun = int((mm**4+mm**3)/2._prec)
    do ii = 1, fun
       read(10,rec=ii) number,lolzix
       i=iand(number,255)
       j=iand(ishft(number,-8),255)
       s=iand(ishft(number,-16),255)
       t=ishft(number,-24)
       intinterm22(i,j,s,t) = lolzix
    end do
    close(10)
  end subroutine read_interm22_file

  real(prec) function inttildef122(mm,nn,oo,pp,ss,tt,intinterm22)
    implicit none
    integer    :: m,n,o,p,s,t,x,y
    integer    :: mm,nn,oo,pp,ss,tt
    integer    :: v1(4),v2(4)
    real(prec) :: tempf12,ttempf12
    real(prec) :: intinterm22(0:,0:,0:,0:)
    inttildef122 = 0._prec
    if (modulo(mm+nn+oo+pp+ss+tt,2).EQ.0) then
       if (ss.lt.tt) then
          s=ss
          t=tt
       else
          s=tt
          t=ss
       end if
       v1(1) = mm
       v1(2) = nn
       v1(3) = oo
       v1(4) = pp
       v2 = sortt(v1)
       m = v2(1)
       n = v2(2)
       o = v2(3)
       p = v2(4)
       tempf12 = 0._prec
       do y = 0, o
          ttempf12 = 0._prec
          do x = 0, m
             ttempf12 = ttempf12 + coeff(x,m,n)*intinterm22(m+n-2*x,o+p-2*y,s,t)
          end do
          tempf12 = tempf12 + coeff(y,o,p)*ttempf12
       end do
       inttildef122 = norm(mm)*norm(nn)*norm(oo)*norm(pp)*norm(ss)*norm(tt)*tempf12
    end if
  end function inttildef122

  ! subroutine make_dzejmu(beta)
  !   implicit none
  !   real(prec) :: beta
  !   integer    :: i,mu,b
  !   b = bet(beta)
  !   dzejmu(b,:)=0._prec
  !   do mu = 0, hn
  !      do i = 1, ggn
  !         dzejmu(b,mu) = dzejmu(b,mu) + hermiteh_gg(b,mu,i) * ggweights(i)
  !      end do
  !      dzejmu(b,mu) = 1._prec/sqrt(beta) * dzejmu(b,mu) * norm(mu)
  !   end do
  ! end subroutine make_dzejmu

  subroutine import_dzejmu()
    implicit none
    integer    :: mu,b
    open(10,file='bas/dz_2.5.dat')
    b = bet(2.5_prec)
    dzejmu(b,:)=0._prec
    do mu = 0, tho
       read(10,*) dzejmu(b,mu)
    end do
    close(10)
    open(10,file='bas/dz_3.dat')
    b = bet(3._prec)
    dzejmu(b,:)=0._prec
    do mu = 0, tho
       read(10,*) dzejmu(b,mu)
    end do
    close(10)
    open(10,file='bas/dz_4.dat')
    b = bet(4._prec)
    dzejmu(b,:)=0._prec
    do mu = 0, tho
       read(10,*) dzejmu(b,mu)
    end do
    close(10)
    open(10,file='bas/dz_5.dat')
    b = bet(5._prec)
    dzejmu(b,:)=0._prec
    do mu = 0, tho
       read(10,*) dzejmu(b,mu)
    end do
    close(10)
    ! print*, abs(dzejmu(1,1)-0.177042).GT.1.e-5
    ! print*, abs(dzejmu(2,1)-0.106225).GT.1.e-5
    ! stop
  end subroutine import_dzejmu

  integer function bet(beta)
    real(prec) :: beta
    if (abs(beta-3._prec).LT.1e-8) then
       bet = 1
    else if (abs(beta-5._prec).LT.1e-8) then
       bet = 2
    else if (abs(beta-2.5_prec).LT.1e-8) then
       bet = 3
    else if (abs(beta-4._prec).LT.1e-8) then
       bet = 4
    else
       print*, "Unsupported beta. Bombing out."
       call kurwout
    end if
  end function bet

  subroutine check(what,intildef12,intildef122,intinterm12,intinterm22,intf12)
    implicit none
    real (prec) :: intildef12(0:newhn,0:newhn),intildef122(0:newhn,0:newhn),intf12(0:,0:,0:,0:)
    real (prec) :: intinterm12(0:newhn,0:newhn,0:newhn,0:newhn),intinterm22(0:newhn,0:newhn,0:newhn,0:newhn)
    character(len=*), intent(in) :: what
    real(prec), parameter :: thr = 1.e-5
    select case (what)
    case ('dzejmu')
       if (abs(dzejmu(1,1)-0.333333).LT.thr) then
          print*, "dzejmu is probably not normalized, issues??"
       else if (abs(dzejmu(1,1)-0.177042).GT.thr.or.abs(dzejmu(2,1)-0.106225).GT.thr) then
          print*, abs(dzejmu(1,1)-0.177042),abs(dzejmu(2,1)-0.106225)
          print*, "normalized dzejmu is wrong"
          call kurwout()
       else
          write(LOUT,'(a)') "    Dzejmu is correct"
       end if
    case ('intf12')
       if (abs(intf12(3,3,3,3)-0.126686006).GT.thr.or.abs(intf12(1,2,3,4)-0.0280728333).GT.thr) then
          print*, "intf12 doesn't work",intf12(3,3,3,3),intf12(1,3,2,4)
          call kurwout()
       else
          write(LOUT,'(a)') "    Intf12 is correct"
       end if
    case ('intf122')
       if (abs(intf122(3,3,3,3,5._prec)-0.03901133673).GT.thr.or.abs(intf122(1,2,3,4,5._prec)-0.011738373).GT.thr) then
          print*, "intf122 doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    Intf122 is correct"
       end if
    case ('imunudf122norm')
       if (abs(imunudf122norm(3,3,5._prec)-0.086242001).GT.thr) then
          print*, "imunudf122norm doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    Imunudf122norm is correct"
       end if
    case ('intdf122')
       if (abs(intdf122(3,3,3,3,5._prec)-0.1662632851030596).GT.thr.or.&
            &abs(intdf122(1,2,3,4,5._prec)-0.03302934860107639).GT.thr) then
          print*, "intf122 doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    Intdf122 is correct"
       end if
    case ('imunutildef122norm')
       if (abs(imunutildef122norm(2,2,4._prec)+0.007343637523333).GT.thr.or.&
            &abs(imunutildef122norm(1,1,4._prec)-0.013847295710).GT.thr) then
          print*, "imunutildef122norm doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    imunutildef122 is correct"
       end if
    case ('imunutildef12norm')
       if (abs(imunutildef12norm(2,2,2.5_prec)+0.020014809331908905).GT.thr.or.&
            &abs(imunutildef12norm(1,1,2.5_prec)-0.03265986323710979).GT.thr) then
          print*, "imunutildef12norm doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    imunutildef12 is correct"
       end if
    case ('intildef12')
       if (abs(intildef12(3,3)+0.020425369523330).GT.thr.or.&
            &abs(intildef12(0,0)-0.3464101615137748).GT.thr) then
          print*, "intildef12 doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    intildef12 is correct"
       end if
    case ('inttildef12')
       if (abs(inttildef12(2,2,2,2,2,2,intinterm12)-0.03541268534).GT.thr.or.&
            &abs(inttildef12(1,5,6,2,7,5,intinterm12)+0.00090355492).GT.thr) then
          print*, "2,2,2,2,2,2",inttildef12(2,2,2,2,2,2,intinterm12)
          print*, "1,5,6,2,7,5",inttildef12(1,5,6,2,7,5,intinterm12)
          print*, "inttildef12 doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    inttildef12 is correct"
       end if
    case ('intildef122')
       if (abs(intildef122(3,3)+0.00826159).GT.thr.or.&
            &abs(intildef122(0,0)-0.117498).GT.thr) then
          print*, "intildef122 doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    intildef122 is correct"
       end if
    case ('inttildef122')
       if (abs(inttildef122(2,2,2,2,2,2,intinterm22)-0.010651100993).GT.thr.or.&
            &abs(inttildef122(1,5,6,2,7,5,intinterm22)+0.00039196641).GT.thr) then
          print*, "inttildef122 doesn't work"
          call kurwout()
       else
          write(LOUT,'(a)') "    inttildef122 is correct"
       end if
    end select
  end subroutine check

  ! function dzejmu(mu,beta)
  !   implicit none
  !   integer    :: i,mu
  !   real(prec) :: dzejmu,beta
  !   dzejmu = 0._prec
  !   if (mu.lt.0) then
  !      dzejmu = 0._prec
  !   else
  !      do i = 1, ggn
  !         dzejmu = dzejmu + hermiteh_gg(mu,i) * ggweights(i)
  !      end do
  !      dzejmu = 1._prec/sqrt(beta) * dzejmu * norm(mu)
  !   end if
  ! end function dzejmu

  subroutine kurwout()
    print*, "__/\\\________/\\\__/\\\________/\\\____/\\\\\\\\\______/\\\______________/\\\_____/\\\\\\\\\____       "
    print*, " _\/\\\_____/\\\//__\/\\\_______\/\\\__/\\\///////\\\___\/\\\_____________\/\\\___/\\\\\\\\\\\\\__      "
    print*, "  _\/\\\__/\\\//_____\/\\\_______\/\\\_\/\\\_____\/\\\___\/\\\_____________\/\\\__/\\\/////////\\\_     "
    print*, "   _\/\\\\\\//\\\_____\/\\\_______\/\\\_\/\\\\\\\\\\\/____\//\\\____/\\\____/\\\__\/\\\_______\/\\\_    "
    print*, "    _\/\\\//_\//\\\____\/\\\_______\/\\\_\/\\\//////\\\_____\//\\\__/\\\\\__/\\\___\/\\\\\\\\\\\\\\\_   "
    print*, "     _\/\\\____\//\\\___\/\\\_______\/\\\_\/\\\____\//\\\_____\//\\\/\\\/\\\/\\\____\/\\\/////////\\\_  "
    print*, "      _\/\\\_____\//\\\__\//\\\______/\\\__\/\\\_____\//\\\_____\//\\\\\\//\\\\\_____\/\\\_______\/\\\_ "
    print*, "       _\/\\\______\//\\\__\///\\\\\\\\\/___\/\\\______\//\\\_____\//\\\__\//\\\______\/\\\_______\/\\\_"
    print*, "        _\///________\///_____\/////////_____\///________\///_______\///____\///_______\///________\///_"
    stop
  end subroutine kurwout

end module ggint
