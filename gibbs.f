       subroutine gibbs( iseed )          
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts,tseq(5,n)
       real*8 d(n),u(n),v(n),e1star(n),e2star(n)
       integer ncase(6),idc(n),ng(5)
       real*8 case(6,n),x(np,n)
       real*8 s(5,ngmax),p(ngmax)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n)
       real*8 alam(5,ngmax),tau,taus
       real*8 bo(5,2)
       integer iseed
       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
c      common variables
       common /vid/idc
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vtseq/tseq
       common /vecd/d
       common /vecu/u
       common /vecv/v
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vtaus/taus
       common /vecbo/bo
       common /vecng/ng       

c       write(*,*) 'before Gtaus'
       call Gtaun( iseed )
       write(*,*) 'after Gtaus'

c       write(*,*) 'before Gomega'       
       call Gomega( iseed )
       write(*,*) 'after Gomega'

c       write(*,*) 'before Galam'
       call Galam( iseed )
       write(*,*) 'after Galam'

c       write(*,*) 'before Ge1star'
       call Ge1starnew( iseed )
       write(*,*) 'after Ge1star'

c       write(*,*) 'before Ge2star'
       call Ge2starnew( iseed )
       write(*,*) 'after Ge2star'

c       write(*,*) 'before Ggam'    
       call Ggam( iseed )
       write(*,*) 'after Ggam' 
       
c       write(*,*) 'before Galpha'         
       call Galpha( iseed )
       write(*,*) 'after Galpha'        

c       write(*,*) 'before Gbeta'
       call Gbeta( iseed )
       write(*,*) 'after Gbeta'

       return
       end
   

       subroutine Galam( iseed )
c      This subroutine samples \lambda_{kj}               
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n)
       real*8 ts(n),tseq(5,n)
       real*8 e1star(n),e2star(n)
       real*8 case(6,n),x(np,n)
       real*8 s(5,ngmax)
       integer ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2)
       real*8 shape(5,ngmax),scale(5,ngmax),rg(1)       
       integer iseed,iaccept(n)
       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vtseq/tseq
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vtaus/taus
       common /vecbo/bo
       common /vecng/ng
c
       do k=1,5
        do j=1,ng(k)
         shape(k,j)=bo(k,1)
         scale(k,j)=bo(k,2)
        enddo
c       
        do 100 i=1,n
         ts(i)=tseq(k,i)
c               
         if (ts(i) .eq. 9999) goto 100
         ein=dotproduct(x(:,i),gam(k,:),np)
         if(k .eq. 1) then
          ein2=case(1,i)+(1.0d0-e1star(i))*case(6,i)
          ein3=case(1,i)
         else if(k .eq. 2) then
          ein2=omega(i)*(case(2,i)+case(3,i)+case(4,i)
     1        +case(5,i)+e1star(i)*case(6,i))
          ein3=case(2,i)+case(3,i)+case(4,i)+case(5,i)
         else if(k .eq. 3) then
          ein2=omega(i)*(case(2,i)+(1.0d0-e2star(i))*case(5,i))
          ein3=case(2,i)
         else if(k .eq. 4) then
          ein2=omega(i)*(case(3,i)+case(4,i)
     1       +e2star(i)*case(5,i))
          ein3=case(3,i)+case(4,i)
         else if(k .eq. 5) then
          ein2=omega(i)*(case(3,i)+case(4,i))
          ein3=case(3,i)
         endif
c
         if (ts(i) .le. s(k,1)) jg=1
         do j1=1,ng(k)-1
          if ( (ts(i) .gt. s(k,j1)) .and.
     1       (ts(i) .le. s(k,j1+1)) ) then
             jg=j1+1
          endif
         enddo
         if(jg .eq. 1) then
          scale(k,1)=scale(k,1)+ts(i)*dexp(ein)*ein2
          shape(k,1)=shape(k,1)+ein3
         else
          scale(k,1)=scale(k,1)+s(k,1)*dexp(ein)*ein2
          shape(k,jg)=shape(k,jg)+ein3
          do j=2,jg-1
           scale(k,j)=scale(k,j)+(s(k,j)-s(k,j-1))*dexp(ein)*ein2
          enddo
          scale(k,jg)=scale(k,jg)+(ts(i)-s(k,jg-1))*dexp(ein)*ein2
         endif
100     continue
c
        do j=1,ng(k)
         shape1=shape(k,j)
         call rnset( iseed )
         call DRNGAM(1,shape1,rg)
         call rnget( iseed )
         alam(k,j)=rg(1)/scale(k,j)
        enddo
       enddo
       return
       end              

    
       subroutine Gomega( iseed )
c      This subroutine sample frailty term \omega_i
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n)
       real*8 ts(n),tseq(5,n)
       real*8 e1star(n),e2star(n)
       real*8 case(6,n),x(np,n)
       integer ng(5)
       real*8 s(5,ngmax)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2),ss(5),xgam,alogss
       real*8 shape(5,ngmax),scale(5,ngmax),rg(1)       
       integer iseed,iaccept(n)
       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vtseq/tseq
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vtaus/taus
       common /vecbo/bo
       common /vecng/ng
       taus=1.0d0/tau
       do 100 i=1,n
        if(case(1,i) .eq. 1.0d0) goto 100   !04-25-2023
        do k=1,5
         ss(k)=0.0d0
        enddo
        a1=(2.0d0+taus)*(case(2,i)+case(4,i))
        a2=(3.0d0+taus)*case(3,i)
        a3=(1.0d0+taus)*e2star(i)*case(5,i)
        a4=(taus)*e1star(i)*case(6,i)       
        shape1=a1+a2+a3+a4
        if(shape1 .eq. 0.0d0) goto 100
c       as not all subjects will have frailty term        
c
        do 101 k=2,5
         ts(i)=tseq(k,i)     
c
         if(ts(i) .ne. 9999) then    
          call alogsurv(k,i,ts(i),gam,alam,xgam,alogss,ahaz)
          ss(k)=-alogss
         endif
101     continue
        b1=(1.0d0+ss(2)+ss(3))*(case(2,i)+(1-e2star(i))*case(5,i))
        b2=(1.0d0+ss(2)+ss(4)+ss(5))*(case(3,i)+case(4,i))
        b3=(1.0d0+ss(2)+ss(4))*(e2star(i)*case(5,i))
        b4=(1.0d0+ss(2))*(e1star(i)*case(6,i))
        scale1=b1+b2+b3+b4
        if(scale1 .eq. 0.0d0) goto 100        
c       not all subjects have frailty term
        call rnset( iseed )
        call DRNGAM(1,shape1,rg)
        call rnget( iseed )
        omega(i)=rg(1)/scale1
100    continue
       write(*,*) 'min omega',minval(omega),'max omega',maxval(omega)
       return
       end       


       subroutine Ggam( iseed )
c      This subroutine sample \gamma parameters              
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)       
       PARAMETER(maxtries=25,im=2,NS=10,in=1000)       
       real*8 emax,step
C ------------------------ Working variables ------------------------------
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1         xlb,xub,rwv(6*NS+15+in)
       real*8 start(1),xmin(1),step1(1)
       integer ifault,iwv(NS+np+in),ifcount
       integer iseed,idum2       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n)
       real*8 e1star(n),e2star(n)
       real*8 case(6,n),x(np,n)
       real*8 s(5,ngmax)
       integer ncase(6),idc(n),ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2)
       integer iaccept(n)
       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
c      common variables
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vecbo/bo
       common /vecng/ng
       common /dummy2/idum2
       common /dummy3/idum3
       common /dummy1C/idumC1
       common /dummy2C/idumC2
c      external  DRNUNF
        external evalg,fgam
       emax=60.0d0
       do 79 k=1,5
        idumC1=k
        DO 80 j=1,np
         idumC2=j
         step=0.5d0
c         nopt = 1
c         reqmin = 1.0d-10
c         konvge = 5
c         kcount = 1000
c         step1(1)=0.2d0
c         start(1)=gam(k,j)
c         call nelmin(fgam, nopt, start, xmin, ynewlo, reqmin, step1,
c     1         konvge, kcount, icount, numres, ifault1) 
        a(1)=gam(k,j)-step
        a(2)=gam(k,j)+step
c         a(1)=xmin(1)-step
c         a(2)=xmin(1)+step
c         CALL evalg(xmin(1), af, adf)
c         write(*,201) 'j',j,'xmin(1)',xmin(1),'af',af,'adf',adf
c201      format(A,1x,I4,1x,A,1x,f16.3,1x,A,1x,f16.3,1x,A,1x,f16.3)         
20       CALL evalg(a(1), ha(1), hpa(1))
         CALL evalg(a(2), ha(2), hpa(2))       
c ------------------ Get starting points for Gilks -----------------------
         IF (hpa(2) .gt. 0.0) THEN
c        Both points to the left of the mode; push a(2) to the right
30        a(2)=a(2)+step
          CALL EVALg(a(2), ha(2), hpa(2))
c         If a(2) still to the left of the mode repeat, else continue
          IF (hpa(2).gt.0.0) THEN
           GOTO 30
          ELSE
           a(1)=a(2)-step
          ENDIF
         ELSE
          IF (hpa(1) .le. 0.0) THEN
c        Both points to the right of the mode
c        Insert new point to the left of a(1)
           a(2)=a(1)
           ha(2)=ha(1)
           hpa(2)=hpa(1)
40         a(1)=a(1)-step
           CALL EVALg(a(1), ha(1), hpa(1))
           IF (hpa(1) .le. 0.0) THEN
            GOTO 40
c          ELSE
c           a(2)=a(1)+step
           ENDIF
          ENDIF
         ENDIF
c         write(*,*) 'k',k,'j',j,'passed1'              
c -----------------------------------------------------------------------
              ifcount=0
50       CALL INITIAL(ns, im, emax, a, ha, hpa, lb, xlb, ub, xub,
     +                ifault, iwv, rwv, evalg)
         IF((ifault.eq.3) .or. (ifault.eq.4)) GOTO 20
C
         CALL SAMPLE(iwv,rwv,evalg,val,ifault,iseed)
         IF(ifault.eq.5) THEN
c          PRINT*, 'ifault =', ifault, ' trying again'
          a(1)=a(1)-step
          a(2)=a(2)+step
          ifcount=ifcount+1
          IF(ifcount.le.maxtries) THEN
           GOTO 50
          ELSE
c             PRINT*, 'Too many tries resetting starting values'
           a(1)=0.0
           a(2)=0.0
           step=1.2d0*step
           GOTO 20
          ENDIF
         ENDIF
         IF(ifault.eq.7) THEN
          a(1)=a(1)-step
          a(2)=a(2)+step
          ifcount=ifcount+1
          IF(ifcount.le.maxtries) THEN
           GOTO 50
          ELSE
c            PRINT*, 'Too many tries resetting starting values'
           a(1)=0.0
           a(2)=0.0
           step=1.2d0*step
           GOTO 20
          ENDIF
         ENDIF
         IF (IFAULT .NE. 0) THEN
c         PRINT*, ' The sampling status is ', IFAULT
         ENDIF
C
         gam(k,j)=val
C
80      continue
79     continue
                 
       return
       end 
                            
       real*8 function fgam(z1)
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)
       PARAMETER(maxtries=25,im=2,NS=10,in=1000)       
       real*8 emax,step
C ------------------------ Working variables ------------------------------
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1         xlb,xub,rwv(6*NS+15+in)
       real*8 start(1),xmin(1),step1(1)
       integer ifault,iwv(NS+np+in),ifcount
       integer iseed,idum2       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n),tseq(5,n)
       real*8 e1star(n),e2star(n)
       real*8 case(6,n),x(np,n)
       real*8 s(5,ngmax),xgam
       integer ncase(6),idc(n),ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax)
       real*8 bo(5,2)
       integer iaccept(n)
       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vtseq/tseq
       common /vcase/case
       common /vxold/xold
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vecbo/bo
       common /vecng/ng
       common /dummy2/idum2
       common /dummy3/idum3
       common /dummy1C/idumC1
       common /dummy2C/idumC2
c      external  DRNUNF
       k=idumC1
       gam(k,idumC2)=z1
c......calculate log likelihood and its derivative
       sum1=0.0d0
       sum2=0.0d0
       do 300 i=1,n
         ts(i)=tseq(k,i)     
c           
        if(ts(i) .ne. 9999) then    
         call alogsurv(k,i,ts(i),gam,alam,xgam,ss,ahaz)
         if(k .eq. 1) then
          a1=case(1,i)
          a2=a1+(1.0d0-e1star(i))*case(6,i)              
          ss1=ss
         else if(k .eq. 2) then
          a1=case(2,i)+case(3,i)+case(4,i)+case(5,i)
          a2=a1+e1star(i)*case(6,i)              
          ss1=ss*omega(i)
         else if(k .eq. 3) then
          a1=case(2,i)
          a2=a1+(1.0d0-e2star(i))*case(5,i)              
          ss1=ss*omega(i)
         else if(k .eq. 4) then
          a1=case(3,i)+case(4,i)
          a2=a1+e2star(i)*case(5,i)
          ss1=ss*omega(i)
         else if(k .eq. 5) then
          a1=case(3,i)
          a2=a1+case(4,i)
          ss1=ss*omega(i)
         endif 
         sum1=sum1+a1*xgam+a2*ss1
        endif
 300   continue
       if((k .eq. 1)) then        
        sum1=sum1-z1**2/(2.0d0*100.0d0)
       else
        sum1=sum1-z1**2/(2.0d0*100.0d0)
       endif              
       fgam=-sum1
       return
       end

              
       subroutine evalg(z1,ahz,ahpz)
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)       
       PARAMETER(maxtries=25,im=2,NS=10,in=1000)       
       real*8 emax,step
C ------------------------ Working variables ------------------------------
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1         xlb,xub,rwv(6*NS+15+in)
       real*8 start(1),xmin(1),step1(1)
       integer ifault,iwv(NS+np+in),ifcount
       integer iseed,idum2       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n),tseq(5,n)
       real*8 e1star(n),e2star(n)
       real*8 case(6,n),x(np,n)
       real*8 s(5,ngmax),xgam
       integer ncase(6),idc(n),ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2)
       integer iaccept(n)
       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vtseq/tseq
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vtaus/taus
       common /vecbo/bo
       common /vecng/ng
       common /dummy2/idum2
       common /dummy3/idum3
       common /dummy1C/idumC1
       common /dummy2C/idumC2
c       external  DRNUNF
        k=idumC1
        gam(k,idumC2)=z1
c......calculate log likelihood and its derivative
       sum1=0.0d0
       sum2=0.0d0
       do 300 i=1,n
        ts(i)=tseq(k,i)
c            
        if(ts(i) .ne. 9999) then    
         call alogsurv(k,i,ts(i),gam,alam,xgam,ss,ahaz)
         if(k .eq. 1) then
          a1=case(1,i)
          a2=a1+(1.0d0-e1star(i))*case(6,i)              
          ss1=ss
         else if(k .eq. 2) then
          a1=case(2,i)+case(3,i)+case(4,i)+case(5,i)
          a2=a1+e1star(i)*case(6,i)              
          ss1=ss*omega(i)
         else if(k .eq. 3) then
          a1=case(2,i)
          a2=a1+(1.0d0-e2star(i))*case(5,i)              
          ss1=ss*omega(i)
         else if(k .eq. 4) then
          a1=case(3,i)+case(4,i)
          a2=a1+e2star(i)*case(5,i)
          ss1=ss*omega(i)
         else if(k .eq. 5) then
          a1=case(3,i)
          a2=a1+case(4,i)
          ss1=ss*omega(i)
         endif 
         sum1=sum1+a1*xgam+a2*ss1
         sum2=sum2+(a1+a2*ss1)*x(idumC2,i)
        endif
300    continue
       if(k .eq. 1) then
        sum1=sum1-z1**2/(2.0d0*100.0d0)
        sum2=sum2-z1/3.0d0
       else         
        sum1=sum1-z1**2/(2.0d0*100.0d0)
        sum2=sum2-z1/100.0d0
       endif
       ahz=sum1
       ahpz=sum2  
       return
       end
  
 
       subroutine Galpha( iseed )
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)       
       PARAMETER(maxtries=25,im=2,NS=10,in=1000)       
       real*8 emax,step
C ------------------------ Working variables ------------------------------
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1         xlb,xub,rwv(6*NS+15+in)
       real*8 start(1),xmin(1),step1(1)
       integer ifault,iwv(NS+np+in),ifcount
       integer iseed,idum2       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n)
       real*8 e1star(n),e2star(n)
       real*8 case(6,n),x(np,n)
       integer ncase(6),idc(n),ng(5)
       real*8 s(5,ngmax),p(ngmax)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2)
       integer iaccept(n)
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vtaus/taus
       common /vecbo/bo
       common /vecng/ng
       common /dummy2C/idumC2
c      external  DRNUNF
       external evala,falpha
       emax=60.0d0
        DO 80 j=1,np+1
         idumC2=j
         step=0.5d0
c         nopt = 1
c         reqmin = 1.0d-10
c         konvge = 5
c         kcount = 1000
c         step1(1)=0.2d0
c         start(1)=alpha(j)
c         call nelmin(falpha, nopt, start, xmin, ynewlo, reqmin, step1,
c     1         konvge, kcount, icount, numres, ifault1) 
        a(1)=alpha(j)-step
        a(2)=alpha(j)+step
c         a(1)=xmin(1)-step
c         a(2)=xmin(1)+step
c         CALL evala(xmin(1), af, adf)
c         write(*,201) 'j',j,'xmin(1)',xmin(1),'af',af,'adf',adf
201      format(A,1x,I4,1x,A,1x,f16.3,1x,A,1x,f16.3,1x,A,1x,f16.3)           
20       CALL evala(a(1), ha(1), hpa(1))
         CALL evala(a(2), ha(2), hpa(2))
c ------------------ Get starting points for Gilks -----------------------
         IF (hpa(2) .gt. 0.0) THEN
c        Both points to the left of the mode; push a(2) to the right
30        a(2)=a(2)+step
          CALL EVALa(a(2), ha(2), hpa(2))
c         If a(2) still to the left of the mode repeat, else continue
          IF (hpa(2) .gt.  0.0) THEN
           GOTO 30
          ELSE
           a(1)=a(2)-step
          ENDIF
         ELSE
          IF (hpa(1) .le. 0.0) THEN
c        Both points to the right of the mode
c        Insert new point to the left of a(1)
           a(2)=a(1)
           ha(2)=ha(1)
           hpa(2)=hpa(1)
40         a(1)=a(1)-step
           CALL EVALa(a(1), ha(1), hpa(1))
           IF (hpa(1) .le. 0.0) THEN
            GOTO 40
c          ELSE
c           a(2)=a(1)+step
           ENDIF
          ENDIF
         ENDIF
c         write(*,*) 'j',j,'passed1'              
c -----------------------------------------------------------------------
              ifcount=0
50       CALL INITIAL(ns, im, emax, a, ha, hpa, lb, xlb, ub, xub,
     +                ifault, iwv, rwv, evala)
         IF((ifault.eq.3) .or. (ifault.eq.4)) GOTO 20
C
         CALL SAMPLE(iwv,rwv,evala,val,ifault,iseed)
         IF(ifault.eq.5) THEN
c          PRINT*, 'ifault =', ifault, ' trying again'
          a(1)=a(1)-step
          a(2)=a(2)+step
          ifcount=ifcount+1
          IF(ifcount.le.maxtries) THEN
           GOTO 50
          ELSE
c             PRINT*, 'Too many tries resetting starting values'
           a(1)=0.0
           a(2)=0.0
           step=1.2d0*step
           GOTO 20
          ENDIF
         ENDIF
         IF(ifault.eq.7) THEN
          a(1)=a(1)-step
          a(2)=a(2)+step
          ifcount=ifcount+1
          IF(ifcount.le.maxtries) THEN
           GOTO 50
          ELSE
c            PRINT*, 'Too many tries resetting starting values'
           a(1)=0.0
           a(2)=0.0
           step=1.2d0*step
           GOTO 20
          ENDIF
         ENDIF
         IF (IFAULT .NE. 0) THEN
c         PRINT*, ' The sampling status is ', IFAULT
         ENDIF
C
         alpha(j)=val
C
80      continue
                 
       return
       end 
                            
       real*8 function falpha(z1)
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)
       PARAMETER(maxtries=25,im=2,NS=10,in=1000)       
       real*8 emax,step
C ------------------------ Working variables ------------------------------
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1         xlb,xub,rwv(6*NS+15+in)
       real*8 start(1),xmin(1),step1(1)
       integer ifault,iwv(NS+np+in),ifcount
       integer iseed,idum2       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n)
       real*8 e1star(n),e2star(n)
       real*8 case(6,n),xx(np+1),x(np,n)
       real*8 s(5,ngmax)
       integer ncase(6),idc(n),ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2)
c       real*8 shape(2,ngmax),scale(2,ngmax),rg(1)       
       integer iaccept(n)
       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vtaus/taus
       common /vecbo/bo
       common /vecng/ng
       common /dummy2/idum2
       common /dummy3/idum3
       common /dummy1C/idumC1
       common /dummy2C/idumC2
c       external  DRNUNF
       alpha(idumC2)=z1
c......calculate log likelihood and its derivative
       sum1=0.0d0
       sum2=0.0d0
       do 300 i=1,n
        xx(1)=1.0d0
        do j=2,np+1
         xx(j)=x(j-1,i)
        enddo
        ein=dotproduct(xx,alpha,np+1)
        a1=case(1,i)+(1.0d0-e1star(i))*case(6,i)
        a2=case(2,i)+case(3,i)+case(4,i)
     1   +case(5,i)+e1star(i)*case(6,i)
        sum1=sum1+a1*dlog(1.0d0/(1.0d0+dexp(ein)))
     1           +a2*dlog(1.0d0/(1.0d0+dexp(-ein)))
300    continue
       if (idumC2 .eq. 1) then
        sum1=sum1-z1**2/(2.0d0*100.0d0)
       else
        sum1=sum1-z1**2/(2.0d0*100.0d0)
       endif
       falpha=-sum1 
       return
       end

 
       subroutine evala(z1,ahz,ahpz)
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)
       PARAMETER(maxtries=25,im=2,NS=10,in=1000)       
       real*8 emax,step
C ------------------------ Working variables ------------------------------
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1         xlb,xub,rwv(6*NS+15+in)
       real*8 start(1),xmin(1),step1(1)
       integer ifault,iwv(NS+np+in),ifcount
       integer iseed,idum2       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n)
       real*8 e1star(n),e2star(n)
       real*8 case(6,n),x(np,n),xx(np+1)
       real*8 s(5,ngmax)
       integer ncase(6),idc(n),ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2)
       integer iaccept(n)
       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
c      common variables
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vcase/case
       common /vxold/xold
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vtaus/taus
       common /vecbo/bo
       common /vecng/ng
       common /dummy2/idum2
       common /dummy3/idum3
       common /dummy1C/idumC1
       common /dummy2C/idumC2
c       external  DRNUNF
       alpha(idumC2)=z1
c......calculate log likelihood and its derivative
       sum1=0.0d0
       sum2=0.0d0
       do 300 i=1,n
        xx(1)=1.0d0
        do j=2,np+1
         xx(j)=x(j-1,i)
        enddo
        ein=dotproduct(xx,alpha,np+1)
        a1=case(1,i)+(1-e1star(i))*case(6,i)
        a2=case(2,i)+case(3,i)+case(4,i)
     1   +case(5,i)+e1star(i)*case(6,i)
        sum1=sum1+a1*dlog(1.0d0/(1.0d0+dexp(ein)))
     1           +a2*dlog(1.0d0/(1.0d0+dexp(-ein)))
        sum2=sum2-a1*xx(idumC2)*dexp(ein)/(1.0d0+dexp(ein))
     1           +a2*xx(idumC2)*dexp(-ein)/(1.0d0+dexp(-ein))        
300    continue
       if (idumC2 .eq. 1) then
        sum1=sum1-z1**2/(2.0d0*100.0d0)
        sum2=sum2-z1/100.0d0
       else
        sum1=sum1-z1**2/(2.0d0*100.0d0)
        sum2=sum2-z1/100.0d0
       endif
       ahz=sum1
       ahpz=sum2  
       return
       end

       subroutine Gbeta( iseed )
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)
       PARAMETER(maxtries=25,im=2,NS=10,in=1000)       
       real*8 emax,step
C ------------------------ Working variables ------------------------------
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1         xlb,xub,rwv(6*NS+15+in)
       real*8 start(1),xmin(1),step1(1)
       integer ifault,iwv(NS+np+in),ifcount
       integer iseed,idum2       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n)
       real*8 e1star(n),e2star(n)
       real*8 case(6,n),x(np,n)
       real*8 s(5,ngmax)
       integer ncase(6),idc(n),ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2)
       integer iaccept(n)
c       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vecbo/bo
       common /vecng/ng
       common /dummy2/idum2
       common /dummy3/idum3
       common /dummy1C/idumC1
       common /dummy2C/idumC2
c      external  DRNUNF
       external evalb,fbeta
       emax=60.0d0
        DO 80 j=1,np+1
         idumC2=j
         step=0.30d0
c         nopt = 1
c         reqmin = 1.0d-10
c         konvge = 5
c         kcount = 1000
c         step1(1)=0.2d0
c         start(1)=beta(j)
c         call nelmin(fbeta, nopt, start, xmin, ynewlo, reqmin, step1,
c     1         konvge, kcount, icount, numres, ifault1) 
        a(1)=beta(j)-step
        a(2)=beta(j)+step
c         a(1)=xmin(1)-step
c         a(2)=xmin(1)+step
c         CALL evalb(xmin(1), af, adf)
c         write(*,201) 'j',j,'xmin(1)',xmin(1),'af',af,'adf',adf
c201      format(A,1x,I4,1x,A,1x,f16.3,1x,A,1x,f16.3,1x,A,1x,f16.3)           
20       CALL evalb(a(1), ha(1), hpa(1))
         CALL evalb(a(2), ha(2), hpa(2))
c        write(*,*) 'a(1)=',a(1),' ha(1)=',ha(1)
c        write(*,*) ' hpa(1)=',hpa(1)
c        write(*,*) 'a(2)=',a(2),' ha(2)=',ha(2)
c        write(*,*) ' hpa(2)=',hpa(2)        
c ------------------ Get starting points for Gilks -----------------------
         IF (hpa(2) .gt. 0.0) THEN
c        Both points to the left of the mode; push a(2) to the right
30        a(2)=a(2)+step
          CALL EVALb(a(2), ha(2), hpa(2))
c         If a(2) still to the left of the mode repeat, else continue
          IF (hpa(2).gt.0.0) THEN
           GOTO 30
          ELSE
           a(1)=a(2)-step
          ENDIF
         ELSE
          IF (hpa(1) .le. 0.0) THEN
c        Both points to the right of the mode
c        Insert new point to the left of a(1)
           a(2)=a(1)
           ha(2)=ha(1)
           hpa(2)=hpa(1)
40         a(1)=a(1)-step
           CALL EVALb(a(1), ha(1), hpa(1))
           IF (hpa(1) .le. 0.0) THEN
            GOTO 40
c          ELSE
c           a(2)=a(1)+step
           ENDIF
          ENDIF
         ENDIF
c         write(*,*) 'j',j,'passed1'              
c -----------------------------------------------------------------------
              ifcount=0
50       CALL INITIAL(ns, im, emax, a, ha, hpa, lb, xlb, ub, xub,
     +                ifault, iwv, rwv, evalb)
         IF((ifault.eq.3) .or. (ifault.eq.4)) GOTO 20
C
         CALL SAMPLE(iwv,rwv,evalb,val,ifault,iseed)
         IF(ifault.eq.5) THEN
c          PRINT*, 'ifault =', ifault, ' trying again'
          a(1)=a(1)-step
          a(2)=a(2)+step
          ifcount=ifcount+1
          IF(ifcount.le.maxtries) THEN
           GOTO 50
          ELSE
c             PRINT*, 'Too many tries resetting starting values'
           a(1)=0.0
           a(2)=0.0
           step=1.2d0*step
           GOTO 20
          ENDIF
         ENDIF
         IF(ifault.eq.7) THEN
          a(1)=a(1)-step
          a(2)=a(2)+step
          ifcount=ifcount+1
          IF(ifcount.le.maxtries) THEN
           GOTO 50
          ELSE
c            PRINT*, 'Too many tries resetting starting values'
           a(1)=0.0
           a(2)=0.0
           step=1.2d0*step
           GOTO 20
          ENDIF
         ENDIF
         IF (IFAULT .NE. 0) THEN
c         PRINT*, ' The sampling status is ', IFAULT
         ENDIF
C
         beta(j)=val
80      continue
       return
       end 
                            
       real*8 function fbeta(z1)
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)
       PARAMETER(maxtries=25,im=2,NS=10,in=1000)       
       real*8 emax,step
C ------------------------ Working variables ------------------------------
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1         xlb,xub,rwv(6*NS+15+in)
       real*8 start(1),xmin(1),step1(1)
       integer ifault,iwv(NS+np+in),ifcount
       integer iseed,idum2       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n)
       real*8 e1star(n),e2star(n)
       real*8 case(6,n),x(np,n),xx(np+1)
       real*8 s(5,ngmax)
       integer ncase(6),idc(n),ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2)
       integer iaccept(n)
       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vecbo/bo
       common /vecng/ng
       common /dummy2/idum2
       common /dummy3/idum3
       common /dummy1C/idumC1
       common /dummy2C/idumC2
c       external  DRNUNF
       beta(idumC2)=z1
c......calculate log likelihood and its derivative
       sum1=0.0d0
       do 300 i=1,n
        xx(1)=1.0d0
        do j=2,np+1
         xx(j)=x(j-1,i)
        enddo
        ein=dotproduct(xx,beta,np+1)
        a1=case(2,i)+(1.0d0-e2star(i))*case(5,i)
        a2=case(3,i)+case(4,i)+e2star(i)*case(5,i)
        sum1=sum1+a1*dlog(1.0d0/(1.0d0+dexp(ein)))
     1           +a2*dlog(1.0d0/(1.0d0+dexp(-ein)))
300    continue
       if (idumC2 .eq. 1) then
        sum1=sum1-z1**2/(2.0d0*10.0d0)
       else
        sum1=sum1-z1**2/(2.0d0*100.0d0)
       endif
       fbeta=-sum1 
       return
       end

 
       subroutine evalb(z1,ahz,ahpz)
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)
       PARAMETER(maxtries=25,im=2,NS=10,in=1000)       
       real*8 emax,step
C ------------------------ Working variables ------------------------------
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1         xlb,xub,rwv(6*NS+15+in)
       real*8 start(1),xmin(1),step1(1)
       integer ifault,iwv(NS+np+in),ifcount
       integer iseed,idum2       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n)
       real*8 e1star(n),e2star(n)
       real*8 case(6,n),x(np,n),xx(np+1)
       real*8 s(5,ngmax)
       integer ncase(6),idc(n),ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2)
       integer iaccept(n)
       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vtaus/taus
       common /vecbo/bo
       common /vecng/ng
       common /dummy2/idum2
       common /dummy3/idum3
       common /dummy1C/idumC1
       common /dummy2C/idumC2
c       external  DRNUNF
       beta(idumC2)=z1
c......calculate log likelihood and its derivative
       sum1=0.0d0
       sum2=0.0d0
       do 300 i=1,n
        xx(1)=1.0d0
        do j=2,np+1
         xx(j)=x(j-1,i)
        enddo
        ein=dotproduct(xx,beta,np+1)
        a1=case(2,i)+(1.0d0-e2star(i))*case(5,i)
        a2=case(3,i)+case(4,i)+e2star(i)*case(5,i)
        sum1=sum1+a1*dlog(1.0d0/(1.0d0+dexp(ein)))
     1           +a2*dlog(1.0d0/(1.0d0+dexp(-ein)))
        sum2=sum2-a1*xx(idumC2)*dexp(ein)/(1.0d0+dexp(ein))
     1           +a2*xx(idumC2)*dexp(-ein)/(1.0d0+dexp(-ein))        
300    continue
       if (idumC2 .eq. 1) then
        sum1=sum1-z1**2/(2.0d0*10.0d0)
        sum2=sum2-z1/10.0d0
       else
        sum1=sum1-z1**2/(2.0d0*100.0d0)
        sum2=sum2-z1/100.0d0
       endif
       ahz=sum1
       ahpz=sum2  
       return
       end       


       subroutine Ge1starnew ( iseed )
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n),tseq(5,n)
       real*8 e1star(n),e2star(n)
       real*8 case(6,n),x(np,n),xx(np+1)
       real*8 s(5,ngmax)
       integer ncase(6),idc(n),ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2),einv(5),ss(2),alogss,xgam
       integer iseed,iaccept(n)
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vtseq/tseq
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vtaus/taus
       common /vecbo/bo
       common /vecng/ng
       external  DRNUNF
       taus=1.0d0/tau
       sume1=0.0d0
       do 100 i=1,n
        if (case(6,i) .eq. 1.0d0) then 
         xx(1)=1.0d0
         do j=2,np+1
          xx(j)=x(j-1,i)     
         enddo
         ein=dotproduct(xx,alpha,np+1)
         pe1=1.0d0/(1.0d0+dexp(-ein))
c         
         do 101 k=1,2
          ts(i)=tseq(k,i)
          if(ts(i) .ne. 9999) then
           call alogsurv(k,i,ts(i),gam,alam,xgam,ss(k),ahaz)
          endif
101      continue
c
         sum1=dlog(1.0d0-pe1)+ss(1)
         sum2=dlog(pe1)-taus*dlog(1.0d0-ss(2))
         pp=dexp(sum2)/(dexp(sum1)+dexp(sum2))         
         call rnset( iseed )
          w=DRNUNF() !generate random number from U(0,1)
         call rnget( iseed )
         if (w .le. pp) then 
          e1star(i)=1.0d0 
         else
          e1star(i)=0.0d0
         endif
         sume1=sume1+e1star(i)
        endif
100    continue  
       return
       end

       subroutine Ge2starnew ( iseed )
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n),tseq(5,n)
       real*8 e1star(n),e2star(n)
       real*8 case(6,n),x(np,n),xx(np+1)
       real*8 s(5,ngmax)
       integer ncase(6),idc(n),ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2),einv(5),ssk(5),alogss,xgam
       integer iseed,iaccept(n)
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vtseq/tseq
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vtaus/taus
       common /vecbo/bo
       common /vecng/ng
       common /dummy2/idum2       
       external  DRNUNF
       taus=1.0d0/tau
       sume2=0.0d0
       do 100 i=1,n
        if (case(5,i) .eq. 1.0d0) then
         xx(1)=1.0d0
         do j=2,np+1
          xx(j)=x(j-1,i)
         enddo
         ein=dotproduct(xx,beta,np+1) 
         pe2=1.0d0/(1.0d0+dexp(-ein))
         do k=2,4
          ts(i)=tseq(k,i)
          if(ts(i) .eq. 9999) goto 100    
          call alogsurv(k,i,ts(i),gam,alam,xgam,alogss,ahaz)
          ssk(k)=-alogss
         enddo
         ss1=1.0d0+ssk(2)+ssk(4)
         ss2=1.0d0+ssk(2)+ssk(3)
         sum1=dlog(pe2)+dlog(ss1)
         sum2=dlog(1.0d0-pe2)+dlog(ss2)
         pp=dexp(sum1)/(dexp(sum1)+dexp(sum2))
c
         call rnset( iseed )
          w=DRNUNF() !generate random number from U(0,1)
         call rnget( iseed )
         if (w .le. pp) then 
          e2star(i)=1.0d0 
         else
          e2star(i)=0.0d0
         endif
         sume2=sume2+e2star(i)
        endif
100    continue  
       return
       end


	   
       subroutine Gtaun( iseed )
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)
       parameter(iprint=3,maxlag=20,nobs=50001)
       PARAMETER(maxtries=25,im=2,NS=10,in=1000)       
       real*8 emax,step
C ------------------------ Working variables ------------------------------
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1         xlb,xub,rwv(6*NS+15+in)
       real*8 start(1),xmin(1),step1(1)
       integer ifault,iwv(NS+np+in),ifcount
       integer iseed,idum2       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n)
       real*8 d(n),u(n),v(n),e1star(n),e2star(n),tseq(5,n)
       real*8 case(6,n)
       integer ncase(6)
       real*8 x(np,n),xbar(np),xbarL(np),xbar2(np),xsd(np)
       real*8 tt1(5,n),tt(5,n),s(5,ngmax),p(ngmax)
       integer idc(n),ntt(5),ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2)
       real*8 shape(2,ngmax),scale(2,ngmax),rg(1)
       real*8 ynewlo,reqmin
       real*8 xiold(1),xinew(1),ximean(1)
       real*8 tol,Rsig(2,2),RV(1,2)
       real*8 sig2inv,sig2inv2,ftaudev2
       integer imean,iseopt       
       integer iaccept(n)
       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vtseq/tseq
       common /vecd/d
       common /vecu/u
       common /vecv/v
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vtaus/taus
       common /vecbo/bo
       common  /vecng/ng
       common  /dummy2/idum2
       common /dummy3/idum3
       common  /dummy1C/idumC1
       common  /dummy2C/idumC2
c      external  DRNUNF
       external ftaun
       emax=60.0d0
       tau=1.0d0/taus
         step=0.1d0
         nopt = 1
         reqmin = 1.0d-10
         konvge = 5
         kcount = 1000
         step1(1)=0.2d0
         start(1)=tau
         call nelmin(ftaun, nopt, start, xmin, ynewlo, reqmin, step1,
     1         konvge, kcount, icount, numres, ifault1) 
        ximean(1)=xmin(1)
        write(*,*) 'ximean(1)',ximean(1)		
c.......calculate Hessian matrix
         call taudev2(ximean(1),ftaudev2)
c         xt=0.52d0
c         call taudev2(xt,testdev2)
         sig2inv=-ftaudev2

         ein=0.0d0
         ein=ein+(xiold(1)-ximean(1))*sig2inv
     1            *(xiold(1)-ximean(1))/2.0d0

         sumold=-ftaun(xiold)+ein
c         write(*,*) 'sumold',sumold   
c         write(*,*) 'sig2inv',sig2inv               
c.......perform metropolis
         do imetr=1,10
          call rnset( iseed )
           nr=1
           std=1/dsqrt(sig2inv)
           rn=DRNNOF() 
          call rnget( iseed )

          xinew(1)=ximean(1)+rn*std

          ein=0.0d0
          ein=ein+(xinew(1)-ximean(1))*sig2inv
     1             *(xinew(1)-ximean(1))/2.0d0

          sumnew=-ftaun(xinew)+ein
c          write(*,*) 'sumnew',sumnew
          ratio=sumnew-sumold
c          write(*,*) 'ratio',ratio
          if (ratio .ge. 0.0d0) then
           xiold=xinew
           tau=xinew(1)
           sumold=sumnew
c          iaccept=iaccept+1
          else
           call rnset( iseed )
            w=DRNUNF()
           call rnget( iseed )
           if (dlog(w) .le. ratio) then
            xiold=xinew
            sumold=sumnew
c           iaccept=iaccept+1
            tau=xinew(1)
           endif
          endif
         enddo  
       return
       end 
                            
       real*8 function ftaun(z1)
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)
       parameter(iprint=3,maxlag=20,nobs=50001)
       PARAMETER(maxtries=25,im=2,NS=10,in=1000)       
       real*8 emax,step
C ------------------------ Working variables ------------------------------
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1         xlb,xub,rwv(6*NS+15+in)
       real*8 start(1),xmin(1),step1(1)
       integer ifault,iwv(NS+np+in),ifcount
       integer iseed,idum2       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n)
       real*8 d(n),u(n),v(n),e1star(n),e2star(n),tseq(5,n)
       real*8 case(6,n),ss(5)
       integer ncase(6)
       real*8 x(np,n),xbar(np),xbarL(np),xbar2(np),xsd(np)
       real*8 tt1(5,n),tt(5,n),s(5,ngmax),p(ngmax)
       integer idc(n),ntt(5),ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2)
       real*8 shape(2,ngmax),scale(2,ngmax),rg(1)       
       integer imean,iseopt       
       integer iaccept(n)
       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vtseq/tseq
       common /vecd/d
       common /vecu/u
       common /vecv/v
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vtaus/taus
       common /vecbo/bo
       common  /vecng/ng
       common  /dummy2/idum2
       common /dummy3/idum3
       common  /dummy1C/idumC1
       common  /dummy2C/idumC2
c       external  DRNUNF
       tau=z1
c......calculate log likelihood and its derivative
       atau=0.01d0
       btau=0.01d0
       sum1=0.0d0
       sum2=0.0d0
       do 300 i=1,n
        do k=2,5
         ts(i)=tseq(k,i)     
c
         if(ts(i) .ne. 9999) then   
          call alogsurv(k,i,ts(i),gam,alam,xgam,alogss,ahaz)
          ss(k)=-alogss
         endif
        enddo
        scum1=dlog(1.0d0+ss(2))
        scum2=dlog(1.0d0+ss(2)+ss(3))
        scum3=dlog(1.0d0+ss(2)+ss(4))
        scum4=dlog(1.0d0+ss(2)+ss(4)+ss(5))
        a1=case(2,i)+case(3,i)+case(4,i)
        ss1a=a1*dlog(1.0d0+tau)
        ss1b=case(3,i)*dlog(1.0d0+2.0d0*tau)
        a3=- 2.0d0*case(2,i) - 3.0d0*case(3,i)
     1   - 2.0d0*case(4,i)-case(5,i)
        ss1c=a3*dlog(tau+0.00001d0)
        a4=-scum1*e1star(i)*case(6,i)
        ss1d=a4/(tau+0.00001d0)
        a5=-(1-e2star(i))*case(5,i)*scum2
     1         -e2star(i)*case(5,i)*scum3
        ss1e=a5*(1.0d0+1.0d0/(tau+0.00001d0))
        a6=-case(2,i)*scum2-case(4,i)*scum4
        ss1f=a6*(2.0d0+1.0d0/(tau+0.00001d0))
        a7=-case(3,i)*scum4
        ss1g=a7*(3.0d0+1.0d0/(tau+0.00001d0))
        sum1=sum1+ss1a+ss1b+ss1c+ss1d+ss1e+ss1f+ss1g      
 300   continue
       sakj=0.0d0
       sbkj=0.0d0
       do k=2,5
        do j=1,ng(k)
         sakj=sakj+ 0.001d0
         sbkj=sbkj+(alam(k,j)*0.001d0)
        enddo
       enddo
       sum1=sum1-(atau+1.0d0+sakj)*dlog(tau+0.00001d0)
     1    -(btau+sbkj)/(tau+0.00001d0)
       ftaun=-sum1 
       return
       end

              
       subroutine taudev2(z1,ftaudev2)
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)
       parameter(iprint=3,maxlag=20,nobs=50001)
       PARAMETER(maxtries=25,im=2,NS=10,in=1000)       
       real*8 emax,step
C ------------------------ Working variables ------------------------------
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1         xlb,xub,rwv(6*NS+15+in)
       real*8 start(1),xmin(1),step1(1)
       integer ifault,iwv(NS+np+in),ifcount
       integer iseed,idum2       
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),ts(n)
       real*8 d(n),u(n),v(n),e1star(n),e2star(n),tseq(5,n)
       real*8 case(6,n),ss(5)
       integer ncase(6)
       real*8 x(np,n),xbar(np),xbarL(np),xbar2(np),xsd(np)
       real*8 tt1(5,n),tt(5,n),s(5,ngmax),p(ngmax)
       integer idc(n),ntt(5),ng(5)
c      model parameters
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus
       real*8 bo(5,2)
       real*8 shape(2,ngmax),scale(2,ngmax),rg(1)       
       integer imean,iseopt       
       integer iaccept(n)
       external dacf,DRNUNF,DRNNOF,DLINDS,DNORDF
       common /vect1/t1
       common /vect2/t2
       common /vect3/t3
       common /vect4/t4
       common /vect0/t0
       common /vtseq/tseq
       common /vecd/d
       common /vecu/u
       common /vecv/v
       common /vcase/case
       common /vecx/x
       common /ve1/e1star
       common /ve2/e2star
       common /vecs/s
       common /valam/alam                  
       common /vgam/gam
       common /vbeta/beta
       common /valpha/alpha
       common /vomega/omega
       common /vtau/tau
       common /vtaus/taus
       common /vecbo/bo
       common  /vecng/ng
       common  /dummy2/idum2
       common /dummy3/idum3
       common  /dummy1C/idumC1
       common  /dummy2C/idumC2
c       external  DRNUNF
       tau=z1
c......calculate log likelihood and its derivative
       atau=0.1d0
       btau=0.1d0
       sum1=0.0d0
       sum2=0.0d0
       do 300 i=1,n
        do k=2,5
         ts(i)=tseq(k,i)
         if(ts(i) .ne. 9999) then   
          call alogsurv(k,i,ts(i),gam,alam,xgam,alogss,ahaz)
          ss(k)=-alogss
         endif 
        enddo
        scum1=dlog(1.0d0+ss(2))
        scum2=dlog(1.0d0+ss(2)+ss(3))
        scum3=dlog(1.0d0+ss(2)+ss(4))
        scum4=dlog(1.0d0+ss(2)+ss(4)+ss(5))
        a1=case(2,i)+case(3,i)+case(4,i)
        ss2a=a1/(1.0d0+tau)
        ss2b=2.0d0*case(3,i)/(1.0d0+2.0d0*tau)
        a3=- 2.0d0*case(2,i) - 3.0d0*case(3,i)
     1   - 2.0d0*case(4,i)-case(5,i)
        ss2c=a3/(tau+0.00001d0)
        a4=-scum1*e1star(i)*case(6,i)
        ss2d=-a4/(tau+0.00001d0)**2
        a5=-(1-e2star(i))*case(5,i)*scum2
     1         -e2star(i)*case(5,i)*scum3
        ss2e=-a5/(tau+0.00001d0)**2
        a6=-case(2,i)*scum2-case(4,i)*scum4
        ss2f=-a6/(tau+0.00001d0)**2
        a7=-case(3,i)*scum4
        ss2g=a7/(tau+0.00001d0)**2
        sum2=sum2+ss2a+ss2b+ss2c+ss2d+ss2e+ss2f+ss2g          
 300   continue
       sakj=0.0d0
       sbkj=0.0d0
       do k=2,5
        do j=1,ng(k)
         sakj=sakj+ 0.001d0
         sbkj=sbkj+(alam(k,j)*0.001d0)
        enddo
       enddo
       sum2=sum2-(atau+1.0d0+sakj)/(tau+0.00001d0)
     1    +(btau+sbkj)/(tau+0.00001d0)**2     
       ftaudev2=sum2
       return
       end

       subroutine alogsurv(idum1,idum2,t,gam,alam,
     1            axgam,aoutss,ahaz)
c      Input:
c        idumH1: k (k=1,2,3,4,5), where 
c                1 indicates t0 (lam(1,)), 
c                2 indicates t1 (lam(2,)),
c                3 indicates t2 (lam(3,)), 
c                4 indicates t3 (lam(4,)),
c                5 indicates t4 (lam(5,))
c        idumH2: i (i=1,2,...,n)
c             t: survival time (tseq(k,i))
c      Output:
c        axgam: linear combination of x'gamma for the ith subject
c               mathbf{x}_i'\boldsymbol{\gamma}_k
c        aouss: logarithm of survival for the ith subject
c               -H_{k0}(t) e^{\mathbf{x}_i'\boldsymbol{\gamma}_k}
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2)
       integer ng(5)
       real*8 x(np,n),s(5,ngmax),t
       real*8 axgam,aoutss
c      model parameters
       real*8 gam(5,np),alam(5,ngmax)
       common /vecx/x
       common /vecs/s
       common /vecng/ng
       k=idum1
       i=idum2
       ein=0.0d0
       do j=1,np
        ein=ein+x(j,i)*gam(k,j)
       enddo
       if (t .le. s(k,1)) jg=1
        do j1=1,ng(k)-1
         if ( (t .gt. s(k,j1)) .and.
     1         (t .le. s(k,j1+1)) ) then
          jg=j1+1
         endif
        enddo
        if (jg .eq. 1) then
         ein1=alam(k,1)*t
        else
         ein1=alam(k,1)*s(k,1)
         do j2=2,jg-1
          ein1=ein1+alam(k,j2)*(s(k,j2)-s(k,j2-1))
         enddo
         ein1=ein1+alam(k,jg)*(t-s(k,jg-1))
        endif
        axgam=ein
        aoutss=-ein1*dexp(ein)
		ahaz=alam(k,jg)
       return
       end

      real*8 function dotproduct(x,y,n)
C      Author: Md Tuhin Sheikh, Yale University      
c      Purpose: This function calculates the linear combination of two real vectors
c      Input:
c          x: a REAL vector of length n
c          y: a REAL vector of length n
c          n: length of x and y
c      Output:
c        x'y: linear combination of x and y
c        Note: THIS FUNCTION IS NOT DESIGNED TO CALCULATE FOR INTEGER VECTORS
       implicit real*8 (a-h,o-z)
       integer n
       real*8 x(n),y(n),result
       integer i
       result=0.0d0
       do i=1,n
        result=result+x(i)*y(i)
       enddo     
       dotproduct=result
       return
      end function dotproduct

