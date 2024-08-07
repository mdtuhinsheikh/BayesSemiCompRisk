      SUBROUTINE initial(ns,m,emax,x,hx,hpx,lb,xlb,ub,xub,ifault,
     1             iwv,rwv,eval)
C                                                                       
C This subroutine takes as input the number of starting values m        
C  and the starting values x(i),hx(i),hpx(i)  i=1,m                     
C As output we have pointer ipt along with ilow and ihigh and the lower 
C and upper hulls defined  by z,hz,scum,cu,hulb,huub stored in working 
C vectors iwv and rwv                   
C Ifault detects wrong starting points or non-concavity                 
C                                                                       
      INTEGER ns,nn,m,ilow,ihigh,ifault,i,iwv(*)
      INTEGER iipt,iz,ihuz,iscum,ix,ihx,ihpx
      LOGICAL ub,lb,horiz
      DOUBLE PRECISION xlb,xub,emax,x(*),hx(*),hpx(*),rwv(*),expon
      DOUBLE PRECISION hulb,huub,eps,cu,alcu,huzmax,zlog,zmax
      DOUBLE PRECISION a,b
      external eval
      zlog(a)=dlog(a)
      zmax(a,b)=dmax1(a,b)
C
C DESCRIPTION OF PARAMETERS and place of storage
C
C     lb   iwv(5) : boolean indicating if there is a lower bound to the
C                    domain
C     ub   iwv(6) : boolean indicating if there is an upper bound
C     xlb  rwv(8) : value of the lower bound
C     xub  rwv(9) : value of the upper bound
C     emax rwv(3) : large value for which it is possible to compute
C                   an exponential, eps=exp(-emax) is taken as a small
C                   value used to test for numerical unstability
C     m    iwv(4) : number of starting points
C     ns   iwv(3) : maximum number of points defining the hulls
C     x    rwv(ix+1)  : vector containing the abscissae of the starting
C                       points
C     hx   rwv(ihx+1) : vector containing the ordinates 
C     hpx  rwv(ihpx+1): vector containing the derivatives
C     ifault      : diagnostic
C     iwv,rwv     : integer and real working vectors
C                                                                       
      eps=expon(-emax,emax)
      ifault=0                                                          
      ilow=1                                                            
      ihigh=1                                                           
      nn=ns+1
C at least one starting point
      IF (m.LT.1) ifault=1
      huzmax=hx(1)                                                      
      IF (.NOT.ub) xub=0.0
      IF (.NOT.lb) xlb=0.0
      hulb=(xlb-x(1))*hpx(1)+hx(1)                                      
      huub=(xub-x(1))*hpx(1)+hx(1)                                      
C if bounded on both sides                                              
      IF ((ub).AND.(lb)) THEN                               
        huzmax=zmax(huub,hulb)                                         
        horiz=(abs(hpx(1)).LT.eps)
        IF (horiz) THEN
          cu=expon((huub+hulb)*0.5-huzmax,emax)*(xub-xlb)
        ELSE
          cu=expon(huub-huzmax,emax)*(1-expon(hulb-huub,emax))/hpx(1)
        ENDIF
      ELSE IF ((ub).AND.(.NOT.lb))THEN                           
C if bounded on the right and unbounded on the left 
        huzmax=huub                                                     
        cu=1.0/hpx(1)
      ELSE IF ((.NOT.ub).AND.(lb))THEN                           
C if bounded on the left and unbounded on the right                     
        huzmax=hulb                                                     
        cu=-1.0/hpx(1)
C if unbounded at least 2 starting points
      ELSE
        cu=0.0
        IF (m.LT.2) ifault=1
      ENDIF                                                             
      IF (cu.GT.0.0) alcu=zlog(cu)
C set pointers 
      iipt=6
      iz=9
      ihuz=nn+iz
      iscum=nn+ihuz
      ix=nn+iscum
      ihx=nn+ix
      ihpx=nn+ihx
C store values in working vectors
      iwv(1)=ilow
      iwv(2)=ihigh
      iwv(3)=ns
      iwv(4)=1
      IF (lb) THEN
        iwv(5)=1
      ELSE
        iwv(5)=0
      ENDIF 
      IF (ub) THEN
        iwv(6)=1
      ELSE
        iwv(6)=0
      ENDIF 
      IF (ns.LT.m) ifault=2
      iwv(iipt+1)=0                                                     
      rwv(1)=hulb
      rwv(2)=huub
      rwv(3)=emax
      rwv(4)=eps
      rwv(5)=cu
      rwv(6)=alcu
      rwv(7)=huzmax
      rwv(8)=xlb
      rwv(9)=xub
      rwv(iscum+1)=1.0
      DO 9 i=1,m
        rwv(ix+i)=x(i)
        rwv(ihx+i)=hx(i)
        rwv(ihpx+i)=hpx(i)
  9   CONTINUE
C create lower and upper hulls                                          
      i=1                                                               
 10   IF (i.LT.m)THEN                                                   
         CALL update(iwv(4),iwv(1),iwv(2),iwv(iipt+1),rwv(iscum+1),
     +           rwv(5),rwv(ix+1),rwv(ihx+1),rwv(ihpx+1),rwv(iz+1),
     +           rwv(ihuz+1),rwv(7),rwv(3),lb,rwv(8),rwv(1),ub,rwv(9),
     +           rwv(2),ifault,rwv(4),rwv(6),eval) 
         i=iwv(4)
         IF (ifault.NE.0) RETURN
      GOTO 10                                                           
      ENDIF                                                             
C test for wrong starting points                                        
      IF ((.NOT.lb).AND.(hpx(iwv(1)).LT.eps)) ifault=3
      IF ((.NOT.ub).AND.(hpx(iwv(2)).GT.-eps)) ifault=4
      RETURN
      END
C                                                                       
C***********************************************************************
C                                                                       
      SUBROUTINE sample(iwv,rwv,eval,beta,ifault,iseed)
      DOUBLE PRECISION beta,rwv(*)
      INTEGER iipt,iz,ns,nn,ihuz,iscum,ix,ihx,ihpx,ifault,iwv(*)
      INTEGER iseed
      LOGICAL ub,lb
      external eval

C
C     set pointers
C
      iipt=6
      iz=9
      ns=iwv(3)
      nn=ns+1
      ihuz=nn+iz
      iscum=nn+ihuz
      ix=nn+iscum
      ihx=nn+ix
      ihpx=nn+ihx
      lb=.FALSE.
      ub=.FALSE.
      IF (iwv(5).EQ.1) lb=.TRUE.
      IF (iwv(6).EQ.1) ub=.TRUE.
C
C     call sampling subroutine
C
      CALL spl1(ns,iwv(4),iwv(1),iwv(2),iwv(iipt+1),rwv(iscum+1),
     +  rwv(5),rwv(ix+1),rwv(ihx+1),rwv(ihpx+1),rwv(iz+1),rwv(ihuz+1),  
     +  rwv(7),lb,rwv(8),rwv(1),ub,rwv(9),rwv(2),eval,beta,ifault,
     +  rwv(3),rwv(4),rwv(6),iseed)    
      RETURN
      END
C                                                                       
C***********************************************************************
C                                                                       
      SUBROUTINE spl1(ns,n,ilow,ihigh,ipt,scum,cu,x,hx,hpx,z,huz,    
     +  huzmax,lb,xlb,hulb,ub,xub,huub,eval,beta,ifault,emax,eps,alcu,
     +  iseed)  
C                                                                       
C this subroutine performs the adaptive rejection sampling, it calls 
C subroutine splhull to sample from the upper hull ,if the sampling 
C involves a function evaluation it calls the updating subroutine
C ifault is a diagnostic of any problem: non concavity, 0 random number 
C    or numerical imprecision
C                                                                       
      INTEGER ns,n,ilow,ihigh,ifault,i,j,ipt(*),n1,l,gencnt
      INTEGER iseed
      LOGICAL ub,lb,sampld
      DOUBLE PRECISION z(*),huz(*),x(*),hx(*),hpx(*),scum(*)
      DOUBLE PRECISION xlb,xub,emax,u1,u2,beta,alu1,fx
      DOUBLE PRECISION hulb,huub,eps,cu,alcu,huzmax,alhl,alhu
      DOUBLE PRECISION zlog,a
      DOUBLE PRECISION DRNUNF
      EXTERNAL eval,DRNUNF
      zlog(a)=dlog(a)
C                                                                       
      sampld=.FALSE.
      gencnt=0.0                                                    
10    IF (.NOT.sampld) THEN
        gencnt=gencnt+1
C        PRINT*, 'Number of generations =',gencnt
	  
        call rnset(iseed)
        u2=DRNUNF() 
        call rnget(iseed)

C test for zero random number                                         
        IF (u2.EQ.0.0) THEN                                             
          ifault=6
          RETURN
        ENDIF                                                           
      CALL splhull(u2,ipt,ilow,lb,xlb,hulb,huzmax,alcu,x,hx,hpx,
     +                              z,huz,scum,eps,emax,beta,i,j)
C      PRINT*, 'Exit from splhull'
C sample u1 to compute rejection                                        
	call rnset(iseed)
        u1=DRNUNF()                                                    
        call rnget(iseed)
        IF (u1.EQ.0.0) ifault=6                                         
        alu1=zlog(u1)                                                   
C compute alhu: upper hull at point u1                                  
        alhu=hpx(i)*(beta-x(i))+hx(i)-huzmax                          
        IF ((beta.GT.x(ilow)).AND.(beta.LT.x(ihigh))) THEN              
C compute alhl: value of the lower hull at point u1                     
          IF (beta.GT.x(i)) THEN                                        
            j=i
            i=ipt(i)
          ENDIF
          alhl=hx(i)+(beta-x(i))*(hx(i)-hx(j))/(x(i)-x(j))-huzmax 
C squeezing test                                                        
          IF ((alhl-alhu).GT.alu1) THEN                                 
             sampld=.TRUE.                                              
          ENDIF                                                         
        ENDIF                                                           
C if not sampled evaluate the function ,do the rejection test and update
        IF (.NOT.sampld) THEN                                           
          n1=n+1                                                        
          x(n1)=beta                                                    
          CALL eval(x(n1),hx(n1),hpx(n1))                               
          fx=hx(n1)-huzmax                                              
          IF (alu1.lt.(fx-alhu)) sampld=.TRUE.                          
C update while the number of points defining the hulls is lower than ns
          IF (n.LT.ns)  THEN
C            PRINT*, 'Entering update'
            CALL update(n,ilow,ihigh,ipt,scum,cu,x,hx,hpx,z,huz,huzmax, 
     +           emax,lb,xlb,hulb,ub,xub,huub,ifault,eps,alcu,eval) 
C            PRINT*, 'Exit from update'
          ENDIF
          IF (ifault.NE.0) RETURN                   
        ENDIF                                                           
        GOTO 10                                                         
      ENDIF            
C      PRINT*, 'Exit from spl1'                                             
      RETURN                                                            
      END                                                               
C                                                                       
C***********************************************************************
C                                                                       
      SUBROUTINE splhull(u2,ipt,ilow,lb,xlb,hulb,huzmax,alcu,x,hx,hpx,
     +                              z,huz,scum,eps,emax,beta,i,j)
C
C this subroutine samples beta from the normalised upper hull
C
      DOUBLE PRECISION z(*),huz(*),x(*),hx(*),hpx(*),scum(*)
      DOUBLE PRECISION expon,xlb,emax,eh,beta,u2,a,logdu,logtg,sign
      DOUBLE PRECISION hulb,eps,alcu,huzmax,zlog
      INTEGER ilow,i,j,ipt(*)
      LOGICAL lb,horiz
      zlog(a)=dlog(a)
C
        i=ilow                                                          
C
C find from which exponential piece you sample                        
 20   IF (u2.gt.scum(i))THEN                                          
        j=i                                                           
        i=ipt(i)                                                      
        GOTO 20                                                       
      ENDIF                                                           
      IF (i.eq.ilow) THEN                                             
C sample below z(ilow),depending on the existence of a lower bound  
        IF (lb) THEN                                            
          eh=hulb-huzmax-alcu                                         
          horiz=(abs(hpx(ilow)).LT.eps)
          IF (horiz) THEN
            beta=xlb+u2*expon(-eh,emax)
          ELSE
           sign=abs(hpx(i))/hpx(i)
           logtg=zlog(abs(hpx(i)))
           logdu=zlog(u2)
           eh=logdu+logtg-eh
           IF (eh.LT.emax) THEN
            beta=xlb+zlog(1.0+sign*expon(eh,emax))/hpx(i)
           ELSE
            beta=xlb+eh/hpx(i)
           ENDIF                                                       
          ENDIF
        ELSE                                                          
C     hpx(i) must be positive , x(ilow) is left of the mode
          beta=(zlog(hpx(i)*u2)+alcu-hx(i)+x(i)*hpx(i)+huzmax)/hpx(i)   
        ENDIF                                                         
      ELSE                                                            
C   sample above(j)                                                   
        eh=huz(j)-huzmax-alcu                                         
        horiz=(abs(hpx(i)).LT.eps)
        IF (horiz) THEN
         beta=z(j)+(u2-scum(j))*expon(-eh,emax)
        ELSE
         sign=abs(hpx(i))/hpx(i)
         logtg=zlog(abs(hpx(i)))
         logdu=zlog(u2-scum(j))
         eh=logdu+logtg-eh
         IF (eh.LT.emax) THEN
          beta=z(j)+(zlog(1.0+sign*expon(eh,emax)))/hpx(i)  
         ELSE                                                        
          beta=z(j)+eh/hpx(i)
         ENDIF
        ENDIF
      ENDIF                                                           
      RETURN
      END
C                                                                       
C***********************************************************************
      SUBROUTINE intersection(x1,y1,yp1,x2,y2,yp2,z1,hz1,eps,ifault)    
C                                                                       
C computes the intersection (z1,hz1) between 2 tangents defined by
C   x1,y1,yp1 and x2,y2,yp2
C
      DOUBLE PRECISION x1,y1,yp1,x2,y2,yp2,z1,hz1,eps,y12,y21,dh
      INTEGER ifault
C
C first test for non-concavity                                          
      y12=y1+yp1*(x2-x1)
      y21=y2+yp2*(x1-x2)
      IF ((y21.lt.y1).OR.(y12.LT.y2)) THEN                              
         ifault=5                                                       
         RETURN                                                         
      ENDIF                                                             
      dh=yp2-yp1                                                        
C  IF the lines are nearly parallel,
C  the intersection is taken at the midpoint
      IF (abs(dh).LE.eps)THEN                                           
        z1=0.5*(x1+x2)                                                  
        hz1=0.5*(y1+y2)                                                 
C  Else compute from the left or the right for greater numerical
C       precision
      ELSE IF (abs(yp1).LT.abs(yp2)) THEN                               
        z1=x2+(y1-y2+yp1*(x2-x1))/dh                                    
        hz1=yp1*(z1-x1)+y1                                              
      ELSE                                                              
        z1=x1+(y1-y2+yp2*(x2-x1))/dh                                    
        hz1=yp2*(z1-x2)+y2                                              
      ENDIF                                                             
C  test for misbehaviour due to numerical imprecision
      IF ((z1.LT.x1).OR.(z1.GT.x2)) ifault=7
      RETURN                                                            
      END                                                               
C                                                                       
C***********************************************************************
C                                                                       
      SUBROUTINE update(n,ilow,ihigh,ipt,scum,cu,x,hx,hpx,z,huz,        
     +  huzmax,emax,lb,xlb,hulb,ub,xub,huub,ifault,eps,alcu,eval)
C                                                                       
C this subroutine increments n and updates all the parameters which
C define the lower and the upper hull
C                                                                       
      INTEGER n,ilow,ihigh,ifault,i,j,ipt(*)
      LOGICAL ub,lb,horiz
      DOUBLE PRECISION xlb,xub,x(*),hx(*),hpx(*),z(*),huz(*),scum(*)
      DOUBLE PRECISION hulb,huub,eps,cu,alcu,huzmax,emax,expon,dh,u
      DOUBLE PRECISION a,b,zlog,zmax
      EXTERNAL eval
      zlog(a)=dlog(a)
      zmax(a,b)=dmax1(a,b)
C
C DESCRIPTION OF PARAMETERS and place of storage
C
C     ilow iwv(1)    : index of the smallest x(i)
C     ihigh iwv(2)   : index of the largest x(i)
C     n    iwv(4)    : number of points defining the hulls
C     ipt  iwv(iipt) : pointer array:  ipt(i) is the index of the x(.) 
C                      immediately larger than x(i)
C     hulb rwv(1)    : value of the upper hull at xlb
C     huub rwv(2)    : value of the upper hull at xub
C     cu   rwv(5)    : integral of the exponentiated upper hull divided
C                      by exp(huzmax)
C     alcu rwv(6)    : logarithm of cu
C     huzmax rwv(7)  : maximum of huz(i); i=1,n
C     z    rwv(iz+1) : z(i) is the abscissa of the intersection between
C                      the tangents at x(i) and x(ipt(i))
C     huz  rwv(ihuz+1): huz(i) is the ordinate of the intersection
C                        defined above
C     scum rwv(iscum): scum(i) is the cumulative probability of the 
C                      normalised exponential of the upper hull 
C                      calculated at z(i)
C     eps  rwv(4)    : =exp(-emax) a very small number 
C
      n=n+1                                                             
C update z,huz and ipt                                                  
      IF (x(n).LT.x(ilow)) THEN                                         
C insert x(n) below x(ilow)                                             
C   test for non-concavity                                              
        IF (hpx(ilow).GT.hpx(n)) ifault=5 
        ipt(n)=ilow
        CALL intersection(x(n),hx(n),hpx(n),x(ilow),hx(ilow),hpx(ilow), 
     +         z(n),huz(n),eps,ifault)                               
        IF (ifault.NE.0) RETURN                     
        IF (lb) hulb=hpx(n)*(xlb-x(n))+hx(n)                      
        ilow=n                                                          
      ELSE                                                              
        i=ilow                                                          
        j=i                                                             
C find where to insert x(n)                                             
  10    IF ((x(n).GE.x(i)).AND.(ipt(i).ne.0))THEN                       
          j=i                                                           
          i=ipt(i)                                                      
          GOTO 10                                                       
         ENDIF                                                          
        IF (x(n).GE.x(i)) THEN                                          
C insert above x(ihigh)                                                 
C   test for non-concavity                                              
          IF (hpx(i).LT.hpx(n)) ifault=5
          ihigh=n                                                       
          ipt(i)=n                                                      
          ipt(n)=0                                                      
          CALL intersection(x(i),hx(i),hpx(i),x(n),hx(n),hpx(n),        
     +         z(i),huz(i),eps,ifault)                               
          IF (ifault.NE.0) RETURN                   
          huub=hpx(n)*(xub-x(n))+hx(n)                                  
          z(n)=0.0                                                      
          huz(n)=0.0
        ELSE                                                            
C insert x(n) between x(j) and x(i)                                     
C   test for non-concavity                                              
          IF ((hpx(j).LT.hpx(n)).OR.(hpx(i).GT.hpx(n))) ifault=5
          ipt(j)=n                                                      
          ipt(n)=i                                                      
C     insert z(j) between x(j) and x(n)                                 
          CALL intersection(x(j),hx(j),hpx(j),x(n),hx(n),hpx(n),        
     +         z(j),huz(j),eps,ifault)                               
          IF (ifault.NE.0) RETURN                   
C     insert z(n) between x(n) and x(i)                                 
          CALL intersection(x(n),hx(n),hpx(n),x(i),hx(i),hpx(i),        
     +         z(n),huz(n),eps,ifault)                               
          IF (ifault.NE.0) RETURN                   
        ENDIF                                                           
      ENDIF                                                             
C update huzmax                                                         
      j=ilow                                                            
      i=ipt(j)                                                          
      huzmax=huz(j)                                                     
 20   IF ((huz(j).LT.huz(i)).AND.(ipt(i).NE.0))THEN                     
        j=i                                                             
        i=ipt(i)                                                        
        huzmax=zmax(huzmax,huz(j))                                     
        GOTO 20                                                         
      ENDIF                                                             
      IF (lb) huzmax=zmax(huzmax,hulb)
      IF (ub) huzmax=zmax(huzmax,huub)
C update cu                                                             
C  scum receives area below exponentiated upper hull left of z(i)       
      i=ilow                                                            
      horiz=(abs(hpx(ilow)).LT.eps)
      IF ((.NOT.lb).AND.(.NOT.horiz)) THEN
        cu=expon(huz(i)-huzmax,emax)/hpx(i)                             
      ELSE IF (lb.AND.horiz) THEN
        cu=(z(ilow)-xlb)*expon(hulb-huzmax,emax)
      ELSE IF (lb.AND.(.NOT.horiz)) THEN
        dh=hulb-huz(i)
        IF (dh.GT.emax) THEN
          cu=-expon(hulb-huzmax,emax)/hpx(i)
        ELSE
          cu=expon(huz(i)-huzmax,emax)*(1-expon(dh,emax))/hpx(i)
        ENDIF
      ELSE
        cu=0
      ENDIF                                                             
      scum(i)=cu                                                        
      j=i                                                               
      i=ipt(i)                                                          
 30   IF (ipt(i).NE.0)THEN                                              
        dh=huz(j)-huz(i)                                                
        horiz=(abs(hpx(i)).LT.eps)
        IF (horiz) THEN
          cu=cu+(z(i)-z(j))*expon((huz(i)+huz(j))*0.5-huzmax,emax)
        ELSE
          IF (dh.LT.emax) THEN                                          
            cu=cu+expon(huz(i)-huzmax,emax)*(1-expon(dh,emax))/hpx(i) 
          ELSE                                                          
            cu=cu-expon(huz(j)-huzmax,emax)/hpx(i)                      
          ENDIF                                                         
        ENDIF
        j=i                                                             
        i=ipt(i)                                                        
        scum(j)=cu                                                      
        GOTO 30                                                         
      ENDIF                                                             
      horiz=(abs(hpx(i)).LT.eps)
C if the derivative is very small the tangent is nearly horizontal
      IF (.NOT.(ub.OR.horiz)) THEN
        cu=cu-expon(huz(j)-huzmax,emax)/hpx(i)                          
      ELSE IF (ub.AND.horiz) THEN
        cu=cu+(xub-x(i))*expon((huub+hx(i))*0.5-huzmax,emax)
      ELSE IF (ub.AND.(.NOT.horiz)) THEN
        dh=huz(j)-huub
        IF (dh.GT.emax) THEN
          cu=cu-expon(huz(j)-huzmax,emax)/hpx(i)
        ELSE
          cu=cu+expon(huub-huzmax,emax)*(1-expon(dh,emax))/hpx(i)        
        ENDIF
      ENDIF                                                             
      scum(i)=cu                                                        
      IF (cu.GT.0) alcu=zlog(cu)
C normalize scum to obtain a cumulative probability while excluding     
C    unnecessary points                                                 
      i=ilow                                                            
      u=(cu-scum(i))/cu                                                 
      IF ((u.EQ.1.0).AND.(hpx(ipt(i)).GT.0)) THEN
        ilow=ipt(i)                                                     
        scum(i)=0.0                                                     
      ELSE                                                              
        scum(i)=1.0-u                                                   
      ENDIF                                                             
      j=i                                                               
      i=ipt(i)                                                          
 40   IF (ipt(i).NE.0) THEN                                             
        j=i                                                             
        i=ipt(i)                                                        
        u=(cu-scum(j))/cu                                               
        IF ((u.EQ.1.0).AND.(hpx(i).GT.0)) THEN 
          ilow=i                                                        
        ELSE                                                            
          scum(j)=1.0-u                                                 
        ENDIF                                                           
        GOTO 40                                                         
      ENDIF                                                             
      scum(i)=1.0
      IF (ub) huub=hpx(ihigh)*(xub-x(ihigh))+hx(ihigh)
      IF (lb) hulb=hpx(ilow)*(xlb-x(ilow))+hx(ilow)                     
C      do 100 ii=1,n
C      100   write(7,110) ii,z(ii),scum(ii)
C      110  format('i= ',i2,2x,2(2x,E10.4))
C       write(7,*) ' huzmax = ',huzmax
      RETURN                                                            
      END                                                               
C                                                                       
C***********************************************************************
C                                                                       
      DOUBLE PRECISION FUNCTION expon(x,emax)   
C                                                                       
C performs an exponential without underflow                             
C                                                                       
      DOUBLE PRECISION x,emax,zexp
      zexp(x)=dexp(x)
C
      IF (x.LT.-emax) THEN                                          
        expon=0.0d0                                                       
      ELSE                                                              
        expon=zexp(min(x,emax))       
      ENDIF                                                             
      RETURN                                                            
      END                                                               
