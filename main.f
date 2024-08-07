       program semicompukb
c      A Semicompeting Risks Model with an Application to UK Biobank 
c      Data to Identify Risk Factors for Diabetes Onset 
c      and Progression
       implicit real*8 (a-h,o-z)
       parameter(n=15000,ngmax=110,np=2,nsim=1,nt=20)       
       parameter(iprint=3,maxlag=20,nobs=50001)
       real*8 t1(n),t2(n),t3(n),t4(n),t0(n),tseq(5,n),ts(n)
       real*8 d(n),u(n),v(n),e1star(n),e2star(n)
       real*8 case(6,n),x(np,n),xx(np+1)
       integer ncase(6),idc(n),ng(5)
       real*8 xold(np,n),xbar(np),xbarL(np),xbar2(np),xsd(np)
       real*8 tt1(5,n),tt(5,n),s(5,ngmax),p(ngmax)
       integer ntt(5)
       real*8 alhaz(5),alss(5),axgamv(5),sdic,adic
c......model parameters 
       real*8 gam(5,np),alpha(np+1),beta(np+1),omega(n),tau
       real*8 alam(5,ngmax),taus,tauold,alamold(5,ngmax)
       real*8 gamold(5,np),alphaold(np+1),betaold(np+1)       
       real*8 bo(5,2),rg(1)
c......posterior estimates
       real*8 egam(5,np),segam(5,np)
       real*8 ealam(5,ngmax),sealam(5,ngmax)
       real*8 ealpha(np+1),sealpha(np+1),ebeta(np+1),sebeta(np+1)
       real*8 etau,setau
       real*8 ee1star(n),ee2star(n),eomega(n)
c......cumulative sums
       real*8 sgam(5,np),s2gam(5,np)
       real*8 salam(5,ngmax),s2alam(5,ngmax)
       real*8 salpha(np+1),s2alpha(np+1),sbeta(np+1),s2beta(np+1) 
       real*8 stau,s2tau,somega(n),se1star(n),se2star(n)
c......for HPD intervals
       real*8 gamupp(5,np),gamlow(5,np)
       real*8 alamupp(5,ngmax),alamlow(5,ngmax)
       real*8 alphaupp(np+1),alphalow(np+1)
       real*8 betaupp(np+1),betalow(np+1),taulow,tauupp
       real*8 alow(2),aupp(2)
c......sequence for hpd computation 
       real*8 seqgam(5,np,50001),seqlam(5,ngmax,50001)
       real*8 seqalpha(np+1,50001),seqbeta(np+1,50001)
       real*8 seqtau(50001)
c      for computing autocorrelation of the MCMC sequence       
       real*8 ac(maxlag+1),acv(maxlag+1),seac(maxlag+1),
     1         xacmean
       real*8 ahpd(50001)
c      for saving computation time
       integer now(3),nowb(3),nowe(3),ntoday(3),ntodayb(3)
       real etime,elapsed(2)
       integer imean,iseopt       
       integer iseed,iaccept(n)
c      functions which are called using IMSL library        
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
c
c......starting time       
       call idate(ntoday)
       ntodayb(1)=ntoday(1)
       ntodayb(2)=ntoday(2)
       ntodayb(3)=ntoday(3)
       call itime(now)
       nowb(1)=now(1)
       nowb(2)=now(2)
       nowb(3)=now(3)

c......number of pieces for baseline hazard functions
       ng(1)=1
       ng(2)=1
       ng(3)=1
       ng(4)=1
       ng(5)=1

c......prior for \lambda       
       do k=1,5
        do j=1,2
         bo(k,j)=0.1d0
        enddo
       enddo

c......read surv data
c      "id", "t1", "t2", "t3", "t4", "t0", "di", "ui",
c      "vi", "case1", "case2", "case3", "case4", "case5", "case6",
c      "x1", "x2"

       open(unit=16,file='exampledat.txt',status='old')
       do isim=1,nsim
        do i=1,n
         read(16,*) idc(i),t1(i),
     1     t2(i),t3(i),t4(i),
     1     t0(i),d(i),u(i),v(i),
     1     (case(j1,i),j1=1,6),(xold(j,i),j=1,np)
        enddo
       enddo
       close(16)


c......initialize cumulative variables to store simulated results
c......from posterior estimation
    
       sdic=0.0d0
       adic=0.0d0
c......initialize random seed
       iseed=9999999

c.....gamma      
     
        do i=1,n
         iaccept(i)=0
        enddo


c.......summary of samples under each case
        do k=1,6
         ncase(k)=0
        enddo
        do k=1,6
         do i=1,n
          ncase(k)=ncase(k)+case(k,i)
         enddo
        enddo

c.......summary of event times
        sum1=0.0d0
        sum2=0.0d0
c........t1      
        nte1=0
        do 175 i=1,n
         if (t1(i) .eq. 9999) goto 175
         nte1=nte1+1 
         sum1=sum1+t1(i)
         sum2=sum2+t1(i)**2
175     continue 
        barte1=sum1/real(nte1)
        sdte1=dsqrt((sum2-real(nte1)*barte1**2)/real(nte1-1))

c.......t2     
        sum1=0.0d0
        sum2=0.0d0 
        ntg1=0
        do 176 i=1,n
         if (t2(i) .eq. 9999) goto 176
         ntg1=ntg1+1 
         sum1=sum1+t2(i)
         sum2=sum2+t2(i)**2
176     continue 
        bartg1=sum1/real(ntg1)
        sdtg1=dsqrt((sum2-real(ntg1)*bartg1**2)/real(ntg1-1))

c.......t3
        sum1=0.0d0
        sum2=0.0d0 
        ntg2=0
        do 177 i=1,n
         if (t3(i) .eq. 9999) goto 177
         ntg2=ntg2+1 
         sum1=sum1+t3(i)
         sum2=sum2+t3(i)**2
177     continue 
        bartg2=sum1/real(ntg2)
        sdtg2=dsqrt((sum2-real(ntg2)*bartg2**2)/real(ntg2-1))       
 
c.......t4   
        sum1=0.0d0
        sum2=0.0d0    
        ntg3=0
        do 178 i=1,n
         if (t4(i) .eq. 9999) goto 178
         ntg3=ntg3+1 
         sum1=sum1+t4(i)
         sum2=sum2+t4(i)**2
178     continue 
        bartg3=sum1/real(ntg3)
        sdtg3=dsqrt((sum2-real(ntg3)*bartg3**2)/real(ntg3-1))       
 
c.......t0       
        sum1=0.0d0
        sum2=0.0d0
        ntd=0
        do 179 i=1,n
         if (t0(i) .eq. 9999) goto 179
         ntd=ntd+1 
         sum1=sum1+t0(i)
         sum2=sum2+t0(i)**2
179     continue 
        bartd=sum1/real(ntd)
        sdtd=dsqrt((sum2-real(ntd)*bartd**2)/real(ntd-1))
 
c.......standardizing covariates
        n1=0
        do j=1,np
         xbar(j)=0.0d0
         xbar2(j)=0.0d0
        enddo
 
        do i=1,n
         n1=n1+1
         do j=1,np
          xbar(j)=xbar(j)+xold(j,i)
          xbar2(j)=xbar2(j)+xold(j,i)**2
         enddo       
        enddo
 
        do j=1,np
         xbar(j)=xbar(j)/real(n1)
         xsd(j)=dsqrt((xbar2(j)-real(n1)*xbar(j)**2)/real(n1-1))
        enddo
       
        do i=1,n
         do j=1,np
          x(j,i)=(xold(j,i)-xbar(j))/xsd(j)
         enddo
        enddo 
 
        do k=1,5
         ntt(k)=0
        enddo
        icounto=0
        icount1=0
        icount2=0
        icount3=0
        icount4=0
 
        do i=1,n
         if(v(i) .eq. 1.0d0) then
          if(t0(i) .ne. 9999) then
          icounto=icounto+1
          tt1(1,icounto)=t0(i)
          endif
         endif
 
         if(d(i) .eq. 1.0d0) then
          if(t1(i) .ne. 9999) then
           icount1=icount1+1
           tt1(2,icount1)=t1(i)
          endif
         endif
 
         if((d(i) .eq. 1.0d0) .and. (v(i) .eq. 1.0d0)) then
          if(t2(i) .ne. 9999) then
           icount2=icount2+1
           tt1(3,icount2)=t2(i)
          endif
         endif
       
         if((d(i) .eq. 1.0d0) .and. (u(i) .eq. 1.0d0)) then
          if(t3(i) .ne. 9999) then
           icount3=icount3+1
           tt1(4,icount3)=t3(i)
          endif
         endif
 
         if((d(i) .eq. 1.0d0) .and. (u(i) .eq. 1.0d0)
     1      .and. (v(i) .eq. 1.0d0)) then
          if(t4(i) .ne. 9999) then
           icount4=icount4+1
           tt1(5,icount4)=t4(i)
          endif
         endif
        enddo
  
        write(*,*) 'icounto=',icounto
        write(*,*) 'icount1=',icount1
        write(*,*) 'icount2=',icount2
        write(*,*) 'icount3=',icount3
        write(*,*) 'icount4=',icount4
 
        ntt(1)=icounto
        ntt(2)=icount1
        ntt(3)=icount2
        ntt(4)=icount3
        ntt(5)=icount4
 
c......set initial values of s
c      sorting the survival times
        do k=1,5
         do i1=1,ntt(k)-1
          do i2=i1+1,ntt(k)
           if (tt1(k,i1) .gt. tt1(k,i2)) then
            einn=tt1(k,i1)
            tt1(k,i1)=tt1(k,i2)
            tt1(k,i2)=einn
           endif
          enddo
         enddo
        enddo

c.....Getting distinct failure times
        do k=1,5
         tt(k,1)=tt1(k,1)
        enddo 
        do k=1,5
         icount=1
         do i1=2,ntt(k)
          if (tt1(k,i1) .gt. tt(k,icount)) then
           tt(k,icount+1)=tt1(k,i1)
           icount=icount+1
          endif
         enddo
         ntt(k)=icount
        enddo
        write(*,*) 'ntt',ntt

        do k1=1,5
         ang=real(ng(k1))
         k=int(dlog(ang)/dlog(2.0d0))
         mm=ng(k1)-2**k
         kk=2**k-1
         do i=1,kk
          p(i)=real(i)/real(kk+1)
         enddo
         do j=1,mm
          p(kk+j)=real(2*j-1)/real(2**(k+1))
         enddo
         do i=1,kk+mm-1
          do j=i+1,kk+mm
           if (p(i) .gt. p(j)) then
            temp=p(i)
            p(i)=p(j)
            p(j)=temp
           endif
          enddo
         enddo
         do i=1,kk+mm
          write(*,*) 'j=',i,' p=',p(i)
         enddo
c.....compute cut points for t0
         do j=1,ng(k1)-1
          antt=real(ntt(k1))
          k=int(antt*p(j))
          einn=antt*p(j)
          pp=einn-real(k)
          if (pp .eq. 0.0d0) then
           s(k1,j)=(tt(k1,k)+tt(k1,k+1))/2.0d0
          else
           s(k1,j)=tt(k1,k+1)
          endif
          write(*,*) 'k1=',k1,' ntt=',ntt(k1),' j=',j,' s=',s(k1,j)
         enddo
         s(k1,ng(k1))=10000.0d0
        enddo

        do k=1,5
         do j2=1,ng(k)
          write(*,*) 's(k,j2)=',s(k,j2)
         enddo          
        enddo

c      storing survival time in a sequence  
        do k=1,5
         if(k .eq. 1) then
          do i=1,n
           tseq(k,i)=t0(i)
          enddo
         else if(k .eq. 2) then
          do i=1,n
           tseq(k,i)=t1(i)
          enddo
         else if(k .eq. 3) then
          do i=1,n
           tseq(k,i)=t2(i)
          enddo
         else if(k .eq. 4) then
          do i=1,n
           tseq(k,i)=t3(i)
          enddo
         else if(k .eq. 5) then
          do i=1,n
           tseq(k,i)=t4(i)
          enddo
         endif
        enddo


c......initialize parameters
c......gamma
       gam(1,1)=-1.0d0; gam(1,2)=0.20d0
       gam(2,1)=-1.00d0; gam(2,2)=0.20d0
       gam(3,1)=-1.00d0; gam(3,2)=0.20d0
       gam(4,1)=-1.00d0; gam(4,2)=0.20d0
       gam(5,1)=-1.00d0; gam(5,2)=0.20d0
       alpha(1)=-1.0d0; alpha(2)=-1.0d0; alpha(3)=0.50d0
       beta(1)=0.5d0; beta(2)=-1.0d0; beta(3)=0.50d0	   
c......lambda       
        do k=1,5
         do j=1,ng(k)
          alam(k,j)=0.50d0
         enddo
        enddo

        tau=0.65d0
        taus=1.0d0/tau

c......initializing e1star and e2star for cases 5 and 6
       call Ge1starnew( iseed )
       call Ge2starnew( iseed )
       call Gomega( iseed )        

c......set number of burn-in, replication, and thinning iteration       
       nbin=0!2000
       nrep=200!10000
       nthin=1!5

c......warm up gibbs

        call rnset( iseed )
        icount=1
        do 51 i=1,nbin
         do ithin=1,nthin
          call gibbs( iseed )
         enddo
         write(*,*) 'warm up i=',i, ' (of ',nbin,')'
         write(*,*) 'alam'
         do k = 1,5
         enddo
         write(*,980) (gam(1,j),j=1,np)
980      format('gamo=',5f12.4)        
         write(*,981) (gam(2,j),j=1,np)
981      format('gam1=',5f12.4)   
         write(*,982) (gam(3,j),j=1,np)
982      format('gam2=',5f12.4) 
         write(*,983) (gam(4,j),j=1,np)
983      format('gam3=',5f12.4) 
         write(*,984) (gam(5,j),j=1,np)
984      format('gam4=',5f12.4)      
         write(*,985) (alpha(j),j=1,np+1)
985      format('alpha=',5f12.4) 
         write(*,986) (beta(j),j=1,np+1)
986      format('beta=',5f12.4) 
         WRITE(*, '(A, F10.4)') 'tau=',tau
51      continue
 
c......initialize cumulative variables
c.......gamma
        do k=1,5
         do j=1,np
          sgam(k,j)=0.0d0
          s2gam(k,j)=0.0d0
         enddo
        enddo
c.......alpha, beta
        do j=1,np+1
         salpha(j)=0.0d0
         s2alpha(j)=0.0d0
         sbeta(j)=0.0d0
         s2beta(j)=0.0d0
        enddo
c.......lambda
        do k=1,5
         do j=1,ng(k)
          salam(k,j)=0.0d0
          s2alam(k,j)=0.0d0
         enddo
        enddo  
c.......tau    
        stau=0.0d0
        s2tau=0.0d0
c......omega
       do i=1,n
        somega(i)=0.0d0
        se1star(i)=0.0d0
        se2star(i)=0.0d0
       enddo

		sdev=0.0d0

c      Gibbs sampling starts after completeing burn-in iterations
        do 111 i=1,nrep
         do ithin=1,nthin
          call gibbs( iseed )
         enddo
          write(*,*) 'gibbs rep i=',i, ' (of ',nrep,')'         
          write(*,*) 'alam='
         do k = 1,5
          write(*,'(I4,5f10.4)') k-1,(alam(k,j),j=1,1) !just printing one
         enddo
         write(*,980) (gam(1,j),j=1,np)
         write(*,981) (gam(2,j),j=1,np)
         write(*,982) (gam(3,j),j=1,np)
         write(*,983) (gam(4,j),j=1,np)
         write(*,984) (gam(5,j),j=1,np)
         write(*,985) (alpha(j),j=1,np+1)
         write(*,986) (beta(j),j=1,np+1)     
         WRITE(*, '(A, F10.4)') 'tau=',tau
         write(*,989) omega(1),omega((n+1)/2),omega(n)
989      format('1st, mid, and last omega=',3f12.4)        
 
c.......cumulative sums and sequence for gamma       
        do k=1,5
         do j=1,np
          sgam(k,j)=sgam(k,j)+gam(k,j)
          s2gam(k,j)=s2gam(k,j)+gam(k,j)**2
          seqgam(k,j,i)=gam(k,j)
         enddo
        enddo
 
c.......cumulative sums and sequence for alpha and beta       
        do j=1,np+1
         salpha(j)=salpha(j)+alpha(j)
         s2alpha(j)=s2alpha(j)+alpha(j)**2
         seqalpha(j,i)=alpha(j)
         sbeta(j)=sbeta(j)+beta(j)
         s2beta(j)=s2beta(j)+beta(j)**2
         seqbeta(j,i)=beta(j)
        enddo
c.......cumulative sums and sequence for lambda
        do k=1,5
         do j=1,ng(k)
          salam(k,j)=salam(k,j)+alam(k,j)
          s2alam(k,j)=s2alam(k,j)+alam(k,j)**2
          seqlam(k,j,i)=alam(k,j)
         enddo
        enddo   
c.......cumulative sums and sequence for tau
        stau=stau+tau
        s2tau=s2tau+tau**2
        seqtau(i)=tau
c

       do 219 i1=1,n
        if(case(1,i1) .eq. 1.0) goto 219
        somega(i1)=somega(i1)+omega(i1)
219    enddo   


111     continue       
        close(300)



c......obtain posterior estimates
c      posterior estimates of gamma
       do k=1,5
        do j=1,np
         egam(k,j)=sgam(k,j)/real(nrep)
         segam(k,j)=dsqrt((s2gam(k,j)-real(nrep)*egam(k,j)**2)
     1       / real(nrep-1) )
        enddo
       enddo

c      posterior estimates of alpha
       do j=1,np+1
        ealpha(j)=salpha(j)/real(nrep)
        sealpha(j)=dsqrt( (s2alpha(j)-real(nrep)*ealpha(j)**2)
     1        / real(nrep-1) )      
       enddo

c      posterior estimates of beta       
       do j=1,np+1
        ebeta(j)=sbeta(j)/real(nrep)
        sebeta(j)=dsqrt( (s2beta(j)-real(nrep)*ebeta(j)**2)
     1       / real(nrep-1) )     
       enddo

c      posterior estimates of lambda       
       do k=1,5
        do j=1,ng(k)
         ealam(k,j)=salam(k,j)/real(nrep)
         sealam(k,j)=dsqrt( (s2alam(k,j)-real(nrep)*ealam(k,j)**2)
     1        / real(nrep-1) )
        enddo
       enddo

c      posterior estimates of tau       
       etau=stau/real(nrep)
       setau=dsqrt((s2tau-real(nrep)*etau**2)/ real(nrep-1) )


c......95% HPD intervals and ACF for gamma
        imean=1
        iseopt=1
        alphahpd=0.05d0
        do k=1,5
         do j=1,np
          do i=1,nrep
           ahpd(i)=seqgam(k,j,i)
          enddo
          call hpd(nrep,alphahpd,ahpd,alow,aupp)
          gamupp(k,j)=aupp(1)
          gamlow(k,j)=alow(1)
         enddo
        enddo

c......95% HPD intervals and ACF for alpha
        imean=1
        iseopt=1
        alphahpd=0.05d0        
        do j=1,np+1
         do i=1,nrep
          ahpd(i)=seqalpha(j,i)
         enddo
         call hpd(nrep,alphahpd,ahpd,alow,aupp)
         alphaupp(j)=aupp(1)
         alphalow(j)=alow(1)        
        enddo

c......95% HPD intervals and ACF for beta
        imean=1
        iseopt=1
        alphahpd=0.05d0        
        do j=1,np+1
         do i=1,nrep
          ahpd(i)=seqbeta(j,i)
         enddo
         call hpd(nrep,alphahpd,ahpd,alow,aupp)
         betaupp(j)=aupp(1)
         betalow(j)=alow(1)
        enddo

c......95% HPD intervals and ACF for tau
       imean=1
       iseopt=1
       alphahpd=0.05d0
       do i=1,nrep
        ahpd(i)=seqtau(i)
       enddo
       call hpd(nrep,alphahpd,ahpd,alow,aupp)
       tauupp=aupp(1)
       taulow=alow(1)

c......95% HPD intervals and ACF for lambda
       imean=1
       iseopt=1
       alphahpd=0.05d0
        do k=1,5
         do j=1,ng(k) 
          do i=1,nrep
           ahpd(i)=seqlam(k,j,i)
          enddo
          call hpd(nrep,alphahpd,ahpd,alow,aupp)
          alamupp(k,j)=aupp(1)
          alamlow(k,j)=alow(1)      
         enddo
        enddo


       open(unit=12,file="output.out",
     1     access='sequential',status='unknown')
        write(12,*) 'Simulation study of semicompeting risks model'
        write(12,*) 'with two nonterminal and one terminal event'
        write(12,*) 'Sample size n=',n
        write(12,*) 'Number of covariates np=',np
        write(12,*) 'Number of Gibbs iterations nrep=',nrep
        write(12,*) 'Number of burn-in iterations, bin=',nbin
        write(12,*) 'Number of thinning, thin',nthin        
        write(12,*) 'Summary of sample size in six cases:'
c        write(12,*) '------------------------------------'
c        write(12,*) 'case   count   proportion'
c        do k=1,6
c         write(12,'(2I8,f10.3)') k,ncase(k),100*real(ncase(k))/real(n)
c        enddo

       write(12,*) '-------------------------------------'
       write(12,*) 'k j     egam         segam      95% HPD Int.'
       do j1=1,5
        do j2=1,np
         write(12,1501) j1-1,j2,egam(j1,j2),segam(j1,j2),
     1  gamlow(j1,j2),gamupp(j1,j2)
1501    format(1x,2I5,2f16.4,'(',f12.4,',',f12.4,')')
        enddo
       enddo

       write(12,*) '-------------------------------------'
       write(12,*) 'j     ealpha         sealpha         95% HPD Int.'
       do j=1,np+1
        write(12,1401) j,ealpha(j),sealpha(j),alphalow(j),alphaupp(j)
1401    format(1x,I4,2f16.4,'(',f12.4,',',f12.4,')')
       enddo

       write(12,*) '-------------------------------------'
       write(12,*) 'j     ebeta         sebeta         95% HPD Int.'
       do j=1,np+1
        write(12,1301) j,ebeta(j),sebeta(j),betalow(j),betaupp(j)
1301    format(1x,I4,2f16.4,'(',f12.4,',',f12.4,')')
       enddo

       write(12,*) '-------------------------------------'
       write(12,*) '     tau         setau         95% HPD Int.'
       write(12,1201) etau,setau,taulow,tauupp
1201   format(1x,2f16.4,'(',f12.4,',',f12.4,')')

       write(12,*) '------------------------------------------'
       write(12,*) 'k  j    ealam       sealam         95% HPD Int.'
        do j1=1,5
         do j2=1,ng(j1)
         write(12,1003) j1-1,j2,ealam(j1,j2),sealam(j1,j2),
     1    alamlow(j1,j2),alamupp(j1,j2)
1003     format(1x,2I5,1x,2f16.4,'(',f12.4,',',f12.4,')')
         enddo
        enddo


c      calculating computation time       
       call idate(ntoday)
       call itime(now)
       write(12,*)
     1  'date: month/day/year and time: hour, minute, and second'
       write(12,*)
     1  'Begining at date=',ntodayb(2),'/',ntodayb(1),'/',ntodayb(3)
       write(12,*) 'time at',nowb(1),':',nowb(2),':',nowb(3)
       write(12,*)
     1   'Ending at date=',ntoday(2),'/',ntoday(1),'/',ntoday(3)
       write(12,*) 'time at',now(1),':',now(2),':',now(3)
       total = etime(elapsed)
       write(12, *) 'elapsed time in minutes'
       write(12,2216) total/60.0d0,elapsed(1)/60.0d0,elapsed(2)/60.0d0
2216   format('end: total=',f12.4,' user=',f12.4,
     1         ' system=', f12.4)
       write(12,*) '--------------------------------------------------'

       close(12)

       stop
       end


       include 'optim1.f'
       include 'gilks2.f'
       include 'hpd.f'
       include 'gibbs.f'
