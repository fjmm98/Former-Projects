c gfortran randomnumber.f -c
c gfortran cuanticavol.f -O3 -c -o cuanticavol.o 
c gfortran cuanticavol.o randomnumber.o -O3 -o cuanticavol		

		program cuanticavol
		use randomnumber
		complex*16::a
		complex*16,parameter::i=(0,1)
		real*8::l,ko,s,pi,pder,pizq,kn,x1,x2,coef,NT,T,d1,i1
		integer*8::N,nciclos,j,k,w,nd,nmax,jprima,p1,p2,cont
		complex*16,ALLOCATABLE::alpha(:),beta(:),phi(:),xhi(:),gama(:)
		complex*16,ALLOCATABLE::phin(:)

c identifico nmax como número de experimentos, NT numero de veces que se detecta a la derecha, T coeficiente de transmision
c identifico p1=4N/5, p2=N/5

c introducimos la semilla para los numeros aleatorios x1,x2
		call dran_ini(1729)

c Hallamos el valor de pi que sera usado.
		pi=acos(-1.d0)
c Pedimos los valores que serán necesarios
		N=2000
		l=5.0
		nciclos=N/5
		nmax=1000
		nd=166
		p1=(4*N)/5
		p2=N/5

c con estos valores realizamos los calculos que serán necesarios
		ko=(2.d0*pi*nciclos)/N
		s=(1.d0)/(4.d0*ko*ko)

c definimos la dimension de los vectores complejos
		allocate(alpha(0:(N-1)),beta(0:(N-1)))
		allocate(phi(0:N),xhi(0:N),phin(0:N))
		allocate(gama(1:(N-1)))
c ponemos las condiciones iniciales
		phi(0)=0.d0
		phi(N)=0.d0
		xhi(0)=0.d0
		xhi(N)=0.d0	
		alpha(N-1)=0.d0
		beta(N-1)=0.d0
		NT=0
c calculamos las alphas y gamas que usaremos que son constantes durante todo el proceso
		do j=N-1,1,-1
			if (j.lt.ceiling(0.6d0*N).and.j.gt.floor(0.4d0*N)) then
			alpha(j-1)=(-1.d0)/(-2-l*ko*ko+alpha(j)+((2.d0*i)/s))
			else
			alpha(j-1)=(-1.d0)/(-2+alpha(j)+((2.d0*i)/s))
			endif
			gama(j)=(-1.d0)/alpha(j-1)
		enddo

c hallamos el valor de phi en cada jota inicial
		do j=1,N-1,1
			a=cos(ko*j)+i*sin(ko*j)
			phi(j)=a*exp((-8.d0*(4*j-N)**2)/(N*N))

		enddo

c guardamos el valor inicial de phin
		do j=0,N
		phin(j)=phi(j)
		enddo
		
		do jprima=1,nmax
		d1=0
		i1=0
				do j=0,N
				phi(j)=phin(j)
				enddo	
		cont=0

			do while(d1.eq.0.and.i1.eq.0)		

			do w=1,nd
c la evolucionamos nd ciclos para que evolucione lo suficiente
c hallamos el valor de beta
				do j=N-1,1,-1
				beta(j-1)=(1.d0/gama(j))*(((4*i*phi(j))/s)-beta(j))
				enddo
c hallamos el valor de xhi
				do j=0,N-2,1
				xhi(j+1)=alpha(j)*xhi(j)+beta(j)
				enddo
c Por ultimo hallamos el valor evolucionado de phi
				do j=1,N-1,1
				phi(j)=xhi(j)-phi(j)
				enddo
			enddo
c hallamos la probabilidad de que esté a la derecha

		pder=0.d0
		do k=p1,N
		pder=pder+abs(phi(k))**2.d0
		enddo
		x1=dran_u()
c realizamos el paso 4 del algoritmo		
		if (x1.lt.pder) then
		NT=NT+1
		d1=1.0
c		write(6,*) 'Se ha detectado a la derecha'
		cycle
		endif
c como no se ha detectado normalizamos	
	
		do k=p1,N
		phi(k)=0.d0
		enddo

c hallamos la constante de normalización

		kn=0.d0
		do k=1,N
		kn=kn+abs(phi(k))**2.d0
		enddo

c normalizamos la función de Onda		

		do k=1,N
		phi(k)=(phi(k)/(sqrt(kn)))
		enddo

c Hallamos la probabilidad de hallar a la izquierda

		pizq=0.d0
		do k=1,p2
		pizq=pizq+abs(phi(k))**2.d0
		enddo

		x2=dran_u()
		
		if(x2.lt.pizq) then
		i1=1.d0
c		write(6,*) 'Se ha detectado a la izquierda'
		cycle
		endif

c Volvemos a normalizar la función de onda	
		do k=1,p2
		phi(k)=0.d0
		enddo

		kn=0.d0

		do k=1,N
		kn=kn+abs(phi(k))**2.d0
		enddo

		do k=1,N
		phi(k)=(phi(k)/(sqrt(kn)))
		enddo

		cont=cont+1

		if(cont.eq.100) then
		i1=1.d0
c		write(6,*) 'Ha sobrepasado el limite de iteraciones '
		cycle
		endif
		
		enddo
		enddo

c una vez realizado todo el proceso obtenemos los datos
		 T=(dble(NT))/(dble(nmax))
		write(6,*)'El coeficiente de transmisión es: '
		write(6,*) T
		write(6,*)'El número de veces que se ha detectado ha sido: '
		write(6,*) NT

		end program
