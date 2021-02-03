c descargar el archivo de la pagina web y usamos las subrutinas que aparecen en la web.
c inicializamos el numero aleatorio call DRAM(1729) solo se usa una vez!!!! y 1729 puede ser cualquier numero x=DRAMU()*7.d0 genera un numero aleatorio 
c gfortran randomnumber.f -c
c gfortran ising.f -O3 -c -o ising.o 
c gfortran ising.o randomnumber.o -O3 -o ising 
		program ising
		use randomnumber
		real*8,ALLOCATABLE::S(:,:)
		real*8::T,p,E,expo,uno
		real*8::mn,en,cn,fi
		integer*8::N,i,j,opcion,k,Maxiter,cont,m,v,contunos
c mn:Magnetizacion magnética, en: energía media, cn: calor específico,fi: función correlación
		cont=0
		k=2
c introducimos la dimension del sistema 
		call dran_ini(1729)
		
		write(6,*)'Introduce la dimensión del sistema: '
		read(5,*) N
		N=N-1
		ALLOCATE(S(0:N,0:N))
		write(6,*)'Introduce el número de pasos Montecarlo: '
		read(5,*) Maxiter

c Paso 0) iniciamos la matriz S prefijando los valores de cada celda en 1. (Forma ordenada)
		write(6,*)'Introduce la temperatura del sistema (0<T<5): '
		read(5,*)T

		write(6,*)'Como incializamos el sistema (Ordenado=1/No.Ord=2)'
		read(5,*)opcion
		
		if (opcion.eq.1) then
			do i=0,N
			do j=0,N
				S(i,j)=1
			enddo
			enddo
c la otra opcion es inicializarlo aleatoriamente con -1,1
		else 
			do i=0,N
			do j=0,N
			S(i,j)=i_dran(k)
			S(i,j)=2*S(i,j)-3
		enddo
		enddo
		endif
c Introducimos el Algoritmo
		open(UNIT=7,FILE="unos.dat",STATUS='OLD')
		open(UNIT=8,FILE="nounos.dat",STATUS='OLD')
		do v=1,Maxiter
			do m=1,(N+1)*(N+1)
c Escogemos un punto al azar de la matriz		
		i=i_dran(N+1)-1
		j=i_dran(N+1)-1
c Evaluamos p 	
		E=2*S(i,j)*(S(modulo(i+1,N),j)+S(modulo(i-1,N),j)+
     &S(i,modulo(j+1,N))+S(i,modulo(j-1,N)))
		expo=exp(-E/T)
		uno=1.d0	
		p=min(uno,expo)
c Generamos un número aleatorio entre 0 y 1. Si ese numero es menor a p entonces cambiamos el espin (i,j) de signo
		random=dran_u()

		if(random.lt.p) then
		S(i,j)=-S(i,j)
		endif

		enddo
		enddo
		
c Volcamos los datos a un fichero para realizar la gráfica con gnuplot
		contunos=0
		do i=0,N
			do j=0,N
				if(S(i,j).gt.0) then
				contunos=contunos+1
				write(7,*)i,j
			else
			write(8,*)i,j
			endif
		enddo
		enddo
		close(7)
		close(8)	
c plot 'unos.dat','nounos.dat'

		endprogram
