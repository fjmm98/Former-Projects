				
		program cuantica
		complex*16::a
		complex*16,parameter::i=(0,1)
		real*8::l,ko,s,pi,probabilidad,aN,C
		integer*8::N,nciclos,j,k,w,nd
		complex*16,ALLOCATABLE::alpha(:),beta(:),phi(:),xhi(:),gama(:)
c nohup planetas & y usamos top ,ps ux para ver los procesos, kill -9 PID: que es el numero que sale del proceso cuando ejecutamos TOP o ns etc.

c para que escriba a continuacion en los writes es open status='APPEND' y así no sobreescribe las líneas si abrimos y cerramos. CALLFLUSH(7) si la unidad es la 7 y haciendo eso
c lo vuelca en disco		
c		x=real(a)
c		y=imag(a)
c asi asociamos numeros a una variable compleja
c		a=complex(x,y)
c el modulo es
c		abs(a)
c el conjugado es 
c		conjg(a)
c ya tiene la sintaxis de los numeros complejos luego si a y b son complejos a/b es la división de números complejos.

c Hallamos el valor de pi que sera usado.
		pi=acos(-1.d0)
c Pedimos los valores que serán necesarios

		N=2000
!		write(6,*)'Introduzca el número de ciclos: '
!		read(5,*) nd=2.3N
		nd=4600
		nciclos=N/4
!1		write(6,*) 'Introduzca el valor de lambda: '
!	read(5,*) l
		l=0.1

c con estos valores realizamos los calculos que serán necesarios
		aN=dble(N)
		ko=(2.d0*pi*nciclos)/N
		s=(1.d0)/(4.d0*ko*ko)

c definimos la dimension de los vectores complejos
		allocate(alpha(0:(N-1)),beta(0:(N-1)))
		allocate(phi(0:N),xhi(0:N))
		allocate(gama(1:(N-1)))
c ponemos las condiciones iniciales
		phi(0)=0.d0
		phi(N)=0.d0
		xhi(0)=0.d0
		xhi(N)=0.d0	
		alpha(N-1)=0.d0
		beta(N-1)=0.d0
		probabilidad=0.d0
c calculamos las alphas y gamas que usaremos que son constantes durante todo el proceso
		do j=N-1,1,-1
			if (j.lt.ceiling(0.6d0*N).and.j.gt.floor(0.4d0*N)) then
			alpha(j-1)=(-1.d0)/(-2-l*ko*ko+alpha(j)+((2.d0*i)/s))
			else
			alpha(j-1)=(-1.d0)/(-2+alpha(j)+((2.d0*i)/s))
			endif
			gama(j)=(-1.d0)/alpha(j-1)
		enddo

c para comprobar los valores de alpha y gama		
c		do j=0,N-1
c		write(6,*)j, alpha(j)
c		enddo	
c		do j=1,N-1
c		write(6,*) j,gama(j)
c		enddo

c hallamos el valor de phi en cada jota inicial
		probabilidad=0.d0
		do j=1,N-1,1
			a=cos(ko*j)+i*sin(ko*j)
			phi(j)=a*exp((-8.d0*(4.d0*dble(j)-aN)**2.d0)/(aN*aN))
		enddo
		
c normalizamos la funcion de onda
		C=0.d0
		do j=1,N
		C=C+abs(phi(j))**2.d0
		enddo
			
	
c para comprobar el valor de phi
	
c		do j=1,N-1
c		write(6,*)j, phi(j)
c		enddo	
		open(UNIT=7,FILE="phi.txt",STATUS='OLD')

		do j=0,N
		write(7,*)abs(phi(j))
		enddo
		write(7,*)
		write(7,*)
		
		do k=1,nd
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


c hallamos la probabilidad
		probabilidad=0.d0
		do j=0,N,1
		probabilidad=probabilidad+(abs(phi(j)))**2
		enddo
		
c		write(6,*) 'La probabilidad es: ',k,probabilidad/C
		
c comprobamos el valor de phi en cada paso
c		do j=0,N
c		write(7,*) real(phi(j)),aimag(phi(j))
c		enddo
		do j=0,N
		write(7,*)(abs(phi(j))**2.d0)/(C)
		enddo
		write(7,*)
		write(7,*)
		 
c para hacer el gif antes poner en el gnuplot set terminal gif animate delay 5,set output "PHI1.gif" y luego el load 'cuanticacomand.dat'		
		enddo
		

		Close(7)
		


		end program
