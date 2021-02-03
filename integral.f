		program integral
		real*8::a,b,dx,a1,f,x1,p,p1,PI,sum,xi,int
		integer*8::i,N
	    write(6,*) 'Introduzca el límite inferior: '
		read(5,*) a
		write(6,*) 'Introduzca el límite superior: '
		read(5,*) b
		write(6,*) 'Introduzca el nº de pasos: '
		read(5,*) N
		dx=(b-a)/dble(N)

c ponemos dble(N) para que lo tome como real de doble precision y ejecute correctamente la división
c para ver el tiempo que tarda ponemos time integral 
	
		sum=0.d0
		do i=1,N-1
c		if(a<0) exit para parar el bucle, pero esto hace que lo compruebe en cada paso, en caso de tener dos ciclos y querer dejar fijo un 
c       valor de uno de ellos ponemos if(a<0) cycle 

		xi=a+dx*i
		sum=sum+f(xi)
		enddo
		int=dx*(f(b)-f(a))/2.d0+dx*sum
		write(6,*) 'La integral es: '
		write(6,*) int

c		x1=3.d0
c		a1=f(x1)
c		PI=4.d0*arcsin(1.d0)

c		write(6,*)a1

c		call G(x1,p1)
c		write(6,*)p1

		end program
		
		real*8 function f(x)
		real*8::x
		f=x**2.d0
		end function
		

		subroutine G(x,p)
		real*8:: x,p
		p=x**2.d0
		return	
		end subroutine G
