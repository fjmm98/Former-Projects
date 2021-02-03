                                         		program cohete
		real*8::nu,rt,rl,Dtl,Lambda,G,Mt,Ml,w,phi,theta,pr,pfi,e,r,energia
		real*8::e1,e2,e3,e4,f01,f02,f03,f04
		real*8::h,hmax,t,tmax,emax,s,rprima,f1,f2,f3,f4,vx,vy,a,b
		real*8::f01bis2,f02bis2,f03bis2,f04bis2,f01bis1,f02bis1,f03bis1,f04bis1

! con el comando -ffree-form los comentarios van con ! y no con c
! phi y theta son las condiciones iniciales

! 	contains y asi las subrutinas usan las variables el contains dentro del cuerpo del programa no después del end program
!		write(6,*) 'Introduce el paso a usar (En minutos): '
!		read(5,*) h
!		write(6,*)'Introduce el tiempo maximo (En minutos): '
!		read(5,*) tmax
!		write(6,*) 'Introduce phi (En radianes*PI): '
!		read(5,*) phi
!		write(6,*) 'Introduce theta (En radianes*PI): '
!		read(5,*) theta

! definimos las constantes
		h=1
		tmax=500
		G=6.67e-11
		Mt=5.9736e24
		Ml=0.07349e24
		w=2.6617e-6
		Dtl=3.844e8
		rt=6.37816e6
		rl=1.7374e6
		Lambda=7.014744e-12
		nu=0.0123024
!		emax=h**5.d0
		emax=1
		r=rt
		pr=5
		n=500

! Introducimos los valores iniciales
! f1=r f2=pr f3=phi f4=pfi

		f03=0.5399
		f04=3.2e-15
		f02=2.9e-5
		f01=0.016592507804370448

! energia inicial
		energia=0
		energia=energia+dtl*dtl*0.5*(f02*f02+f04*f04)
		energia=energia-6.67e-11/dtl*(7.349e22+5.9736e24/f01)
		energia=energia-w*dtl*dtl*f04
		WRITE(6,*)'Energía inicial: ',energia		


! Realizamos el programa para h fija sin adaptar

		t=0.d0
		open(UNIT=8,FILE="hfija.dat",STATUS='NEW')

		do while(t.lt.tmax)
		CALL RK(t,h,f01,f02,f03,f04)
		write(8,*)f01,f02,f03,f04
		t=t+h
		enddo
		close(8)

		energia=0
		energia=energia+dtl*dtl*0.5*(f02*f02+f04*f04)
		energia=energia-6.67e-11/dtl*(7.349e22+5.9736e24/f01)
		energia=energia-w*dtl*dtl*f04
		WRITE(6,*)'Energía final (h fija): ',energia		



! limpiamos las variables
! f1=r f2=pr f3=phi f4=pfi

		f03=0.5399
		f04=3.2e-15
		f02=2.9e-5
		f01=0.01659
! Realizamos el programa para h adaptada
		t=0.d0

		open(UNIT=9,FILE='hadaptada.dat',STATUS='NEW')
!		do while(t.lt.tmax)
		do i=1,500
! f1=r f2=pr f3=phi f4=pfi
! cambiar los f0 por sus valores iniciales
!		f01bis1=f01
!		f02bis1=f02
!		f03bis1=f03
!		f04bis1=f04

		f01bis2=f01
		f02bis2=f02
		f03bis2=f03
		f04bis2=f04

		CALL RK(t,h,f01,f02,f03,f04)

		CALL RK(t,(h/2.d0),f01bis2,f02bis2,f03bis2,f04bis2)

		CALL RK(t,(h/2.d0),f01bis2,f02bis2,f03bis2,f04bis2)

! guardamos los datos del primero y luego hacemos los dos pasos con h/2

		e1=(16.d0/15.d0)*abs(f01-f01bis2)
		e2=(16.d0/15.d0)*abs(f02-f02bis2)
		e3=(16.d0/15.d0)*abs(f03-f03bis2)
		e4=(16.d0/15.d0)*abs(f04-f04bis2)

! usamos el error maximo de estos

		e=max(e1,e2,e3,e4)
		b=(e/emax)**0.2
		s=max(a,0.00000001)
		hmax=h/s
		write(6,*)h
		if(s.gt.2) then
		h=hmax
		cycle
		endif

!		if(s.lt.2) then
		t=t+h
		f01=f01bis2
		f02=f02bis2
		f03=f03bis2
		f04=f04bis2
!		endif


		if(h.lt.hmax) then
		h=2.d0*h
		cycle
		endif

		write(9,*)f01,f02,f03,f04

		t=t+h

		enddo
	
		energia=0
		energia=energia+dtl*dtl*0.5*(f02*f02+f04*f04)
		energia=energia-6.67e-11/dtl*(7.349e22+5.9736e24/f01)
		energia=energia-w*dtl*dtl*f04
		WRITE(6,*)'Energía final (h adaptada): ',energia	

		close(9)

		end program



! f1=r f2=pr f3=phi f4=pfi
! definimos las funciones como las ecuaciones diferenciales, f1 es la ecuacion diferencial de r y asi sucesivamente
		real*8 function f1(r,pr,phi,pfi,t)
		real*8::r,pr,phi,pfi,t
		f1=pr
		end function f1

		real*8 function f2(r,pr,phi,pfi,t)
		real*8::pfi,r,rprima,phi,t,Lambda,nu,w,a,pr
		Lambda=7.014744e-12
		nu=0.0123024
		w=2.6617e-6
		a=r-cos(phi-w*t)
		rprima=sqrt(1+(r**2.d0)-2.d0*r*cos(phi-w*t))
		bis=Lambda*(((1/(r**2.d0))+((nu*a)/(rprima**3.d0))))
		f2=((pfi**2.d0)/(r**3.d0))-bis
		end function f2

		real*8 function f3(r,pr,phi,pfi,t)
		real*8::r,pr,phi,pfi,t
		f3=(pfi/(r**2.d0))
		endfunction f3

! puede haber fallos en las funciones por usar cosas que no estan bien definidas

		real*8 function f4(r,pr,phi,pfi,t)
		real*8::rprima,t,phi,nu,Lambda,w,r,pr,pfi
		Lambda=7.014744e-12
		nu=0.0123024
		w=2.6617e-6
		rprima=sqrt(1+(r**2.d0)-2.d0*r*cos(phi-w*t))
		f4=((-Lambda*nu*r*sin(phi-w*t))/(rprima**3.d0))
		end function f4

! definimos la subrutina para hacer el runge kutta de 4º orden

		subroutine RK(t,h,f01,f02,f03,f04)
		real*8::t,h,f01,f02,f03,f04,f1,f2,f3,f4,r,pr,phi,pfi
		real*8,DIMENSION(4)::k1,k2,k3,k4
! introducimos las condiciones iniciales
! calculamos k1
		k1(1)=h*f01
		k1(2)=h*f02
		k1(3)=h*f03
		k1(4)=h*f04
! calculamos k2
		k2(1)=h*f1(f01+(k1(1)/2.d0),f02+(k1(2)/2.d0),f03+
     &(k1(3)/2.d0),f04+(k1(4)/2.d0),t+(h/2.d0))
		k2(2)=h*f2(f01+(k1(1)/2.d0),f02+(k1(2)/2.d0),f03+
     &(k1(3)/2.d0),f04+(k1(4)/2.d0),t+(h/2.d0))
		k2(3)=h*f3(f01+(k1(1)/2.d0),f02+(k1(2)/2.d0),f03+
     &(k1(3)/2.d0),f04+(k1(4)/2.d0),t+(h/2.d0))
		k2(4)=h*f4(f01+(k1(1)/2.d0),f02+(k1(2)/2.d0),f03+
     &(k1(3)/2.d0),f04+(k1(4)/2.d0),t+(h/2.d0))
! calculamos k3
		k3(1)=h*f1(f01+(k2(1)/2.d0),f02+(k2(2)/2.d0),
     &f03+(k2(3)/2.d0),f04+(k2(4)/2.d0),t+(h/2.d0))
		k3(2)=h*f2(f01+(k2(1)/2.d0),f02+(k2(2)/2.d0),
     &f03+(k2(3)/2.d0),f04+(k2(4)/2.d0),t+(h/2.d0))
		k3(3)=h*f3(f01+(k2(1)/2.d0),f02+(k2(2)/2.d0),
     &f03+(k2(3)/2.d0),f04+(k2(4)/2.d0),t+(h/2.d0))
		k3(4)=h*f4(f01+(k2(1)/2.d0),f02+(k2(2)/2.d0),
     &f03+(k2(3)/2.d0),f04+(k2(4)/2.d0),t+(h/2.d0))
! calculamos k4
		k4(1)=h*f1(f01+k3(1),f02+k3(2),f03+k3(3),f04+k3(4),t+h)
		k4(2)=h*f2(f01+k3(1),f02+k3(2),f03+k3(3),f04+k3(4),t+h)
		k4(3)=h*f3(f01+k3(1),f02+k3(2),f03+k3(3),f04+k3(4),t+h)
		k4(4)=h*f4(f01+k3(1),f02+k3(2),f03+k3(3),f04+k3(4),t+h)
! Hallamos la f evolucionada
		f01=f01+(1.d0/6.d0)*(k1(1)+2.d0*(k2(2)+k3(1))+k4(1))
		f02=f02+(1.d0/6.d0)*(k1(2)+2.d0*(k2(2)+k3(2))+k4(2))
		f03=f03+(1.d0/6.d0)*(k1(3)+2.d0*(k2(3)+k3(3))+k4(3))
		f04=f04+(1.d0/6.d0)*(k1(4)+2.d0*(k2(4)+k3(4))+k4(4))
! terminado el runge-kutta

		end subroutine RK

