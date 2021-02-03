        program planetas
        real*8, allocatable::x(:),y(:),vx(:),vy(:),ax(:),ay(:),m(:)
        real*8, allocatable::d(:),wx(:),wy(:)
        real*8::h,xx,yy,vv,bb,mm,den,tiempo,t,ti,vi
        integer*8::i,a,j,k
c       lleno los vectores con los datos de posicion velocidad y masa
        write(6,*)'hola, que le damos a la h'
        read(5,*)h
        write(6,*)'hola, que le damos al tiempo'
        read(5,*)tiempo
        OPEN(UNIT=7,FILE='datos.txt',STATUS='OLD')
        READ(7,*)a
        READ(7,*)
        READ(7,*)
        allocate(x(a),y(a),vx(a),vy(a),ax(a),ay(a),m(a),d(a))
        allocate(wx(a),wy(a))
        do i=1,a
        READ(7,*)xx,yy,vv,bb,mm
        x(i)=xx
        y(i)=yy
        vx(i)=vv
        vy(i)=bb
        m(i)=mm  
        enddo
        CLOSE(7)
c       Reescalamientos del tiempo y velocidad
        ti= 0.01719792
        h=h*ti
        tiempo=tiempo*ti
        vi=0.00003358037
        vy=vi*vy
c       1.calculamos las aceleraciones
c       calculamos las distancias
        ax=0.
        ay=0.
        do i=1,a
        do j=1,a
            if(i.eq.j)cycle
        den=sqrt(((x(i)-x(j))**2.0+(y(i)-y(j))**2.0)**3.0)             
        ax(i)=ax(i)-m(j)*(x(i)-x(j))*(1.0/den)
        ay(i)=ay(i)-m(j)*(y(i)-y(j))*(1.0/den)
   
        enddo
        enddo  
c       calculamos las vueltas de iteracciones
        t=0.
        do while(t.le.tiempo)
c       2.Calculamos las nuevas posiciones y las w


        x=x+h*vx+((h*h)/2.0)*ax
        y=y+h*vy+((h*h)/2.0)*ay
        wx=vx+(h/2.0)*ax
        wy=vy+(h/2.0)*ay
        
c       3.evaluamos de nuevo la aceleracion con nuevas posiciones
        ax=0.
        ay=0.
        do i=1,a
        do j=1,a
            if(i.eq.j)cycle
        den=sqrt(((x(i)-x(j))**2.0+(y(i)-y(j))**2.0)**3.0)             
        ax(i)=ax(i)-m(j)*(x(i)-x(j))*(1.0/den)
        ay(i)=ay(i)-m(j)*(y(i)-y(j))*(1.0/den)
 
        enddo
        enddo  

c       4.evaluamos las nuevas velocidades
        vx=0.
        vy=0.

        vx=wx+(h/2)*ax
        vy=wy+(h/2)*ay
 
        t=t+h
        write(6,*)x(4),y(4)
        enddo
        read(5,*)mm

        end
    

