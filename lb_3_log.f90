program lb_3
	use omp_lib
	implicit none
	real:: t1, t2, h, tau, maxT
	integer:: steps_n, steps_tau, i, j, k

	real, allocatable:: u(:,:,:), u_reference(:,:,:), x(:), y(:), t(:), f(:,:,:)


	write(*,*) "Дано:"
	write(*,*) "U(x,y,t)=lg((t+1)/2)*cos(x)*cos(y)"
	write(*,*)
	
	write(*,*) "Находим f:"
	write(*,*) "f(x,y,t)=cos(x)*cos(y)*(1/(ln(10.0)*(t+1))+2*lg((t+1)/2))"
	write(*,*)
	
	write(*,*) "Начальные условия:"
	write(*,*) "U(x,y,0)=lg(0.5)*cos(x)*cos(y)"
	write(*,*)
	
	write(*,*) "Граничные условия: "
	write(*,*) "U(0,y,t)=log10((t+1)/2)*cos(y)"
	write(*,*) "U(1,y,t)=log10((t+1)/2)*cos(1.0)*cos(y)"
	write(*,*) "U(x,0,t)=log10((t+1)/2)*cos(x)"
	write(*,*) "U(x,1,t)=log10((t+1)/2)*cos(x)*cos(1.0)"
	write(*,*)
	
	write(*,*) "Введите h: "
	read(*,*) h
	write(*,*)

	! Количество шагов по осям x и y
	steps_n = INT(1.0/h)+1;
	
	write(*,*) "Введите ограничение по времени maxT: "
	read(*,*) maxT
	write(*,*)

	write(*,*) "Введите tau: "
	read(*,*) tau
	write(*,*)

	! Количество шагов по оси t
	steps_tau = INT(maxT/tau)+1;
	
	allocate(u(steps_n,steps_n,steps_tau), u_reference(steps_n,steps_n,steps_tau), & 
				x(steps_n), y(steps_n), t(steps_tau), f(steps_n,steps_n,steps_tau))
	
	! Заполняем точки по осям x и y
	forall(i=1:steps_n)
		x(i)=(i-1)*h;
		y(i)=(i-1)*h;
	end forall

	! Заполняем точки по оси времени
	forall(k=1:steps_tau)
		t(k)=(k-1)*tau;
	end forall

	! Заполняем значения f
	forall(i=1:steps_n, k=1:steps_tau)
		f(i,j,k)=cos(x(i))*cos(y(i))*(1/(log(10.0)*(t(k)+1))+2*log10((t(k)+1)/2));
	end forall

	! Заполняем значения эталонной функции
	forall(i=1:steps_n, k=1:steps_tau)
		u_reference(i,j,k)=log10((t(k)+1)/2)*cos(x(i))*cos(y(i));
	end forall

	! Заполняем начальное граничное условие по времени
	forall(i=1:steps_n)
		u(i,j,1)=cos(x(i))*cos(y(i));
	end forall

	! Заполняем граничные условия
	forall(i=1:steps_n, k=2:steps_tau)
		u(1,i,k)=log10((t(k)+1)/2)*cos(y(i));
		u(steps_n,i,k)=log10((t(k)+1)/2)*cos(1.0)*cos(y(i));
		u(i,1,k)=log10((t(k)+1)/2)*cos(x(i));
		u(i,steps_n,k)=log10((t(k)+1)/2)*cos(x(i))*cos(1.0);
	end forall

	call cpu_time(t1)
	do k=2,steps_tau
		do i=2,steps_n-1
			do j=2,steps_n-1
				! Заполняем оставшиеся значения
				u(i,j,k)=((u(i+1,j,k-1)-4*u(i,j,k-1)+u(i-1,j,k-1)+u(i,j+1,k-1)+u(i,j-1,k-1))/(h**2)+f(i,j,k-1))*tau+u(i,j,k-1);
			end do
		end do
	end do
	call cpu_time(t2)
		
	OPEN(10,FILE='lb_3_plot_data.dat')

	WRITE(10,*) h
	WRITE(10,*) steps_n
	WRITE(10,*) tau
	WRITE(10,*) steps_tau
	
	do k=1,steps_tau
		do i=1,steps_n
			do j=1,steps_n
				WRITE(10,*) u(i, j, k)
			end do
		end do
	end do

	CLOSE(10)

	write(*,*) "Последовательное вычисление:"

	write(*,*) "Затраченное время: ", t2-t1, " секунд"
	write(*,*) "Максимальное абсолютное отклонение от эталона", &
				MAXVAL(ABS(u_reference(:,:,:)-u(:,:,:)))
	write(*,*)
	
	u(:,:,:)=0.0;

	! Заполняем начальное граничное условие по времени
	forall(i=1:steps_n)
		u(i,j,1)=cos(x(i))*cos(y(i));
	end forall

	! Заполняем граничные условия
	forall(i=1:steps_n, k=2:steps_tau)
		u(1,i,k)=log10((t(k)+1)/2)*cos(y(i));
		u(steps_n,i,k)=log10((t(k)+1)/2)*cos(1.0)*cos(y(i));
		u(i,1,k)=log10((t(k)+1)/2)*cos(x(i));
		u(i,steps_n,k)=log10((t(k)+1)/2)*cos(x(i))*cos(1.0);
	end forall
	
	call cpu_time(t1)
	do k=2,steps_tau
		!$OMP PARALLEL SHARED(u,f,tau,h)
		!$OMP WORKSHARE
		forall (i=2:steps_n-1,j=2:steps_n-1)
				u(i,j,k)=((u(i+1,j,k-1)-4*u(i,j,k-1)+u(i-1,j,k-1)+u(i,j+1,k-1)+u(i,j-1,k-1))/(h**2)+f(i,j,k-1))*tau+u(i,j,k-1);		
		end forall
		!$OMP END WORKSHARE
		!$OMP END PARALLEL
	end do
	call cpu_time(t2)

	write(*,*) "Параллельное вычисление:"

	write(*,*) "Затраченное время: ", t2-t1, " секунд"
	write(*,*) "Максимальное абсолютное отклонение от эталона", &
				MAXVAL(ABS(u_reference(:,:,:)-u(:,:,:)))

	deallocate(u, u_reference, x, y, t, f)

end program lb_3
