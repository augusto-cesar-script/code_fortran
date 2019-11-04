PROGRAM GFDM_diff_explicito
REAL, DIMENSION(2,2)::A 
INTEGER N, i
REAL:: h(4), x(4), w(4)

N = 4
x(1) = 1
x(2) = -1.2
x(3) = 2.5
x(4) = -2
xo = 0

Do i= 1,N
    h(i) = x(i) - xo
    w(i) = 1/h(i)
ENDDO
  
A(1,1) = 0
A(2,1) = 0
A(2,2) = 0
DO i = 1, N
    A(1,1) = A(1,1) + h(i)**2*w(i)**2
    A(2,1) = A(2,1) + h(i)**3*w(i)**2
    A(2,2) = A(2,2) + h(i)**4*w(i)**2/4
ENDDO

A(1,2)=A(2,1)

PRINT*, 'Matriz A: ', A

l(1,1) = SQRT(A(1,1))
l(2,1) = A(1,2)/SQRT(A(1,1))
l(2,2) = SQRT(A(2,2) - l(1,2)**2)
L2(1,1) = l(1,1)
L2(2,1) = l(2,1)
L2(1,2) = 0
L2(2,2) = l(2,2)
L2T = TRANSPOSE(L2)

FUNCTION forward_substitution(L2,b)
REAL, DIMENSION(2,1)::Y 
REAL sumat
INTEGER i, j, N
Y(1) = b(1)/L2(1,1) ! defining the [Y] vector
! solving for the remaining N-1 unknow
DO i= 1, N
    sumat = 0
    DO j= 1, i-1
        sumat = sumat + L2(i,j)*Y(j)
    ENDDO
    Y(i) = (b(i) - sumat)/L(i,i)
    forward_substitution = Y
ENDDO
return
END


END
    
