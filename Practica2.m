clear all; clc;
%echo on 
%                          Intento 2
%
%                   ECUACION LINEAL DE CALOR

% Solucion de la ecuacion de calor mediante diferencias finitas.
% METODO IMPLICITO

%       Ut - Uxx=0        con a=1

%       U(x,0) = 2*x        si 0<=x<=1/2
%                2 - 2*x    si 1/2<x<=1
%       U(0,t)=  U(1,t)= 0    O
%       U(0,t)=sin(2wt);    U(1,t)=cos(wt/2)


%Definicion de parametros
t0  = 0.0;
tf  = 0.078;
%tf  = 0.15;
J   = 40;
dx  = 1.0/J;   
dt  = 0.00036;
N   = round(tf/dt)
mu  = dt/(dx)^2
w   = 40;

%  Matriz triangular
   D1=(1+2*mu)*ones(J-1,1);
   D2=(-mu)*ones(J-2,1);
   A=diag(D1)+diag(D2,1)+diag(D2,-1);
   
% Condiciones iniciales
x  = 0:dx:1;
U0(1:J/2+1)   = 2.0*x(1:J/2+1);
U0(J/2+2:J+1) = 2.0 - 2.0*x(J/2+2:J+1);
U0=U0';
figure(1)
plot(x,U0,'b-*')
   xlabel('x','FontSize',12);
   ylabel('U','FontSize',12,'VerticalAlignment','bottom');
    title('time = 0');


%disp('Presione cualquier tecla para continuar'); pause
clc ; 

%Algoritmo iterativo
%U0=U;
t=t0;
while(t<=tf)
  % Condiciones iniciales
  %U1
  % Condiciones de frontera  
  %U1(1,1)=0.0;
  %U1(J+1,1)=0.0;
  U1(1,1)=sin(2*w*t);
  U1(J+1,1)=cos(w/2 * t);
  % Construccion del termino independiente b
  b=U0(2:J,1);
  b(1)=b(1) + mu*U1(1,1);
  b(J-1)=b(J-1) + mu*U1(J+1,1);
  % Resolver el sistema AU1=b
  %U1(2:J,1)=A\b ;
  U1(2:J,1)=Thomas(D2,D1,D2,b);
  % Grafica de resultados al tiempo t
   figure(1)
   plot(x,U1,'b-*');
   axis([0 1 -1 1]);  
   axis(axis);
   xlabel('x','FontSize',12);
   ylabel('U','FontSize',12,'VerticalAlignment','bottom');
   title(sprintf('time = %f',t));
   pause(0.30);
  %Avanzar al siguiente paso de tiempo
   t=t+dt;
   U0=U1;
end


%Thomas Algorithm as function in matlab by Erick Rios
% Solves linear algebraic equation when the coefficient matrix is 
% tridiagonal. 
% Parameters :
% ------------
%   lowerDiagonal  : array
%                    contains the elements of the lower diagonal of matrix A.
%   mainDiagonal   : array 
%                    contains the elements of the main diagonal of the matrix
%                    A.
%   upperDiagonal  : array 
%                    contains the elements of the upper diagonal of the matrix
%                    A.
%   answerMatrix   : array 
%                    contains the elements of the answer matrix of the
%                    linear system. Can be represented by b.
% Returns    :
% ------------
%   solutionThomas : array 
%                    contains the solution of solve the linear equations
%                    using the Thomas Algorithm.
%                    
function solutionThomas = Thomas(lowerDiagonal,mainDiagonal,upperDiagonal,answerMatrix)

lengthMainDiagonal = length(mainDiagonal);
%Here we normalize the elements of the matrix creating
%two row vector
w = zeros(lengthMainDiagonal, 1) ; g = zeros(lengthMainDiagonal, 1);
w(1) = upperDiagonal(1)/mainDiagonal(1) ; g(1) = answerMatrix(1)/mainDiagonal(1);
%Creation of a column vector
if isrow(upperDiagonal)
    upperDiagonal = upperDiagonal';
end
if isrow(lowerDiagonal)
    lowerDiagonal = lowerDiagonal' ;
end

%Begins the loop
upperDiagonal = [upperDiagonal; 0]; 
lowerDiagonal = [0; lowerDiagonal] ;
for i=2:lengthMainDiagonal
    w(i) = upperDiagonal(i)/(mainDiagonal(i)-lowerDiagonal(i)*w(i-1)) ;
    g(i) = (answerMatrix(i)-lowerDiagonal(i)*g(i-1))/(mainDiagonal(i)-lowerDiagonal(i)*w(i-1)) ;
end

solutionThomas = zeros(lengthMainDiagonal, 1) ;
solutionThomas(lengthMainDiagonal) = g(lengthMainDiagonal) ;
for i=lengthMainDiagonal-1:-1:1
    solutionThomas(i) = -w(i)*solutionThomas(i+1)+g(i) ;
end
end

%echo off;
