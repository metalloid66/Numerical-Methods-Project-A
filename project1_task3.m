%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Task 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                  Start of Matrix And Vector set up

A = [7 2 1 -2; 4 -10 3 -1; 2 -1 -8 -1; 5 -2 1 -8];
b = [6 8 0 2]';

%                  End of Matrix And Vector set up




%                       The main program test

[L,U,D] = decomposeLUD(A);
jacobi(A,b,L,U,D);
gaussSeidel(A,b,L,U,D);


%                    End of The main program test




%                    Start of Matrix Decomposition Function

function [L,U,D] = decomposeLUD(vecA)
%getting lower triangular matrix L
[i,j] = size(vecA);
rowL = 1;
colL = 1;
L = vecA;
    if rowL == 1 
        while colL <= j 
         L(rowL,colL) = 0;
         colL = colL+1;
        end
    end
    rowL = rowL + 1;
    colL = rowL;
    if rowL >= 2 
        while rowL<=i
          while colL<=j
        L(rowL,colL)=0;
        colL = colL + 1;
          end
         rowL = rowL + 1;
         colL = rowL;
        end
    end
%getting upper triangular matrix U
colU = i;
rowU = j;
U = vecA;
    if rowU == j
        while colU >= 1 
         U(rowU,colU) = 0;
         colU = colU-1;
        end
    end
    rowU = rowU - 1;
    colU = rowU;
    if rowU <= 3 
        while rowU>=1
          while colU>=1
        U(rowU,colU)= 0;
        colU = colU - 1;
          end
         rowU = rowU - 1;
         colU = rowU;
        end
    end
%getting Diagonal matrix D 
colD = 1;
rowD = 1;
D = vecA;
        while rowD<=i
          while colD<=j 
              if rowD == colD
              D(rowD,colD) = D(rowD,colD);
              else 
              D(rowD,colD) = 0;
              end
              colD = colD + 1;
          end
        rowD = rowD + 1;
        colD = 1;
        end
end

%                    End of Matrix Decomposition Function




%                    Start of Jacobi Iteration Function

function jacobi(matA,vecB,L,U,D)
[i,j] = size(matA);
Atemp = matA;
Xr = zeros(i)'; %starting guess vector
k = 1; %number of iterations
kVec = zeros(25,1); %vector for number of iterations
errNormVec = zeros(25,1); %vector for norms of errors
while 1
Xr = (inv(-D)*((L+U)*Xr) + (inv(D)*vecB));
rowA = 1;
colA = 1;
while rowA <=i
    while colA<=j
        Atemp(rowA,colA) = matA(rowA,colA)*Xr(colA,1);
        colA = colA + 1;
    end
    colA = 1;
    rowA = rowA + 1;
end
AtempSum = sum(Atemp,2);
err = abs(vecB - AtempSum);
errNorm = norm(err);

kVec(k,1) = k;
errNormVec(k,1) = errNorm;

if errNorm <= 10^-10
    break;
end
k = k+1;
end
figure(1)
plot(kVec,errNormVec);
xlabel('Number of iterations K');
ylabel('Norm of the solution err');
title('Jacobi Iteration');
end

%                    End of Jacobi Iteration Function




%                    Start of Gauss-Seidel Iteration Function

function gaussSeidel(matA,vecB,L,U,D)
[i,j] = size(matA);
Atemp = matA;
Xr = zeros(i)'; %starting guess vector
k = 1; %number of iterations
kVec = zeros(25,1); %vector for number of iterations
errNormVec = zeros(25,1); %vector for norms of errors
while 1
Xr = -(inv(D+L)*(U*Xr)) + (inv(D+L) * vecB);
rowA = 1;
colA = 1;
while rowA <=i
    while colA<=j
        Atemp(rowA,colA) = matA(rowA,colA)*Xr(colA,1);
        colA = colA + 1;
    end
    colA = 1;
    rowA = rowA + 1;
end
AtempSum = sum(Atemp,2);
err = abs(vecB - AtempSum);
errNorm = norm(err);

kVec(k,1) = k;
errNormVec(k,1) = errNorm;
if errNorm <= 10^-10
    break;
end
k = k+1;
end
figure(2);
plot(kVec,errNormVec);
xlabel('Number of iterations K');
ylabel('Norm of the solution err');
title('Gauss-Seidel Iteration');
end

%                    End of Gauss-Seidel Iteration Function   










