%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Task 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                  Matrix set up
A = [14 26 22 16 22;
    26 50 46 28 40; 
    22 46 50 20 32;
    16 28 20 20 26; 
    22 40 32 26 35];


%                    The main program test
AeigValuesShift = EigvalQRshifts(A,10^-6,1000);
AeigValuesNoShift = EigvalQRNoShift(A,10^-6,1000);

%                    Start of the QR factorization Function

function [Q,R]=QRfactor(A)
[m,n]=size(A);
R=A; %Start with R=A
Q=eye(m); %Set Q as the identity matrix
for k=1:m-1
    x=zeros(m,1);
    x(k:m,1)=R(k:m,k);
    g=norm(x);
    v=x; v(k)=x(k)+g;
    %Orthogonal transformation matrix that eliminates one element
    %below the diagonal of the matrix it is post-multiplying:
    s=norm(v);
    if s~=0, w=v/s; u=2*R'*w;
        R=R-w*u'; %Product HR
        Q=Q-2*Q*(w*w'); %Product QR
    end
end
end

%                    End of the QR factorization Function


%                    Start of roots of a quadratic polynomial function

function [x1,x2] = quadpolynroots(A,B,C)
    l1 = -B + sqrt(B*B - 4*A*C);
    l2 = -B - sqrt(B*B - 4*A*C);
    %we choose a counter with a larger module
    if abs(l1) > abs(l2)
       counter = l1;
    else
        counter = l2;
    end
    x1 = counter/(2*A);
    % second root is calculated from Viete's formula
    x2 = ((-B)/A) - x1;
end

%                    End of roots of a quadratic polynomial function


%                    Start of finding Eigenvalue with no shift function

function eigenvalues=EigvalQRNoShift(D,tol,imax)
%tol - tolerance (upper bound) for nulled elements
%imax - maximal number of iterations
n=size(D,1);
i=1;
itrNum = 0;
while i <= imax && max(max(D-diag(diag(D)))) > tol
[Q1,R1]=QRfactor(D);
D=R1*Q1; %transformed matrix
i=i+1;
itrNum = i;
end
fprintf('iterations needed to force all off-diagonal \n elements below abs of tolerance value: \n');
disp(itrNum);
disp('final matrix');
disp(D);
if i > imax
error('imax exceeded program terminated');
end
eigenvalues=diag(D);
end

%                    End of finding Eigenvalue with no shift function


%                    Start of finding Eigenvalue with shift function

function eigenvalues=EigvalQRshifts(A,tol,imax)
n=size(A,1);
eigenvalues=diag(ones(n));
INITIALsubmatrix=A;
itrNum = 0;
for k=n:-1:2
    DK=INITIALsubmatrix;
    i=0;
    disp(DK);
    while i<=imax && max(abs(DK(k,1:k-1)))>tol
        DD=DK(k-1:k,k-1:k);
        [ev1,ev2]=quadpolynroots(1,-(DD(1,1)+DD(2,2)),DD(2,2)*DD(1,1)-DD(2,1)*DD(1,2));
        if abs(ev1-DD(2,2)) < abs(ev2-DD(2,2))
            shift=ev1;
        else 
            shift=ev2;
        end
        DP=DK-eye(k)*shift;
        [Q1,R1]=QRfactor(DP);
        DK=R1*Q1+eye(k)*shift;
        i=i+1;        
    end
    itrNum = itrNum + i;
%     disp('iterations needed to force all off-diagonal elements below abs of tolerance value: ');
    
    if i > imax
        error('imax exceeded program terminated');
    end
    
    eigenvalues(k)=DK(k,k);
    if k > 2
     INITIALsubmatrix=DK(1:k-1,1:k-1);
    else
        eigenvalues(1)=DK(1,1);
    end
end
fprintf('iterations needed to force all off-diagonal \n elements below abs of tolerance value: \n');
disp(itrNum);
disp('final matrix');
disp(DK);
end
%                    End of finding Eigenvalue with shift function
