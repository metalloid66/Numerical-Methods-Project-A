%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Task 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% to calculate condition number ...  norm(inv(A))*norm(A);

%                    start of matrix And Vector set up

n = 2; %starting number of equations 
lastN = 15; %last number of equations
vecOfNs = zeros(lastN,1); % a vector to keep count of what number of equations we're at, to use for plot later
vecOfvecOfResPartA = zeros(lastN,1); %a vector to keep the vector of residuum for each run of number of equations
vecOfvecOfResPartB = zeros(lastN,1);
vecOfvecOfResPartTest = zeros(lastN,1);
while n <= lastN
    %setting up part a) of task 2
    
    A = zeros(n,n);
    [Ai,Aj] = size(A);
    rowA = 1;
    while rowA <= Ai
        colA = 1;
        while colA <= Aj
            if rowA == colA
                A(rowA,colA) = 4;
            elseif rowA == colA -1
                A(rowA,colA) = -2;           %looping according to the given task
            elseif rowA == colA +1           %a[i,j] = 4 for i = j
                A(rowA,colA) = -2;           %a[i,j] = -2 for i=j-1 or i=j+1
            else                             %a[i,j] = 0 for all other cases
                A(rowA,colA) = 0 ;
            end
            colA = colA + 1;
        end
        rowA = rowA +1;
    end
    
    b = zeros(n,1);
    [Bi,Bj] = size(b);
    rowB = 1;
    while rowB <= Bi
        b(rowB,1) = -1 + (0.4 *rowB);
        rowB = rowB + 1;                     %looping according to the given task
        if rowB > n                          %b[i,1] = -1 + 0.4i
            break;
        end
    end
    
    Atest = [0 1 2 3; 0 1 4 12; 1 1 1 1; 1 2 4 8]; % a functions texting matrix (mainly for gewpp)
    btest = [-0.2 0.8 1.5 1.2]';
    
    
    % setting up part b) of task 2
    A2 = zeros(n,n);
    [Ai,Aj] = size(A2);
    rowA = 1;
    while rowA <=Ai
        colA = 1;                          %looping according to the given task
        while colA <= Aj                   %a[i,j] = 2/[5(i+j+1)]
            A2(rowA,colA) = 2 /(5 *(rowA + colA + 1));
            colA = colA + 1;
        end
        rowA = rowA + 1;
    end
    
    b2 = zeros(n,1);
    [Bi,Bj] = size(b2);
    rowB = 1;
    while rowB <= Bi
        k = 0;
        evenOrOdd = 1;
        while k < rowB
            evenOrOdd = evenOrOdd * -1;
            k = k+1;                           %looping according to the given task
        end                                 %b[i,1] = 8/(7*i) for i even
        if evenOrOdd == -1 %odd             %b[i,1] = 0 for i odd
            b2(rowB,1) = 0;
        else               %even
            b2(rowB,1) = 8/(7*rowB);
        end
        rowB = rowB + 1;
    end
    %                    End of Matrix and Vector set up
    
    %                     The main program test
    
    
    resultsPartA = gewpp(A,b);
    resultsPartB = gewpp(A2,b2);
    resultsTest = gewpp(Atest,btest);
    
    vecOfResPartA = errFun(resultsPartA,A,b);
    vecOfResPartB = errFun(resultsPartB,A2,b2);
    vecOfResPartTest = errFun(resultsTest,Atest,btest);
    
    vecOfNs(n,1) = n;
    vecOfvecOfResPartA(n,1) = norm(vecOfResPartA);
    vecOfvecOfResPartB(n,1) = norm(vecOfResPartB);
    vecOfvecOfResPartTest(n,1) = norm(vecOfResPartTest);
    
    n = n+1;
    if n == 11
        resCorrectionFun(resultsPartA,A,b,vecOfResPartA,15);
        resCorrectionFun(resultsPartB,A2,b2,vecOfResPartB,15);
%       resCorrectionFun(resultsTest,Atest,btest,vecOfResPartTest);
    end
    
end
%plotting

 
figure(1)
plot(vecOfNs,vecOfvecOfResPartA);
xlabel('Number of equations n');
ylabel('Euclidean norm of the vector of residuum for n');
title('Task2 a)');

figure(2)
plot(vecOfNs,vecOfvecOfResPartB);
xlabel('Number of equations n');
ylabel('Euclidean norm of the vector of residuum for n');
title('Task2 b)');




%                    End Of The main program test




%        Start of Gaussian Elimination with Partial Pivoting Function

%start of forward elimiation step

function results = gewpp(Amat,Bvec) %Matrix A, vector b, number of iterations n

D = [Amat,Bvec];
[i,j] = size(D);
col = 1;
while col <= j
    k = D(col,col); %K11, K22, K33 ....
    row = col;
    while 1
        v = D(row,col);
        if abs(k) < abs(v)
            %          t = ['value bigger than k',num2str(col),' is: ',num2str(v), ' procceding to swap the values row with Ks value row'];
            %           disp(t);
            D([col row],:) = D([row col],:); %swapping the rows
            %        disp(['D is now swapped: ']);
            %        disp(D);
            k = D(col,col);
            %        disp(['k is now: ', num2str(k)]);
            for m=col+1:i
                D(m,:) = D(m,:) - ((D(m,col)/k) * D(col,:));
            end
            break;
        end
        if row >= i
            %             disp('Could not find a value bigger than k.. procceding to formula');
            for m=col+1:i
                D(m,:) = D(m,:) - ((D(m,col)/k) * D(col,:));
            end
            break;
        end
        row = row+1;
    end
    col = col+1;
    %disp(D);
    if col >= j
        %         disp('last column has been reached');
        break;
    end
end

%end of forward elimination step


%start of backward substitution step

U = D(:,1:j-1);
C = D(:,j:j);
[Ci,Cj] = size(C);
[Ui,Uj] = size(U);

rowU = Ui;
colU = Uj;
rowC = Ci;
X = ((C(rowC,1) - (U(rowU,Ui) * 0)) / U(rowU,colU));
results = zeros(Ci,1);
results(rowC,1) = X;
while rowC >=1
    rowC = rowC - 1;
    rowU = rowU - 1;
    colU = colU - 1;
    Ujjer = Uj;
    p = 0;
    while 1
        p = p - U(rowU,Ujjer)*results(Ujjer,1);
        Ujjer = Ujjer - 1;
        if Ujjer == rowU
            break;
        end
    end
    Xi = ((C(rowC,1) - -p ) / U(rowU,colU));
    results(rowC,1) = Xi;
    if rowU == 1
        break;
    end
end
end
% disp(['The results for n = ',num2str(n),' are:']);
% disp(results);
%end of backward substitution step...results ready

%        End of Gaussian Elimination with Partial Pivoting Function




%                    Start of Solution Error Function
function vecOfRes = errFun(results,matA,vecB)
[i,j] = size(results);
vecOfRes = zeros(i,1);
rowA = 1;
rowB = 1;
counter = 0;
while rowB <= i
    counter = counter + 1;
    colA = 1;
    rowR = 1;
    temp = 0;
    
    while rowR <= i
        temp = temp + (matA(rowA,colA)*results(rowR,1));
        colA = colA + 1;
        rowR = rowR + 1;
    end
    
    vecOfRes(rowB,1) = abs(temp - vecB(counter,1));
    rowB = rowB + 1;
    rowA = rowA + 1;
end
end
%                    End of Solution Error Function

%                    Start of the residual correction function

function resCorrectionFun(results,matA,vecB,vecOfRes,itr)
tempvecOfRes = vecOfRes;
disp('The vector of residuum r^(1) for n = 10 is: ');
disp(abs(vecOfRes));
disp('The norm of the vector of residuum r^(1) for n = 10 is: ');
disp(norm((vecOfRes)));
counter = 1;
while counter <= itr
deltaX = gewpp(matA,tempvecOfRes);
results = results - deltaX;
tempvecOfRes = errFun(results,matA,vecB);
disp('Residual Correction Performed');
disp(['The norm of the vector of residuum after ',num2str(counter),' residual correction(s) is: ']);
disp(norm(tempvecOfRes));
counter = counter + 1;
end
end
%                    End of the residual correction function