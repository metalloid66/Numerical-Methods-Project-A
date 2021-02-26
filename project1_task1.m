%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Task 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Task#1: Finding macheps');
macheps = 1;
while 1.000000000 + (macheps/2) > 1.000000000
    macheps = macheps ./ 2;
    disp('progress of showing macheps:')
    disp(macheps)
end
format long;
disp('Last value before the termination of the loop:');
disp(macheps);
disp('Verfication using matlab eps function: ')
matlabeps = eps;
disp(matlabeps)
if matlabeps == macheps
    disp('Our macheps is verified')
end