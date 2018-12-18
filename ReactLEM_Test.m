%Testing ReactLEM.m creates local reation element vector correctly

%Set up initial parameters needed for the test
%%TEST FOR LINEAR AND QUADRATIC CASES
%Create Quadratic mesh
order = 2;  %Quadratic order
xmin = 0;
xmax = 0.01;
Ne = 50;    % Number of elements
meshQ = OneDimLinearMeshGen(xmin, xmax, Ne, order); %Create mesh
%Get Material and source coefficients
bloodflow = true; % no blood flow
meshQ = setMatCoeffVectors(meshQ, bloodflow, order);

%Create Linear mesh
order = 1;  %Linear order
xmin = 0;
xmax = 0.01;
Ne = 50;    % Number of elements
meshL = OneDimLinearMeshGen(xmin, xmax, Ne, order); %Create mesh

%Get Material and source coefficients
bloodflow = true; % no blood flow
meshL = setMatCoeffVectors(meshL, bloodflow, order);


%%%%% TESTS %%%%    

%% Test 1: test symmetry of the matrix for linear mesh
% Test that this matrix is symmetric
tol = 1e-14;
eID=1; %element ID
elemat = ReactLEM(meshL, 2, eID, 1); %Calculate Linear Element Matrix

assert(abs(elemat(1,2) - elemat(2,1)) <= tol && abs(elemat(1,1) - elemat(2,2))  <= tol, ...
    'Local element matrix is not symmetric!')

%% Test 2: Repeat Test 1 for quadratic funciton
% Test that this matrix is symmetric
tol = 1e-14;
eID=1; %element ID
elemat = ReactLEM(meshQ, 3, eID, 2); %Calculate Linear Element Matrix

assert(abs(elemat(1,3) - elemat(3,1)) <= tol && abs(elemat(1,1) - elemat(3,3))  <= tol, ...
    'Local element matrix is not symmetric!')
%% Test 3: test 2 different elements of the same size produce same matrix
%Test that for two elements of an equispaced mesh, as described in the
%lectures, the element matrices calculated are the same
tol = 1e-14;
eID=1; %element ID
meshL = OneDimLinearMeshGen(0,1,6,1);
meshL.RCvec = ones(1, length(meshL.nvec)); %returns random diffusion coefficient between 0 and 4

elemat1 = ReactLEM(meshL, 2, eID, 1);	%calculate element 1 matrix

eID=2; %element ID

elemat2 = ReactLEM(meshL, 2, eID, 1); %calculate element 2 matrix

diff = elemat1 - elemat2;
diffnorm = sum(sum(diff.*diff));
assert(abs(diffnorm) <= tol)

%% Test 4: Repeat Test 3 for Quadratic Basis Functions
%Test that for two elements of an equispaced mesh, as described in the
%lectures, the element matrices calculated are the same
tol = 1e-14;
eID=1; %element ID
meshQ = OneDimLinearMeshGen(0,1,6,2);   %Quadratic mesh
meshQ.RCvec = ones(1, length(meshQ.nvec)); %returns random diffusion coefficient between 0 and 4

elemat1 = ReactLEM(meshQ, 2, eID, 2);	%calculate element 1 matrix

eID=2; %element ID

elemat2 = ReactLEM(meshQ, 2, eID, 2); %calculate element 2 matrix

diff = elemat1 - elemat2;
diffnorm = sum(sum(diff.*diff));
assert(abs(diffnorm) <= tol)

%% Test 5: test that one matrix is evaluted correctly
% % Test that element 1 of the three element mesh problem described 
% in Tutorial 3 Question 2c has the element matrix evaluated correctly
tol = 1e-14;

eID=1; %element ID
mshL = OneDimLinearMeshGen(0,1,6,2);
mshL.RCvec = ones(1, length(mshL.nvec)); %diffusion coefficient

elemat1 = ReactLEM(mshL, 2, 2, 1); %calculate element 1 matrix

elemat2 = [ (1/18), (1/36); (1/36), (1/18)];    % This is the known result
diff = elemat1 - elemat2; %calculate the difference between the two matrices
SummedDiff = sum(sum(diff)); %calculates the sum of the elements in the diff matrix
assert(abs(SummedDiff) <= tol) %Checks the error of the summed differences is below tol

%% Test 6: test main diagonal values are double antidiagonal values for linear
%basis functions
% using random inputs for number of elemnts and lambda
tol = 1e-14;
eID=1; %element ID
ne = randi([4 100], 1,1);     %sets number of elements to random integer between 4 and 100
meshL = OneDimLinearMeshGen(0,1,ne,1);
meshL.RCvec = 5*abs(rand(1,1))*ones(1, length(meshL.nvec)); %returns random diffusion coefficient between 0 and 4

elemat1 = ReactLEM(meshL, 2, eID, 1);
elem = elemat1(1,1)/elemat1(2,1); %Gets ratio between 1,1 and 2,1 elements
assert(abs(elem - 2) <= tol)     %Checks ratio is 2 within the tolerance



