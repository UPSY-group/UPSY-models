clc
clear all
close all

% Calculate exact solutions for all test matrix equations
% listed in this directory

filenames_base = list_all_test_matrix_equations;

for fi = 1: length( filenames_base)
  filename_base = filenames_base{ fi};
  solve_test_matrix_equation( filename_base);
end

%% Functions

function filenames_base = list_all_test_matrix_equations

filenames_base = [];

henk = dir;

for i = 1: length( henk)
  filename = henk( i).name;
  if endsWith( filename,'_b.nc')
    filename_base = filename( 1:end-5);
    filenames_base{ end+1} = filename_base;
  end
end

end

function solve_test_matrix_equation( filename_base)

disp(['Calculating exact solution for test matrix equation ' filename_base '...'])

% Read matrix A and right-hand side b
filename_A      = [filename_base '.nc'];
filename_b      = [filename_base '_b.nc'];

A = read_CSR_from_NetCDF( filename_A);
b = ncread( filename_b, [filename_base '_b']);

% Display some info about the matrix
is_symmetric = issymmetric( A);

try chol(A)   % See: https://nl.mathworks.com/help/matlab/math/determine-whether-matrix-is-positive-definite.html
  is_spd = true;
catch ME
  is_spd = false;
end

disp(['  Symmetric: ' num2str( is_symmetric) ', symmetric positive-definite: ' num2str( is_spd)])

% Solve Ax=b using Matlab's direct solver
x = A\b;

% Write x to NetCDF
filename_x = [filename_base '_x.nc'];
if exist( filename_x,'file')
  delete( filename_x)
end

var_name = [filename_base '_x'];
f = ncinfo( filename_b);
f.Variables(1).Name = var_name;

ncwriteschema( filename_x, f);
ncwrite( filename_x, var_name, x);

end

function A = read_CSR_from_NetCDF( filename)

% Get matrix size
f = ncinfo( filename);
A.m   = [];
A.n   = [];
A.nnz = [];
for di = 1: length( f.Dimensions)
  if     strcmpi(f.Dimensions(di).Name,'m')
    A.m   = f.Dimensions(di).Length;
  elseif strcmpi(f.Dimensions(di).Name,'n')
    A.n   = f.Dimensions(di).Length;
  elseif strcmpi(f.Dimensions(di).Name,'nnz')
    A.nnz = f.Dimensions(di).Length;
  end
end

% Safety
if isempty(A.m) || isempty(A.n) || isempty(A.nnz)
  error('Couldnt find matrix size in file!')
end

% Read matrix data
A.ptr = ncread( filename, 'ptr');
A.ind = ncread( filename, 'ind');
A.val = ncread( filename, 'val');

% Convert to Matlab sparse matrix
A = CSR_to_sparse( A);

end

function A = CSR_to_sparse( A)

nnz = A.ptr(end)-1;

Ai = zeros( nnz, 1);
Aj = zeros( nnz, 1);
Av = zeros( nnz, 1);

k = 0;
for i = 1: A.m
  for ii = A.ptr( i): A.ptr( i+1)-1
    j = A.ind( ii);
    v = A.val(   ii);
    k = k+1;
    Ai( k) = i;
    Aj( k) = j;
    Av( k) = v;
  end
end

A = sparse(Ai,Aj,Av,A.m,A.n);

end


