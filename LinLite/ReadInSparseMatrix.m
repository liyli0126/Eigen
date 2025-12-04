function A = ReadInSparseMatrix(fn)
%% Read in a SparseMatrix data file written by C++ code "LinLite"
% James Liu, ColoState; 2012/07--2016/05

fp = fopen(fn, 'rt');

m = fscanf(fp, '%lf', 1);
n = fscanf(fp, '%lf', 1);
nnz = fscanf(fp, '%lf', 1);

II = zeros(nnz,1);
JJ = zeros(nnz,1);
val = zeros(nnz,1);

for k=1:nnz
  II(k) = fscanf(fp, '%lf', 1);
  JJ(k) = fscanf(fp, '%lf', 1);
  tmp = fscanf(fp, '%lf', 1);
  val(k) = fscanf(fp, '%lf', 1);
end

fclose(fp);

A = sparse(II,JJ,val);

return;