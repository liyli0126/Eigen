function A = ReadInFullMatrix(fn) 
%% Read in a FullMatrix data file written by C++ code "LinLite"
% James Liu, ColoState; 2012/07--2016/05

fp = fopen(fn, 'rt');

m = fscanf(fp, '%lf', 1);
n = fscanf(fp, '%lf', 1);

A = zeros(m,n);

for i=1:m
  for j=1:n
    tmp = fscanf(fp, '%lf', 1);
    tmp = fscanf(fp, '%lf', 1);
    A(i,j) = fscanf(fp, '%lf', 1);
  end
end

fclose(fp);

return;