function b = ReadInVector(fn)
%% Read in a Vector data file written by C++ code "LinLite"
% James Liu, ColoState; 2012/07--2016/05

fp = fopen(fn, 'rt');

m = fscanf(fp, '%lf', 1);

b = zeros(m,1);

for i=1:m
  tmp = fscanf(fp, '%lf', 1);
  b(i) = fscanf(fp, '%lf', 1);
end

fclose(fp);

return;