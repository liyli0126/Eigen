function A = ReadInSparseBlockMatrix(fn)
%% Read in a SparseBlockMatrix data file written by C++ code "LinLite"
% James Liu, ColoState; 2012/07--2016/05

fp = fopen(fn, 'rt');

m = fscanf(fp, '%lf', 1);
n = fscanf(fp, '%lf', 1);
nnz = fscanf(fp, '%lf', 1);

mb = fscanf(fp, '%lf', 1);
nb = fscanf(fp, '%lf', 1);
nnzb = fscanf(fp, '%lf', 1);

A = sparse(m,n);
BlkA = cell(nnzb,1);

BlkRowPos = zeros(nnzb,1);
BlkColPos = zeros(nnzb,1);
BlkBgnRow = zeros(nnzb,1);
BlkBgnCol = zeros(nnzb,1);
BlkRowSz = zeros(nnzb,1);
BlkColSz = zeros(nnzb,1);

for k=1:nnzb
  BlkRowPos(k) = fscanf(fp, '%lf', 1);
  BlkColPos(k) = fscanf(fp, '%lf', 1);
  BlkBgnRow(k) = fscanf(fp, '%lf', 1);
  BlkBgnCol(k) = fscanf(fp, '%lf', 1);
  BlkRowSz(k) = fscanf(fp, '%lf', 1);
  BlkColSz(k) = fscanf(fp, '%lf', 1);
  B = zeros(BlkRowSz(k),BlkRowSz(k));
  for i=1:BlkRowSz(k)
    for j=1:BlkColSz(k)
      B(i,j) = fscanf(fp, '%lf', 1);
    end
  end
  BlkA{k} = B;
end

fclose(fp);

II = zeros(nnz,1);
JJ = zeros(nnz,1);
val = zeros(nnz,1);

loc = 1;
for k=1:nnzb
  B = BlkA{k};
  for i=1:BlkRowSz(k)
    for j=1:BlkColSz(k)
      II(loc) = BlkBgnRow(k) + (i-1);
      JJ(loc) = BlkBgnCol(k) + (j-1);
      val(loc) = B(i,j);
      loc = loc + 1;
    end
  end
end

A = sparse(II,JJ,val);

return;