function [c, newH] = makeParityChk(dSource, H, strategy)
% Generate parity check vector bases on LDPC matrix H using sparse LU decomposition
%  用LU分解法生成校验向量
%  dSource : Binary source (0/1)
%  H       : LDPC matrix
%  strategy: Strategy for finding the next non-zero diagonal
%  elements
%  寻找下一个非零斜线向量
%            {0} First  : First non-zero found by column search
%                           首先用列寻找非零量
%            {1} Mincol : Minimum number of non-zeros in later columns
%            {2} Minprod: Minimum product of:
%                         - Number of non-zeros its column minus 1
%                         - Number of non-zeros its row minus 1
%           
%  c       : Check bits 校验位
%
%
% Copyright Bagawan S. Nugroho, 2007 
% http://bsnugroho.googlepages.com

% Get the matric dimension
[M, N] = size(H);
% Set a new matrix F for LU decomposition
F = H;
% LU matrices M行N-M列零矩阵
L = zeros(M, N - M);
U = zeros(M, N - M);

% Re-order the M x (N - M) submatrix
for i = 1:M

   % strategy {0 = First; 1 = Mincol; 2 = Minprod}
   switch strategy
      
      % Create diagonally structured matrix using 'First' strategy
      case {0}
         
         % Find non-zero elements (1s) for the diagonal
         % find用于查询非零元素行列的位置，返回该两个列矩阵，分别代表行r与列位置c，
         % 先向下找非零元素，再向右找，重复之后，不再找第i-1列前的元素，以第i列为新的第1列
         [r, c] = find(F(:, i:end));
         
         % Find non-zero diagonal element candidates 找出非零元素中r=i的在整个行矩阵中的行坐标
         rowIndex = find(r == i);
            
         % Find the first non-zero column第i行第1个非零元素在整个F矩阵中的列坐标
         chosenCol = c(rowIndex(1)) + (i - 1);
            
      % Create diagonally structured matrix using 'Mincol' strategy
      case {1}
         
         % Find non-zero elements (1s) for the diagonal
         [r, c] = find(F(:, i:end));
         colWeight = sum(F(:, i:end), 1);
         
         % Find non-zero diagonal element candidates
         rowIndex = find(r == i);
         
         % Find the minimum column weight
         [x, ix] = min(colWeight(c(rowIndex)));
         % Add offset to the chosen row index to match the dimension of the... 
         % original matrix F
         chosenCol = c(rowIndex(ix)) + (i - 1);
             
      % Create diagonally structured matrix using 'Minprod' strategy   
      case {2}
            
         % Find non-zero elements (1s) for the diagonal
         [r, c] = find(F(:, i:end));
         colWeight = sum(F(:, i:end), 1) - 1;
         rowWeight = sum(F(i, :), 2) - 1;
         
         % Find non-zero diagonal element candidates
         rowIndex = find(r == i);
            
         % Find the minimum product
         [x, ix] = min(colWeight(c(rowIndex))*rowWeight);
         % Add offset to the chosen row index to match the dimension of the... 
         % original matrix F
         chosenCol = c(rowIndex(ix)) + (i - 1);
         
      otherwise
         fprintf('Please select columns re-ordering strategy!\n');
      
   end % switch

   % Re-ordering columns of both H and F
   % F与H是已经进行列置换的的矩阵，其前半矩阵对角线元素为1，矩阵任意两行或列交换不影响其值
   tmp1 = F(:, i);
   tmp2 = H(:, i);
   F(:, i) = F(:, chosenCol);
   H(:, i) = H(:, chosenCol);
   F(:, chosenCol) = tmp1;
   H(:, chosenCol) = tmp2;
                     
   % Fill the LU matrices column by column
   % L是下三角矩阵，U是上三角矩阵
   L(i:end, i) = F(i:end, i);
   U(1:i, i) = F(1:i, i);
         
   % There will be no rows operation at the last row
   if i < M           
            
      % Find the later rows with non-zero elements in columni
      % 寻找第i列下从第i+1到最后的非零元素坐标
      [r2, c2] = find(F((i + 1):end, i));          
      % Add current row to the later rows which have a 1 in column i
      % 将F第i行元素取出来，并且向下复制length(r2)这么多次，然后与F中的i+r2行，
      % 也就是其第i列含有1元素的行进行相加，然后除2求余，放回第i+r2行,
      F((i + r2), :) = mod(F((i + r2), :) + repmat(F(i, :), length(r2), 1), 2);
                                                           
   end % if
         
end % for i

% Find B.dsource 乘得结果后，需要把结果化成是0,1的二进制形式，因为
% 在这两个矩阵的元素相乘后，得到的结果是是第1矩阵的行与第2矩阵的列相乘后
% 再相加的结果，因为在matlab中运算的是十进制，而我们运算的是二进制的，所以
% 要对结果取对2求余，把其化回事二进制，我们所有的运算都是模2的运算，matlab
% 这里的运算都是十进制的运算，这是需要注意的
z = mod(H(:, (N - M) + 1:end)*dSource, 2);

% Parity check vector found by solving sparse LU
c = mod(U\(L\z), 2); 

% Return the rearrange H 
newH = H;

fprintf('Message encoded.\n');
