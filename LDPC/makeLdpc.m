function H = makeLdpc(M, N, method, noCycle, onePerCol)
% Create R = 1/2 low density parity check matrix
%             行重已经确定，是2
%  M        : Number of row
%  N        : Number of column
%  method   : Method for distributing non-zero element
%             {0} Evencol : For each column, place 1s uniformly at random  
%                           在每列随机均匀排列1的位置
%             {1} Evenboth: For each column and row, place 1s uniformly at random
%  noCyle   : Length-4 cycle
%             {0} Ignore (do nothing)
%             {1} Eliminate
%  onePerCol: Number of ones per column
%
%  H        : Low density parity check matrix                   
%
%
% Copyright Bagawan S. Nugroho, 2007 
% http://bsnugroho.googlepages.com


% Number of ones per row (N/M ratio must be 2)
if N/M ~= 2
   fprintf('Code rate must be 1/2\n');
end
%列位置总的1和行的总的1是相等的
onePerRow = (N/M)*onePerCol;

fprintf('Creating LDPC matrix...\n');

switch method
   % Evencol
   case {0}
      % Distribute 1s uniformly at random within column
      %randperm(M)产生1到M的无重复随机数
      %生成1个随机矩阵，行先不变，列就从1到第M列，总共生成N行，然后去其共轭矩阵
      for i = 1:N
         onesInCol(:, i) = randperm(M)';
      end      
      % Create non zero elements (1s) index
      %生成非零元素索引
      %将oneINcol矩阵中依次从第1列的第1行到onepercol行，一直到最后一列元素取出，然后按列读取，组成一个含有N*onepercol个元素的每行1个元素，纵向排列矩阵
      r = reshape(onesInCol(1:onePerCol, :), N*onePerCol, 1);
      %repmat复制并展开命令，是将1:N这个向量复制，并且在行方向展开onepercol甘多次，在列方向展开1次
      tmp = repmat([1:N], onePerCol, 1);
      %C是一个【1，1,2,2,3,3····这样的序列，这就保证了列重】
      c = reshape(tmp, N*onePerCol, 1);
      
      % Create sparse matrix H
      % full 把稀疏矩阵转换为满阵显示 sparse 创建稀疏矩阵
      % 生成M*N阶稀疏矩阵H在以向量r和c为坐标的位置上的对应元素值为向量1的对应值
      H = full(sparse(r, c, 1, M, N));
      
   % Evenboth
   case {1}
      % Distribute 1s uniformly at random within column
      %生成1个随机矩阵，从第i列第1个元素一直到最后，一共有M行，包含M个元素的随机排列
      for i = 1:N
         onesInCol(:, i) = randperm(M)';
      end
        
      % Create non zero elements (1s) index
      r = reshape(onesInCol(1:onePerCol, :), N*onePerCol, 1);
      %repmat复制并展开命令，是将1:N这个向量复制，并且在行方向展开onepercol甘多次，在列方向展开1次
      tmp = repmat([1:N], onePerCol, 1);
      c = reshape(tmp, N*onePerCol, 1);
     
      % Make the number of 1s between rows as uniform as possible     
      
      % Order row index
      %将原来r矩阵中的元素放入ix中，然后将r中的元素进行升序排列放入到r中
      [r, ix] = sort(r);
      
      % Order column index based on row index
      %得到重排后行元素对应的列元素
      for i = 1:N*onePerCol
         cSort(i, :) = c(ix(i));
      end
      
      % Create new row index with uniform weight
      % 生成固定行重的行元素索引
      tmp = repmat([1:M], onePerRow, 1);
      r = reshape(tmp, N*onePerCol, 1);
      
      % Create sparse matrix H
      % Remove any duplicate non zero elements index using logical AND
      S = and(sparse(r, cSort, 1, M, N), ones(M, N));
      H = full(S);     
      
end % switch

% Check rows that have no 1 or only have one 1
for i = 1:M
   
   n = randperm(N);
   % Add two 1s if row has no 1
   if length(find(r == i)) == 0
      H(i, n(1)) = 1;
      H(i, n(2)) = 1;
   % Add one 1 if row has only one 1   
   elseif length(find(r == i)) == 1
      H(i, n(1)) = 1;
   end

end % for i

% If desired, eliminate any length-4 cycle
if noCycle == 1
   
   for i = 1:M
      % Look for pair of row - column
      for j = (i + 1):M         
         w = and(H(i, :), H(j, :));
         c1 = find(w);
         lc = length(c1);
         if lc > 1
                       
            % If found, flip one 1 to 0 in the row with less number of 1s
            if length(find(H(i, :))) < length(find(H(j, :)))
               % Repeat the process until only one column left 
               for cc = 1:lc - 1
                  H(j, c1(cc)) = 0;
               end
            else
               for cc = 1:lc - 1
                  H(i, c1(cc)) = 0;
               end
            end % if            
         
         end % if
      end % for j
   end % for i
  
end % if

fprintf('LDPC matrix is created.\n');
