% Bit error rate of BPSK modulated LDPC codes under AWGN channel
%
%
% Copyright Bagawan S. Nugroho, 2007 
% http://bsnugroho.googlepages.com

clc; 
clear all;

% LDPC matrix size, rate must be 1/2
% Warning: encoding - decoding can be very long for large LDPC matrix!
M = 1000;
N = 2000;

% Method for creating LDPC matrix (0 = Evencol列均衡分布; 1 = Evenboth行均衡分布)
%evencol和evenboth是Neal的构造方法吧，分别指的是“所有列中1的个数均相同”和“所有列中1的个数相同，并且所有行中1的个数相同”。

method = 1;

% Eliminate length-4 cycle消除长度为4的闭环
noCycle = 1;

% Number of 1s per column for LDPC matrix校验矩阵每列有多少个1
onePerCol = 3;

% LDPC matrix reorder strategy (0 = First; 1 = Mincol; 2 = Minprod)
strategy = 0;

% EbN0 in dB     
EbN0 = [0 : 0.5 : 1.5];


% Number of iteration;
iter = 5;

% Number of frame (N bits per frame)一帧的字节数
frame = 10;

% Make the LDPC matrix
H = makeLdpc(M, N, 1, 1, onePerCol);

for i = 1:length(EbN0)
   
   ber1(i) = 0;
   ber2(i) = 0;
   
   % Make random data (0/1) 
   % 信源信号，生成M*FRAME的随机数据矩阵，然后四舍五入，变成一个0、1矩阵
   dSource = round(rand(M, frame));
   
   for j = 1:frame
      fprintf('Frame : %d\n', j);
      
      % Encoding message 编码，返回两个参数，校验位c和新的校验矩阵H，编码的时候
      % 只取了dSource的其中一列
      [c, newH] = makeParityChk(dSource(:, j), H, strategy);
      % u是带校验位的信号，不带噪声，未调制，c是校验位，再取回信源信号的一列，码率
      % 是二分之一
      u = [c; dSource(:, j)];
%----------------------------------------------------------------------
      % BPSK modulation BPSK调制
      bpskMod = 2*u - 1;
      % Additional white gaussian noise   AWGN高斯白噪声 N0单边功率谱密度，
      % 由于在理论分析时，0是中频，而在实际中，没有负频率的关系，双边功率谱密度是单边功率谱密度的一半
      N0 = 1/(exp(EbN0(i)*log(10)/10));
      %tx是带噪声信号randn产生平均值是0，方差是1，标准差是1的正态分布的随机数矩阵
      %而现在就变成了方差是N0/2，均值是0的高斯白噪声
      tx = bpskMod + sqrt(N0/2)*randn(size(bpskMod));

%----------------------------------------------------------------------
      % Decoding (select decoding method)     译码
      %vhat = decodeProbDomain(tx, H, newN0, iter);
     vhat1 = decodeLogDomain(tx, newH, N0, iter);
      vhat2 = decodeLogDomainSimple(tx, newH, iter);
     % vhat = decodeBitFlip(tx, newH,iter);
   
      % Get bit error rate (for brevity, BER calculation includes parity bits)
      [num1, rat1] = biterr(vhat1', u);
      ber1(i) = (ber1(i) + rat1);
      
      [num2, rat2] = biterr(vhat2', u);
      ber2(i) = (ber2(i) + rat2);
      
      
   end % for j
   
   % Get average of BER
   % ber1(i) = ber1(i)/frame;
   % ber2(i) = ber2(i)/frame;
   
end % for i

% Plot the result
semilogy(EbN0, ber1, '*--')
hold
semilogy(EbN0, ber2, 'o--')
grid on
hold off
