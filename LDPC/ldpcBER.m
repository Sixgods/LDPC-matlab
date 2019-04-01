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

% Method for creating LDPC matrix (0 = Evencol�о���ֲ�; 1 = Evenboth�о���ֲ�)
%evencol��evenboth��Neal�Ĺ��췽���ɣ��ֱ�ָ���ǡ���������1�ĸ�������ͬ���͡���������1�ĸ�����ͬ��������������1�ĸ�����ͬ����

method = 1;

% Eliminate length-4 cycle��������Ϊ4�ıջ�
noCycle = 1;

% Number of 1s per column for LDPC matrixУ�����ÿ���ж��ٸ�1
onePerCol = 3;

% LDPC matrix reorder strategy (0 = First; 1 = Mincol; 2 = Minprod)
strategy = 0;

% EbN0 in dB     
EbN0 = [0 : 0.5 : 1.5];


% Number of iteration;
iter = 5;

% Number of frame (N bits per frame)һ֡���ֽ���
frame = 10;

% Make the LDPC matrix
H = makeLdpc(M, N, 1, 1, onePerCol);

for i = 1:length(EbN0)
   
   ber1(i) = 0;
   ber2(i) = 0;
   
   % Make random data (0/1) 
   % ��Դ�źţ�����M*FRAME��������ݾ���Ȼ���������룬���һ��0��1����
   dSource = round(rand(M, frame));
   
   for j = 1:frame
      fprintf('Frame : %d\n', j);
      
      % Encoding message ���룬��������������У��λc���µ�У�����H�������ʱ��
      % ֻȡ��dSource������һ��
      [c, newH] = makeParityChk(dSource(:, j), H, strategy);
      % u�Ǵ�У��λ���źţ�����������δ���ƣ�c��У��λ����ȡ����Դ�źŵ�һ�У�����
      % �Ƕ���֮һ
      u = [c; dSource(:, j)];
%----------------------------------------------------------------------
      % BPSK modulation BPSK����
      bpskMod = 2*u - 1;
      % Additional white gaussian noise   AWGN��˹������ N0���߹������ܶȣ�
      % ���������۷���ʱ��0����Ƶ������ʵ���У�û�и�Ƶ�ʵĹ�ϵ��˫�߹������ܶ��ǵ��߹������ܶȵ�һ��
      N0 = 1/(exp(EbN0(i)*log(10)/10));
      %tx�Ǵ������ź�randn����ƽ��ֵ��0��������1����׼����1����̬�ֲ������������
      %�����ھͱ���˷�����N0/2����ֵ��0�ĸ�˹������
      tx = bpskMod + sqrt(N0/2)*randn(size(bpskMod));

%----------------------------------------------------------------------
      % Decoding (select decoding method)     ����
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
