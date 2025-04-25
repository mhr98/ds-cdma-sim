
% DS_CDMA_Sim: This program simulate the BER of DS-CDMA system at various conditions  
% using predefined m-sequence spreading codes.
% Copyright (C) 2025  Mohammad Safa
% GitHub Repository: https://github.com/mhr98/ds-cdma-sim
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%% spreading codes

mseq31=[1,1,0,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0,0,1,0,0,0,0,1,0,1,0,1;
        0,1,0,0,1,1,0,0,0,0,1,0,1,1,0,1,0,1,0,0,0,1,1,1,0,1,1,1,1,1,0;
        0,1,0,0,1,0,0,0,1,0,1,1,1,1,1,0,1,1,0,0,1,1,1,0,0,0,0,1,1,0,1];
    
mseq127=[1,1,0,1,0,0,1,0,1,0,1,0,1,1,0,1,1,1,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,0,1,1,1,1,0,1,0,1,0,0,0,0,0,1,0,1,1,1,1,0,0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,0,0,1,0,1,1,0,0,0,1,1,0,0,1,1,0,1,1,0,1,0,0,0,1,1,1,1,1,0,0,0,0,1,1,1,0,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0,1,1,0,1,0,1;
         0,1,0,0,1,0,0,1,0,1,1,0,1,0,1,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,0,1,1,1,0,1,0,1,1,1,0,0,1,1,1,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,1,1,0,1,0,0,1,1,0,1,0,0,0,1,0,1,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,1,0,1,0,1,1,0,0,0,1,0,0,1,1,1,1,0,0,1;
         0,1,0,1,1,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0,1,1,0,1,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,1,0,0,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0,0,1,1,1,0,1,0,1,0,1,0,0,1,0,1,0,1,1,1,1,0,1,0,0,1,1,0,0,1,1,1,0,0,1,1,0,1,0,1,1,0,0,0,1,0,0,0,1,0,1,1,1,0,1,1,0,0,1,0,0,0,0,1,1,1,1,0];

%% codes characteristics
sp_codes = mseq127; % select the spreading code (mseq31, mseq127)
code1= sp_codes(1,:);
code2= sp_codes(2,:);
code3= sp_codes(3,:);

L=length(code1);

% converting the sequence to bipolar
code1(code1==0)=-1;
code2(code2==0)=-1;
code3(code3==0)=-1;

% results vectors
autoV_c1=[];
autoV_c2=[];
crossV_12=[];

for i= -(L-1):(L-1)
    auto_c1=sum(code1.*circshift(code1',i)');
    auto_c2=sum(code2.*circshift(code2',i)');
    cross_12=sum(code1.*circshift(code2',i)');

    autoV_c1=[autoV_c1 auto_c1];
    autoV_c2=[autoV_c2 auto_c2];
    crossV_12=[crossV_12 cross_12];
end
tau= -(L-1):(L-1);

%ploting the auto and cross correlation
figure(1)
plot(tau,autoV_c1)
grid on
title('autocorrelation code1')
xlabel('tau')
ylabel('Rc1')

figure(2)
plot(tau,autoV_c2)
grid on
title('autocorrelation of code2')
xlabel('tau')
ylabel('Rc2')

figure(3)
plot(tau,crossV_12)
grid on
title('crosscorrelation code1 and code2')
xlabel('tau')
ylabel('Rc12')

%% BER calc (AWGN)

% define SNR vector
EbNodB=0:4;
BE1=[];
BE2=[];
BE3=[];
BE_th=[];

%%%% starting BER calc
for i=1:length(EbNodB)
    
    EbNo=10^(EbNodB(i)/10); % convert from dB to linear
    var_n=1/(EbNo*2*L);  % L from the prossacing gain
    n_sigma = sqrt(var_n);
    
    ne1=0;
    ne2=0;
    ne3=0;
    n_fram=0;
    
    while (ne1+ne2+ne3)/3<100   
      %transmitter
      data1=rand>0.5;
      data2=rand>0.5;
      data3=rand>0.5;
      
      % convert the data to bipolar
      Bipolar_data1=2*(data1)-1;
      Bipolar_data2=2*(data2)-1;
      Bipolar_data3=2*(data3)-1; 
      
      % spreading the data
      s1=rectpulse(Bipolar_data1,L)'.*code1;
      s2=rectpulse(Bipolar_data2,L)'.*code2;
      s3=rectpulse(Bipolar_data3,L)'.*code3;
      
      % reciver
      r=(1/L)*(s1+s2+s3)+n_sigma*randn(1,L);
      
      % despread the data using each user's code and make decision
      z1=(sum(r.*code1))>0;
      z2=(sum(r.*code2))>0;
      z3=(sum(r.*code3))>0;
      
      % evaluate the error
      ne1= ne1+(z1 ~= data1);
      ne2= ne2+(z2 ~= data2);
      ne3= ne3+(z3 ~= data3);   

      n_fram=n_fram+1;
    end
    
    % computing BER for the given EbNodB
    BE1=[BE1 ne1/n_fram];
    BE2=[BE2 ne2/n_fram];
    BE3=[BE3 ne3/n_fram];
    
    % computing the theoretical value
    BE_th=[BE_th qfunc(sqrt(2*EbNo))];
      
end

% plotting the results
figure(4)
semilogy(EbNodB, BE1,':r');
hold on
semilogy(EbNodB,BE_th, '-b');
xlabel('Eb/N0 [dB]');ylabel('BE of user 1');legend('Simulated', 'Theory');grid on;hold off;

figure(5)
semilogy(EbNodB, BE2,':r');
hold on
semilogy(EbNodB,BE_th, '-b');
xlabel('Eb/N0 [dB]');ylabel('BE of user 2');legend('Simulated', 'Theory');grid on;hold off;

figure(6)
semilogy(EbNodB, BE3,':r');
hold on
semilogy(EbNodB,BE_th, '-b');
xlabel('Eb/N0 [dB]');ylabel('BE of user 3');legend('Simulated', 'Theory');grid on;hold off;

%% BER calc (Multipath, conv. reciver)

% define SNR vector
EbNodB=0:8;
BE1=[];
BE2=[];
BE3=[];
BE_th=[];

P = 4; %no. of pathes
m_sigma = 1/sqrt(2*P); % each path's power

%starting BER calc
for i=1:length(EbNodB)
    
    EbNo=10^(EbNodB(i)/10); %convert from dB to linear
    var_n=1/(EbNo*2*L);  %L from the prossacing gain
    n_sigma = sqrt(var_n);
    
    ne1=0;
    ne2=0;
    ne3=0;
    n_fram=0;
    
    while (ne1+ne2+ne3)/3<500
      
      %transmitter
      data1=rand>0.5;
      data2=rand>0.5;
      data3=rand>0.5;
      
      %convert the data to bipolar
      Bipolar_data1=2*(data1)-1;
      Bipolar_data2=2*(data2)-1;
      Bipolar_data3=2*(data3)-1; 
      
      %%spreading the data
      s1=rectpulse(Bipolar_data1,L)'.*code1;
      s2=rectpulse(Bipolar_data2,L)'.*code2;
      s3=rectpulse(Bipolar_data3,L)'.*code3;
      
      Tx = (1/L)*(s1+s2+s3);
      
      %channel
      H = abs(m_sigma*(randn(1,P) + 1i*randn(1,P)));
      
      %reciver
      Rx = conv(Tx, H);
      
      r = Rx(1:L) + n_sigma*randn(1,L);
      
      %%% despread the data using each user code and make disicion
      z1=(sum(r.*code1))>0;
      z2=(sum(r.*code2))>0;
      z3=(sum(r.*code3))>0;
      
      %%% evaluate the error
      ne1= ne1+(z1 ~= data1);
      ne2= ne2+(z2 ~= data2);
      ne3= ne3+(z3 ~= data3);   

      n_fram=n_fram+1;
    end
    
    % computin BER for the given EbNodB
    BE1=[BE1 ne1/n_fram];
    BE2=[BE2 ne2/n_fram];
    BE3=[BE3 ne3/n_fram];
    
    % computing the theoritical value
    BE_th=[BE_th 0.5*(1-sqrt(((1/P)*EbNo)/(1+(1/P)*EbNo)))];
      
end

% plotting the results
figure(7)
semilogy(EbNodB, BE1,':r');
hold on
semilogy(EbNodB,BE_th, '-b');
xlabel('Eb/N0 [dB]');ylabel('BE of user 1');legend('Simulated', 'Theory');grid on;hold off;

figure(8)
semilogy(EbNodB, BE2,':r');
hold on
semilogy(EbNodB,BE_th, '-b');
xlabel('Eb/N0 [dB]');ylabel('BE of user 2');legend('Simulated', 'Theory');grid on;hold off;

figure(9)
semilogy(EbNodB, BE3,':r');
hold on
semilogy(EbNodB,BE_th, '-b');
xlabel('Eb/N0 [dB]');ylabel('BE of user 3');legend('Simulated', 'Theory');grid on;hold off;

%% BER calc (Multipath, Rake reciver)
% define SNR vector
EbNodB=0:8;
BE1=[];
BE2=[];
BE3=[];
BE_th_R=[];
BE_th=[];

P = 4; %no. of pathes
m_sigma = 1/sqrt(2*P);

% starting BER calc
for i=1:length(EbNodB)
    
    EbNo=10^(EbNodB(i)/10); %convert from dB to linear
    var_n=1/(EbNo*2);  %L from the prossacing gain
    var_n = var_n/L;
    n_sigma = sqrt(var_n);
    
    ne1=0;
    ne2=0;
    ne3=0;
    n_fram=0;
    
    while (ne1+ne2+ne3)/3<500 
    
      % transmitter
      data1=rand>0.5;
      data2=rand>0.5;
      data3=rand>0.5;
      
      % convert the data to bipolar
      Bipolar_data1=2*(data1)-1;
      Bipolar_data2=2*(data2)-1;
      Bipolar_data3=2*(data3)-1; 
      
      % spreading the data
      s1=rectpulse(Bipolar_data1,L)'.*code1;
      s2=rectpulse(Bipolar_data2,L)'.*code2;
      s3=rectpulse(Bipolar_data3,L)'.*code3;
      
      Tx = (1/L)*(s1+s2+s3);
      
      % channel
      H = abs(m_sigma*(randn(1,P) + 1i*randn(1,P)));
      
      % reciver
      Rx = conv(Tx, H);
      
      r = Rx + n_sigma*randn(1,length(Rx));
      
      % despread the data using each user codes
      
      for k=0:P-1
          r1(k+1) = sum(r((1+k):(L+k)).*code1);
          r2(k+1) = sum(r((1+k):(L+k)).*code2);
          r3(k+1) = sum(r((1+k):(L+k)).*code3);
      end
        
      % combining and decision
      z1=sum(H.*r1)>0;
      z2=sum(H.*r2)>0;
      z3=sum(H.*r3)>0;
      
      % evaluate the error
      ne1= ne1+(z1 ~= data1);
      ne2= ne2+(z2 ~= data2);
      ne3= ne3+(z3 ~= data3);   

      n_fram=n_fram+1;
    end
    
    % computing BER for the given EbNodB
    BE1=[BE1 ne1/n_fram];
    BE2=[BE2 ne2/n_fram];
    BE3=[BE3 ne3/n_fram];
    
    % computing the theoretical value
    BE_th_R=[BE_th_R 0.5*(1-sqrt(((1/P)*EbNo)/(1+(1/P)*EbNo)))];  %single tap theoritical value
    
    mu=sqrt((EbNo/P)/(1+(EbNo/P)));
    pe=0;
    for l=0:P-1
        pe = pe + nchoosek(P-1+l,l)*((1+mu)/2)^l;
    end
    pe=pe*((1-mu)/2)^P;
    BE_th=[BE_th pe];  % theoretical value with multipath diversity
      
end

% plotting the results
figure(10)
semilogy(EbNodB, BE1,':r');
hold on
semilogy(EbNodB,BE_th, '-b');
semilogy(EbNodB,BE_th_R, '-.k');
xlabel('Eb/N0 [dB]');ylabel('BE of user 1');legend('Simulated', 'Theory','Theory single tap');grid on;hold off;

figure(11)
semilogy(EbNodB, BE2,':r');
hold on
semilogy(EbNodB,BE_th, '-b');
semilogy(EbNodB,BE_th_R, '-.k');
xlabel('Eb/N0 [dB]');ylabel('BE of user 2');legend('Simulated', 'Theory','Theory single tap');grid on;hold off;

figure(12)
semilogy(EbNodB, BE3,':r');
hold on
semilogy(EbNodB,BE_th, '-b');
semilogy(EbNodB,BE_th_R, '-.k');
xlabel('Eb/N0 [dB]');ylabel('BE of user 3');legend('Simulated', 'Theory','Theory single tap');grid on;hold off;

%% BER calc (AWGN) single user
% define SNR vector
EbNodB=0:4;
BE1=[];
BE_th=[];

% starting BER calc
for i=1:length(EbNodB)
    
    EbNo=10^(EbNodB(i)/10); %convert from dB to linear
    var_n=1/(EbNo*2*L);  %L from the prossacing gain
    n_sigma = sqrt(var_n);
    
    ne1=0;
    n_fram=0;
    
    while ne1<100
      
      %transmitter
      data1=rand>0.5;
      
      %convert the data to bipolar
      Bipolar_data1=2*(data1)-1;
      
      %%spreading the data
      s1=rectpulse(Bipolar_data1,L)'.*code1;
      
      %reciver
      r=(1/L)*(s1)+n_sigma*randn(1,L);
      
      %%% despread the data using each user code and make disicion
      z1=(sum(r.*code1))>0;
      
      %%% evaluate the error
      ne1= ne1+(z1 ~= data1);  

      n_fram=n_fram+1;
    end
    
    %computin BER for the given EbNodB
    BE1=[BE1 ne1/n_fram];
    
    %computing the theoritical value
    BE_th=[BE_th qfunc(sqrt(2*EbNo))];
      
end

%plotting the results
figure(13)
semilogy(EbNodB, BE1,':r');
hold on
semilogy(EbNodB,BE_th, '-b');
xlabel('Eb/N0 [dB]');ylabel('BE of user 1');legend('Simulated', 'Theory');grid on;hold off;

%% BER calc (Multipath, conv. reciver) single user
% define SNR vector
EbNodB=0:10;
BE1=[];
BE_th=[];

P = 3; %no. of pathes
m_sigma = 1/sqrt(2*P);

% starting BER calc
for i=1:length(EbNodB)
    
    EbNo=10^(EbNodB(i)/10); %convert from dB to linear
    var_n=1/(EbNo*2*L);  %L from the prossacing gain
    n_sigma = sqrt(var_n);
    
    ne1=0;
    n_fram=0;
    
    while ne1<1000
      
      % transmitter
      data1=rand>0.5;
      
      % convert the data to bipolar
      Bipolar_data1=2*(data1)-1;
      
      % spreading the data
      s1=rectpulse(Bipolar_data1,L)'.*code1;
      
      Tx = (1/L)*(s1);
      
      % channel
      H = abs(m_sigma*(randn(1,P) + 1i*randn(1,P)));
      
      % reciver
      Rx = conv(Tx, H);
      
      r = Rx(1:L) + n_sigma*randn(1,L);
      
      % despread the data using each user code and make disicion
      z1=(sum(r.*code1))>0;
      
      % evaluate the error
      ne1= ne1+(z1 ~= data1);

      n_fram=n_fram+1;
    end
    
    % computing BER for the given EbNodB
    BE1=[BE1 ne1/n_fram];
    
    % computing the theoritical value
    BE_th=[BE_th 0.5*(1-sqrt(((1/P)*EbNo)/(1+(1/P)*EbNo)))];
      
end

% plotting the results
figure(14)
semilogy(EbNodB, BE1,':r');
hold on
semilogy(EbNodB,BE_th, '-b');
xlabel('Eb/N0 [dB]');ylabel('BE of user 1');legend('Simulated', 'Theory');grid on;hold off;

%% BER calc (Multipath, Rake reciver) single user
% define SNR vector
EbNodB=0:5;
BE1=[];
BE_th_R=[];
BE_th=[];

P = 3; %no. of pathes
m_sigma = 1/sqrt(2*P);

% starting BER calc
for i=1:length(EbNodB)
    
    EbNo=10^(EbNodB(i)/10); %convert from dB to linear
    var_n=1/(EbNo*2*L);  %L from the prossacing gain
    n_sigma = sqrt(var_n);
    
    ne1=0;
    n_fram=0;
    
    while ne1<500
      
      % transmitter
      data1=rand>0.5;
      
      % convert the data to bipolar
      Bipolar_data1=2*(data1)-1; 
      
      % spreading the data
      s1=rectpulse(Bipolar_data1,L)'.*code1;
      
      Tx = (1/L)*(s1);
      
      % channel
      H = abs(m_sigma*(randn(1,P) + 1i*randn(1,P)));
      
      % reciever
      Rx = conv(Tx, H);
      
      r = Rx + n_sigma*randn(1,length(Rx));
      
      % despread the data using each user codes
      
      for k=0:P-1
          r1(k+1) = sum(r((1+k):(L+k)).*code1);
      end
        
      % combining and decisin
      z1=sum(H.*r1)>0;
      
      % evaluate the error
      ne1= ne1+(z1 ~= data1);  

      n_fram=n_fram+1;
    end
    
    % computing BER for the given EbNodB
    BE1=[BE1 ne1/n_fram];
    
    % computing the theoritical value
    BE_th_R=[BE_th_R 0.5*(1-sqrt(((1/P)*EbNo)/(1+(1/P)*EbNo)))];  %single tap theoritical value
    
    mu=sqrt((EbNo/P)/(1+(EbNo/P)));
    pe=0;
    for l=0:P-1
        pe = pe + nchoosek(P-1+l,l)*((1+mu)/2)^l;
    end
    pe=pe*((1-mu)/2)^P;
    BE_th=[BE_th pe];  %theoritical value with multipath diversity
      
end

% plotting the results
figure(15)
semilogy(EbNodB, BE1,':r');
hold on
semilogy(EbNodB,BE_th, '-b');
semilogy(EbNodB,BE_th_R, '-.k');
xlabel('Eb/N0 [dB]');ylabel('BE of user 1');legend('Simulated', 'Theory','Theory single tap');grid on;hold off;
