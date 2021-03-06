C=3e8;                            %光速
B=20e6;                           %带宽20MHz
F0=16e9;                          %载波频率16Ghz
Fs=2 * B;                         %采样频率40MHz
Tr=75e-6;                         %时宽75us
tao=0.6e-6;                          %脉宽0.6us

v=100;
lemda=C/F0;
fd=2*v/lemda;                     %多普勒频率

R=1000;
delay=2*R/C;


%产生基带、中频信号
K=B/tao;
n1=(-tao/2:1/Fs:tao/2);
n2=(0:1/Fs:Tr);
sn0=exp(1j*(2*pi*F0*n1+pi*K*n1.^2));%中频
% sn=exp(1j*(pi*K*n1.^2));%基带

%产生回波信号并计算频谱
N1=round(tao*Fs);
N2=round(Tr*Fs);
N0=size(sn0,2);
sn=exp(1j*(2*pi*(fd+F0)*n1+pi*K*n1.^2));
Ndelay=round(delay*Fs);
sr=[zeros(1,Ndelay) sn zeros(1,N2-N0-Ndelay+1)];
figure(1);
dbw_sr=10*log10((sum(abs(sn0.^2)))/length(sn0));
sr=awgn(sr,30,dbw_sr);
subplot(211);plot(n2,real(sr));
xlabel( {'t / s'}, 'Interpreter', 'LaTex');
ylabel( {'I / dB'}, 'Interpreter', 'LaTex');
%title( 'MF time', 'Interpreter', 'LaTex');
title( '中频时域', 'Interpreter', 'LaTex');

spectrum=abs(fft(real(sr)));
spec=[spectrum(end/2:end) spectrum(1:end/2)];
x=(-Fs/2:Fs/(size(sr,2)-1):Fs/2);
subplot(212);plot(x,spec);
xlabel( {'f / Hz'}, 'Interpreter', 'LaTex');
ylabel( {'I / dB'}, 'Interpreter', 'LaTex');
%title( 'MF frequence', 'Interpreter', 'LaTex');
title( '中频频域', 'Interpreter', 'LaTex');
set(gca,'XTick',[-100e6 -75e6 -50e6 -25e6 0 25e6 50e6 75e6 100e6]);
saveas(gcf,'../fig/one-MF-matlab.png')

%下变频后滤波输出基带信号
hunpin=sr.*cos(2*pi*F0*n2);
Hk=filter(10e6,15e6,400e6);        %调用滤波器函数
sjidai=conv(hunpin,Hk);
figure(2);
subplot(211);plot((0:Tr/(size(sjidai,2)-1):Tr),sjidai);
xlabel( {'t / s'}, 'Interpreter', 'LaTex');
ylabel( {'I / dB'}, 'Interpreter', 'LaTex');
%title( 'baseband time', 'Interpreter', 'LaTex');
title( '基带时域', 'Interpreter', 'LaTex');

sjidai_down=sjidai(1:2:end);
spectrum=abs(fft(sjidai_down));
spec=[spectrum(end/2:end) spectrum(1:end/2)];
x=(-Fs/4:Fs/2/(size(spec,2)-1):Fs/4);
subplot(212);plot(x,spec);
xlabel( {'f / Hz'}, 'Interpreter', 'LaTex');
ylabel( {'I / dB'}, 'Interpreter', 'LaTex');
%title( 'baseband frequence', 'Interpreter', 'LaTex');
title( '基带频域', 'Interpreter', 'LaTex');
saveas(gcf,'../fig/one-baseband-matlab.png')

%匹配滤波
ht=exp(-1j*pi*K*n1.^2);
srmy=conv(ht,sjidai);
n3=(0:1/Fs:(size(srmy,2)-1)/Fs);
figure(3);
subplot(211);plot(n3,real(srmy));
xlabel( {'t / s'}, 'Interpreter', 'LaTex');
ylabel( {'I / dB'}, 'Interpreter', 'LaTex');
%title( 'pulse compression time', 'Interpreter', 'LaTex');
title( '脉压时域', 'Interpreter', 'LaTex');

srmy=srmy./max(srmy);
srmy_db=20*log10(abs(srmy));
subplot(212);plot(n3,srmy_db);
xlabel( {'f / Hz'}, 'Interpreter', 'LaTex');
ylabel( {'I / dB'}, 'Interpreter', 'LaTex');
%title( 'pulse compression frequence', 'Interpreter', 'LaTex');
title( '脉压频域', 'Interpreter', 'LaTex');
saveas(gcf,'../fig/one-pulse-matlab.png')

