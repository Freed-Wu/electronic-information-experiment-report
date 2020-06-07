clear;
T=75e-6;
B=20e6;
K=B/T;
Fs=2*B;
Ts=1/Fs;
fc=16e9;
N=round(T/Ts);
rep_mat=32;
t=linspace(-0.004*T,0.004*T,0.008*N);
x=1000;%目标距离
v=40;%目标速度
%发射信号
ze=zeros(1,1488);
St=exp(j*pi*K*t.^2);
St0=[ze,St,ze];
St_1=real(St);
%回波信号
fd=2*v*fc/3e8;
i=1:N*rep_mat;
Dop=exp(2*j*pi*fd*Ts*i);
num_x=round(2*x/3e8/Ts);
ze_left=zeros(1,1488+num_x);
ze_right=zeros(1,1488-num_x);
St=[ze_left,St,ze_right];
St2=repmat(St,1,rep_mat);
St_he_rep=St2.*Dop;
snr=30;%信噪比30
St_he=awgn(St_he_rep,snr,'measured');%加高斯白噪声
%脉压处理
pipei=fliplr(St0);%匹配滤波器
pipei=conj(pipei);
St_he_maiya=conv(pipei,St_he_rep);%卷积
St_he_abs=abs(St_he_maiya);
St_he_log=20*log10(St_he_abs);
%距离门重排
for r=1:rep_mat
	for h=1:N
		St_he_chongpai(h,r)=St_he_maiya((r-1)*N+h);
end
end
St_he_ch=abs(St_he_chongpai);
%FFT
for h=1:N
	St_he_fft(h,:)=abs(fft(St_he_chongpai(h,:)));
end
figure;
mesh(1:rep_mat,1:N,St_he_fft);
xlabel( {'$ v $ / (m / s)'}, 'Interpreter', 'LaTex');
ylabel( {'$ d $ / m'}, 'Interpreter', 'LaTex');
zlabel( {'$ I $ / dB'}, 'Interpreter', 'LaTex');
title( '单目标动态检测', 'Interpreter', 'LaTex');
saveas(gcf,'../fig/one-MTD-matlab.png')

