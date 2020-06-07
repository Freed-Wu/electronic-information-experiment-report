clear;
T=75e-6;%ʱ��
B=20e6;%����
K=B/T;
Fs=2*B;%����Ƶ��
Ts=1/Fs;
fc=16e9;
N=round(T/Ts);
rep_mat=32;
t=linspace(-0.004*T,0.004*T,0.008*N);
x1=400;%��һ��Ŀ�����
x2=400;%�ڶ���Ŀ�����
v1=40;%��һ��Ŀ���ٶ�
v2=45;%�ڶ���Ŀ���ٶ�
%�����ź�
ze=zeros(1,1488);
St=exp(j*pi*K*t.^2);
St0=[ze,St,ze];
St_1=real(St);
%�ز��ź�
fd1=2*v1*fc/3e8;
fd2=2*v2*fc/3e8;
i=1:N*rep_mat;
Dop1=exp(2*j*pi*fd1*Ts*i);
Dop2=exp(2*j*pi*fd2*Ts*i);
num_x1=round(2*x1/3e8/Ts);
num_x2=round(2*x2/3e8/Ts);
ze_left=zeros(1,1488+num_x1);
ze_right=zeros(1,1488-num_x1);
St1=[ze_left,St,ze_right];
ze2_left=zeros(1,1488+num_x2);
ze2_right=zeros(1,1488-num_x2);
St2=[ze2_left,St,ze2_right];
St3=repmat(St1,1,rep_mat);
St4=repmat(St2,1,rep_mat);
St_he_rep=St3.*Dop1+St4.*Dop2;
snr=30;%�����30
St_he=awgn(St_he_rep,snr,'measured');%�Ӹ�˹������
%��ѹ����
pipei=fliplr(St0);%ƥ���˲���
pipei=conj(pipei);
St_he_maiya=conv(pipei,St_he_rep);%���
St_he_abs=abs(St_he_maiya);
St_he_log=20*log10(St_he_abs);
%����������
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
title( '˫Ŀ�궯̬���', 'Interpreter', 'LaTex');
saveas(gcf,'../fig/two-MTD-matlab.png')

