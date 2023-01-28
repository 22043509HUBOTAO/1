plot(data)
ecg=data(2500:4400)
fs=100;%Hz
fc =20;%Hz
forder=3;
[b,a]=butter(forder,fc/(fs/2));
ecg_lf=filtfilt(b,a,data);
figure,plot(ecg_lf)
fc =15;
[b,a]=butter(forder,fc/(fs/2),'high');
ecg_hf=filtfilt(b,a,ecg);
figure,plot(ecg_hf)
a = 8;
b = [-1 0 -2 0 0 0 2 0 1];
ecg_diff=filter(b,a,ecg');%
figure,plot(ecg_diff)
ecg_diff_2 =ecg_diff.^2;
a = 1;
b = ones(1,0.15*500);
ecg_diff_int = filter(b,a,ecg_diff_2);
figure,plot(ecg_diff_int)

fs = 256;
f = 50;
T =0.1;
t =0:1/fs:T;
x=sin(2*pi*f*t);
[Pxx,F]=periodogram(x,rectwin(length(x)),length(x),fs);
figure
plot(F,10*log10(Pxx))

T=2;
t= 0:1/fs:T;
x =sin(2*pi*f*t);
[Pxx,F]=periodogram(x,rectwin(length(x)),length(x),fs);
figure
plot(F,10*log10(Pxx))

T=10;
t= 0:1/fs:T;
x =sin(2*pi*f*t);
[Pxx,F]=periodogram(x,rectwin(length(x)),length(x),fs);
figure
plot(F,10*log10(Pxx))

fs = 256;
%サンプリングレート
f = 50;
%正弦波の周波数
T = 0.1;
%データの持続時間
t = 0:1/fs:T; %サンプルの時刻
x = sin(2*pi*f*t); %正弦波データ
[Pxx,F] = periodogram(x,hamming(length(x)),length(x),fs);
figure
plot(F,10*log10(Pxx))

fs = 256;
%サンプリングレート
f = 50;
%正弦波の周波数
T = 0.1;
%データの持続時間
t = 0:1/fs:T; %サンプルの時刻
x = sin(2*pi*f*t); %正弦波データ
[Pxx,F] = periodogram(x,rectwin(length(x)),length(x),fs);
figure
plot(F,10*log10(Pxx))

x = randn(1024, 1);
fs = 1;
A = [1 2 3 4 3 2 1]';
X = filter(1, A, x);
order = 2;
[Pxx,F] = pyulear(X,order,1024,fs);
plot(F,10*log10(Pxx))
[Pxx,F] = periodogram(X,rectwin(length(X)),length(X),fs);
figure
plot(F,10*log10(Pxx))

M = 50;
num = length(X)-M;
for i=1:M
[A,E]= arburg(X,i); % A:係数, E:残差分散
AIC(i)=num*(log(2*pi)+1)+num*log(E)+2*(i+1);
end
subplot(1,2,1);plot(AIC,'-bo');
xlabel('次数 order');ylabel('AIC');
title('ARモデルの次数とAIC')
[order min_at] = min(AIC);
order = round(min_at);
plot(AIC);
plot(AIC(2:15));

[Pxx,F] = pburg(X,order,1024,fs);
plot(F,10*log10(Pxx))
[Pxx,F] = pyulear(X,order,1024,fs);
figure,plot(F,10*log10(Pxx))
[Pxx,F] = pburg(X,order,1024,fs);
figure,plot(F,10*log10(Pxx))