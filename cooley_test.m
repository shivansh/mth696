clc; clf; close all; clear all;
N = 50;
rows = 5;
columns=10;
data = linspace(0, 2*pi, N);
data = sin(data)+0.1*sin(pi*data);
plot(data);
count=1;
matrix=(reshape(data, 10, 5))'; % Read in row wise
test = matrix;
%Perform column-wise dft
for x=1:columns
    matrix(:,x)=dft(matrix(:,x));
end
%Perform Twiddles
for x=1:columns
    for y=1:rows
        matrix(y,x)=matrix(y,x)*exp((-2*pi*(y-1)*(x-1)/N)*1i);
        % matrix(y,x)=matrix(y,x)*exp((-2*pi*y*x/N)*1i);
    end
end
%Perform FFT on rows
for y=1:rows
    matrix(y,:)=dft(matrix(y,:));
end
out=20*log10(abs(reshape(matrix,1,N)));
figure
subplot(211);
stem(out); title('Using cooley-tukey');
subplot(212);
stem(20*log10(abs(fft(data)))); title('Using fft');