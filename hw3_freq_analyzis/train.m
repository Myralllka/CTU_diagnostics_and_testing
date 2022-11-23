N = 512;  
periods = 2;
dx = 2*periods*pi/N;
x = -periods*pi:dx:periods*pi-dx;
c = zeros(1,13);
y = zeros(size(x));
for n = 1:2:25
    c(n) = 4/pi/n;
    y = y + c(n) * sin(n*x);
end
c = fft(y);
mag = abs(c);
mag = mag/N;
del = 1/periods;
f = -N/periods/2:del:N/periods/2-del;
amp = [mag(1) mag(2:N/2)+mag(end:-1:N/2+2)];
plot(f(N/2+1:end), amp,'.b');
xlim([-30 30]);

figure(2)
plot(f,fftshift(mag),'o');
grid on
xlim([-30 30]);

figure(3)
covN = c/N;
plot(f,real(fftshift(covN)),f,imag(fftshift(covN)),'o');
grid on
xlim([-30 30]);