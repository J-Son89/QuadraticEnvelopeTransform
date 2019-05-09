function [corrupted_signal, K, N, total_iterations] = config()
C =[ 0.94737646+0.57831985j,  -0.93834136+0.80723275j, 0.93370208-0.74539468j,  0.04446746-0.61503949j, 0.68611731+0.94646068j, -0.82304000+0.39125639j,0.233645-0.36850826j, -0.99867433+0.95581703j];
Zeta = [0.96607614 + 22.94274544j, 0.96607614-22.94274544j,-0.4+38.39425339j,  -0.4-38.39425339j,0.70902471+12.43801711j,  0.70902471-12.43801711j ,0+ 28.85837933j,  0-28.85837933j];
a=-1;
b=1;
K = 8; %number of exponentials
N=127; %matrix dimension, these should be odd
x = linspace(a,b,2*N+1);
snr =1;
signal = formSignal(x,C,Zeta,K);
noise = awgn(x,snr,'measured');
corrupted_signal = signal;
%+ 1*noise;
%corrupted_signal = sampleSignal().';
total_iterations=1000;
end
