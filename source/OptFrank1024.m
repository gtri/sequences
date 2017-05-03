% Example program to show application of gradient descent to optimization
% of the peak sidelobe of a 1024-element Frank sequence.  The approach is
% as described in a paper submitted to IEEE Transactions on Aerospace and
% Electronic Systems.
tic
% Generate 1024-element Frank sequence
k=0:31;
kk=k'*k;
c=exp(1i*2*pi*kk(:)/32)';
Frank1024_Start = c;
 
N=length(c);
zstart = conv(c,conj(c(end:-1:1)));
 
w=ones(1,2*N-1); % compression length is 2N-1
% Add suppression near the main peak.  Weight can be higher than 20 dB,
% but stability can be a problem when using the ad-hoc scale factor below.
%w(N-150:N+150)=10; % 20 dB null region
w(N)=0; % set peak weight to zero for numeric accuracy
 
for p=[1 35] % iterate for ISL then PSL
    for m=1:5000
        % Compute gradient: From equations and code sample in paper.
        fa=fft(c,2*N-1);
        x = fftshift(ifft(abs(fa).^2));
        g=w.^(2*p) .* x .* abs(x).^(2*p-2);
        gar = fftshift(ifft(fft(g,2*N-1).*fa));
        grad = -4*p*imag(c.*conj([gar(end) gar(1:N-1)]));
        
        % Normalize gradient to reduce numeric accuracy problems.
        grad = grad/max(abs(grad));
        
        % Rough estimate of scale factor for normalized gradient.  Ad-hoc
        % for this example; see [4] for a better equation.
        s=-0.05 / m^0.5;  
        
        % Apply the gradient
        c = c .* exp(1i*s*grad);
    end
end
Frank1024_CPM=c;
zcpm = conv(c,conj(c(end:-1:1)));
psl = max(abs(zcpm(1:N-1)));
fprintf('PSL = %5.3f linear = %5.3f dB\n', psl, 20*log10(psl));
 
% Plot the compressions
figure;
xvals =(-N+1):(N-1);
plot(xvals, 20*log10(abs(zstart)));
hold on;
plot(xvals, 20*log10(abs(zcpm)));
xlim([-N N]);
ylim([-20 65]);
legend('Frank 1024', 'Gradient Descent CPM');
xlabel('Correlation bin');
ylabel('Correlation (dB)');
grid on
toc
