function h=rayleigh_chan(C_samples, Fdmax, N, hindex)
%%The above function simulates the multipath Rayleigh Fading channel using the
%%Jakes's Fading model. It generates a normalised channel response h(t)
%%such that its average power is 1 for simplicity of comparison. The
%%simulator generates C_samples of channel samples at the nomralised
%%Doppler frequency Fdmax, to generate N multipath signals.For frequency
%%selective fading, change the M in gn to hindex. {Refer Notes}

%% Jake's fading Simulator
M=0.5*(0.5*N-1);        %Number of oscillators required to generate N signals
n=1:M;
H=hadamard(M);          %Hadamard Matrix for increasing independence between multipaths
theta=2*pi*n/N;  
fd=Fdmax*cos(theta);    %Doppler Frequency - Normalised
bn=pi*n/(M+1);
gn=(2*pi*(M+1)*n)/(M+1);
A=H(hindex,:);
        
for index=1:C_samples
    %Channel Response
    h_t(index)=abs(sum(A.*cos(2*pi*fd*index+gn)*(cos(bn)+sin(bn)*1i)')); 
end
Pavg=sum(h_t.^2)/(C_samples);   %Average Energy
h=h_t/sqrt(Pavg);               %Normalised Channel Response
