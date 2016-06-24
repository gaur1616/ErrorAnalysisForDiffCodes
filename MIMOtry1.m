clc;
clear all;
close all;

%% Transceiver Parameters:
Mt=2;                       %No. of Transmitters=2 
T=2;                        %No. of Time Slots=2 
n=10^5;                     %Number of bits to be transferred
rhodB=1:20;                 %SNR Range in dB
rho=10.^(0.1.*rhodB);       %SNR Randge 

%% Codeword Parameters (Cyclic STBC):
u1=1;                       %Optimized value for u1
u2=1;                       %Optimized value for u1
l=0:((Mt*T)-1);             
L=Mt*T;                     %No. of Codewords
theta=2*pi/L;               
V=[(exp(1j*u1*theta)) 0; 0 (exp(1j*u2*theta))]; %Generator Matrix for Cyclic STBC (C1, codeword corresponding to l=1)

%% Message to be transferred for the MIMO simulation:
msg=randi([0 max(l)], n, 1); 

%% Channel attributes (For Doppler effect consideration):
speedMPH=60;              %speed in miles per hour
speed=speedMPH*0.44704;   %speed conversion to m/sec, 1 miles/hour = 0.44704 m/sec
fc= 900e6 ;               %carrier frequency= 900Mhz
B= 30e3;                  %Bandwidth=30kHz
Rs= B;                    %Symbol Rate for Normalisation
c=3e8;                    %speed of light
N=34;                     %Number of Signals
C_samples=n;              %Channel Samples

%% Calculation of Maximum Doppler frequency & normalised Doppler frequency=fdmax*Ts:
wavlen=c/fc;
fdmax=speed/wavlen;       %Maximum Doppler Frequency
Fdmax=fdmax/Rs;           %Normalise Doppler Frequency

%% Flat Fading Channel Response using Jakes Simulator (4 Independent Channels):
for hi=1:L
    h(hi,:)=rayleigh_chan(C_samples, Fdmax, N, hi);
end
hdB=10*log10(h);

%% Channel Response Plot for the first 400 samples:
figure(1);
t=1:400;
plot(t, hdB(1, 1:400), t, hdB(2, 1:400), t, hdB(3, 1:400), t, hdB(4, 1:400));
grid on;
legend('Channel 1','Channel 2','Channel 3','Channel 4');
title('Flat Fading Channel Response (400 Samples)');
xlabel('Number of Samples');
ylabel('Channel Amplitude |h(t)| in dB');

%Initializations:
p_err=zeros(2, length(rhodB));    %Simulated probability of error
recmsg=zeros(n,1);                %Received Message(Decoded)
distance=zeros(1,L);

%% Transmission and Reception for different SNR values, for different number of receivers:
for Mr=1:2
    for index=1:length(rhodB)
        S=sqrt(Mt)*eye(Mt);    %Differential space time Code, ?=0
        Y=ones(T, Mr);         %Received Signal matrix corressponding to msg
        
        %% Transmission of message signal using Differential ST scheme:
        for i=2:n+1
            %Mapping of information bits to corressponding codewords Cl:
            in=msg(i-1);
            Cl=sqrt(Mt)*(V^(in));
            %Differential ST Coding:
            S=(1/sqrt(Mt))*Cl*S;
            %Flat fading channel response:
            H=h(1:(Mt*Mr),i-1);
            H=reshape(H,Mt,Mr);
            %AWGN ~ CN(0,1):
            var=1;
            Noise=sqrt(var)*randn(Mt,Mr) + sqrt(var)*1j*randn(Mt,Mr);
            %Received Signal at the Receiver for the 2 Time Slots:
            y=sqrt(rho(index)/Mt)*(S*H)+Noise;
            Y=horzcat(Y,y);
        end

        %% Maximum likelihood Detection of the received codeword:
        for j=Mr+1:Mr:Mr*(n+1)
          yt=Y(:,j:j+Mr-1);     %Received signal for time ?
          yt_1=Y(:,j-Mr:j-1);   %Received signal for time ?-1
          for k=1:length(l)
             Cd=sqrt(Mt)*(V^(k-1));
             dist= yt-((1/sqrt(Mt))*Cd*yt_1);
             distance(k)=(norm(dist,'fro'))^2; %Maximum Likelihood detection for differential ST scheme
          end
          recmsg((j-1)/Mr)=find(distance==min(distance))-1;
        end

    %Experimental SER Computation from the Simulation:
    no_errors=length(find(msg~=recmsg));
    p_err(Mr,index)=no_errors/n;
    end
end

%% Plot of Symbol Error Rate as a function of SNR: 
figure(2);
semilogy(rhodB, p_err(1,:),rhodB, p_err(2,:));
grid on;
legend('Mr=1','Mr=2');
title('SER vs. SNR curve for Differential codes');
xlabel('SNR(dB)');
ylabel('Symbol Error Rate (SER)');