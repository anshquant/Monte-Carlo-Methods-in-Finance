%General Parameters
T=1;                %time to maturity in years
n=252;              %number of days contract will run
dt=T/n;             %time step
N=10^6;             %number of simulations

%GBM Model Parameters
r=0.0007;   %risk-free interest rate
sigZ=0.73;  %annualized volatility of Zoom
sigN=0.42;  %annualized volatility of Netflix
rho=0.29;   %correlation between stock prices of Zoom and Netflix
zZ=normrnd(0,1,N,n);         %generating correlated Standard Normal 
zN=rho*zZ+sqrt(1-rho^2)*normrnd(0,1,N,n);   % distributions for GBM


%Initializing Variables
Sz0=554.58;         %Initial Level for Zoom
Sn0=331.28;         %Initial Level for Netflix
Bz=305.0190;        %Barrier Level for Zoom
Bn=182.2040;        %Barrier level for Netflix
CRz=3.0186;         %Conversion Ratio for Zoom
CRn=1.8032;         %Conversion Ratio for Netflix
Sz=Sz0*ones(N,n+1); %Stock trajectories of Zoom
Sn=Sn0*ones(N,n+1); %Stock trajectories of Netflix
I=1000;             %Denomination USD 1,000
C=0.11;             %Coupon Rate to be paid semi-annually
C=(C/2)*exp(-r*0.5)+(C/2)*exp(-r*1);  %Discounting Coupon for Present Value
P=C*I*ones(N,1);

%Simulating GBM
for i=1:N
    for j=1:n
        Sz(i,j+1)= Sz(i,j)*exp((r-sigZ^2/2)*dt+sigZ*sqrt(dt)*zZ(i,j));
        Sn(i,j+1)= Sn(i,j)*exp((r-sigN^2/2)*dt+sigN*sqrt(dt)*zN(i,j));
        
    end
    
end

%Determining Scenarios
X=[Sz(:,n+1)/Sz0,Sn(:,n+1)/Sn0];    %Performance of Stock relative to initial level
for i=1:N
    if(all(Sz(i,:)>Bz) & all(Sn(i,:)>Bn)) %#ok<AND2>
        %Scenario (a) i.e Payoff(a)
        P(i)=P(i)+exp(-r*T)*I;
    
    else
        %Scenario (b) i.e Payoff(b)
        [x,k] = min(X(i,:)); %#ok<*ASGLU>
        if(k==1)
            %ZOOM is lower performing
            if(Sz(i,n+1)>Sz0)
                P(i)=P(i)+exp(-r*T)*I;
            else
                P(i)=P(i)+exp(-r*T)*CRz*Sz(i,n+1);
            end
        else
            %Netflix is lower performing
            if(Sn(i,n+1)>Sn0)
                P(i)=P(i)+exp(-r*T)*I;
            else
                P(i)=P(i)+exp(-r*T)*CRn*Sn(i,n+1);
            end
        end
    end
end
disp("Standard Monte Carlo Method- Price & Variance")
[mean(P) var(P)/N]

%%%%%%ANTITHETIC VARIATES METHOD%%%%%%%%%%%
%Re-initializing Variables
Sz1=Sz0*ones(N/2,n+1);
Sz2=Sz0*ones(N/2,n+1);
Sn1=Sn0*ones(N/2,n+1);
Sn2=Sn0*ones(N/2,n+1);
zZ1=normrnd(0,1,N/2,n); %generating Z+
zZ2=-zZ1;               %generating z- for antithetic variates
zN1=rho*zZ1+sqrt(1-rho^2)*normrnd(0,1,N/2,n);
zN2=-zN1;               %doing the same for correlated random variables
I=1000;
c=0.11;
c=(c/2)*exp(-r*0.5)+(c/2)*exp(-r*1);
P=c*I*ones(N/2,2);        

for i=1:N/2
    for j=1:n
        Sz1(i,j+1)= Sz1(i,j)*exp((r-sigZ^2/2)*dt+sigZ*sqrt(dt)*zZ1(i,j));
        Sz2(i,j+1)= Sz2(i,j)*exp((r-sigZ^2/2)*dt+sigZ*sqrt(dt)*zZ2(i,j));
        Sn1(i,j+1)= Sn1(i,j)*exp((r-sigN^2/2)*dt+sigN*sqrt(dt)*zN1(i,j));
        Sn2(i,j+1)= Sn2(i,j)*exp((r-sigN^2/2)*dt+sigN*sqrt(dt)*zN2(i,j));
        
    end
    
end


Sz=Sz1;
Sn=Sn1;
X=[Sz(:,n+1)/Sz0,Sn(:,n+1)/Sn0];

for i=1:N/2
    if(all(Sz(i,:)>Bz) & all(Sn(i,:)>Bn)) %#ok<AND2>
        P(i,1)=P(i,1)+exp(-r*T)*I;
    
    else
        [x,k] = min(X(i,:));
        if(k==1)
            %ZOOM is lower performing
            if(Sz(i,n+1)>Sz0)
                P(i,1)=P(i,1)+exp(-r*T)*I;
            else
                P(i,1)=P(i,1)+exp(-r*T)*CRz*Sz(i,n+1);
            end
        else
            %Netflix is lower performing
            if(Sn(i,n+1)>Sn0)
                P(i,1)=P(i,1)+exp(-r*T)*I;
            else
                P(i,1)=P(i,1)+exp(-r*T)*CRn*Sn(i,n+1);
            end
        end
    end
end

Sz=Sz2;
Sn=Sn2;
X=[Sz(:,n+1)/Sz0,Sn(:,n+1)/Sn0];

for i=1:N/2
    if(all(Sz(i,:)>Bz) & all(Sn(i,:)>Bn)) %#ok<AND2>
        P(i,2)=P(i,2)+exp(-r*T)*I;
    
    else
        [x,k] = min(X(i,:));
        if(k==1)
            %ZOOM is lower performing
            if(Sz(i,n+1)>Sz0)
                P(i,2)=P(i,2)+exp(-r*T)*I;
            else
                P(i,2)=P(i,2)+exp(-r*T)*CRz*Sz(i,n+1);
            end
        else
            %Netflix is lower performing
            if(Sn(i,n+1)>Sn0)
                P(i,2)=P(i,2)+exp(-r*T)*I;
            else
                P(i,2)=P(i,2)+exp(-r*T)*CRn*Sn(i,n+1);
            end
        end
    end
end
Pa = (P(:,1)+P(:,2))/2; %averaging payoffs from z+ and z-
disp("Antithetic Variates Method")
[mean(Pa) 2*var(Pa)/N] %#ok<*NOPTS>

