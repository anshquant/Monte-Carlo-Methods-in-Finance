%The function is a copy of the standard MC applied in the
%RBC_CreditSuisse.m file and the function is being called in the
%RBC_CreditSuisse_MultiCoupons file.
function Payoff=RBC_CS(C)
%General Parameters
T=1;                %time to maturity in years
n=252;              %number of days contract will run
dt=T/n;             %time step
N=10^5;             %number of simulations

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
Payoff=mean(P);