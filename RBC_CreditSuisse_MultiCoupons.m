%Plotting Price of contract as a function of the Coupon Rate from 6% to 21%
C=ones(15,1);
Payoff=ones(15,1);
c=0.06;
i=1;
while(c<0.21)
    C(i)=c;
    Payoff(i)=RBC_CS(c);
    i=i+1;
    c=c+0.01;
end
plot(C, Payoff)