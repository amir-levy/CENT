function [sigma_A, surfaceCharge, mu, Eb]=fit1DTransport(cc,j, L, R,initGuess,clr, lB0, maxC0, minC0)

if nargin<6
    clr='k';
end

Na=6.022e23*1e-30;
lB=700;
minC=0.1;
maxC=5000;

if nargin>6
    lB=lB0;
    minC=minC0;
    maxC=maxC0;
end

gamma = lB/L*2*asinh(L/(2*R));
c=2*cc*Na*L*pi*R^2;

if length(initGuess)==2
    myErr=@(sigma_A) myFitErr(c,j,gamma, L*2*pi*R*sigma_A(1),sigma_A(2));
    sigma_A=fminsearch(myErr, initGuess);
    
    myErr=@(sigma_A) myFitErr(c*sigma_A(3),j,gamma, L*2*pi*R*sigma_A(1),sigma_A(2));
    sigma_A=fminsearch(myErr, [sigma_A 0.01]);
else
    myErr=@(sigma_A) myFitErr(c*initGuess(3),j,gamma, L*2*pi*R*sigma_A(1),sigma_A(2));
    sigma_A=fminsearch(myErr, initGuess(1:2));
    sigma_A(3)=initGuess(3); 
end


%cPlot=c; 
cPlot=logspace(log10(minC*2*Na*L*pi*R^2),log10(2*maxC*Na*L*pi*R^2),100); 
sigma=L*2*pi*R*sigma_A(1);

q=findNanoPoreCharge(cPlot*gamma.*sigma_A(3), sigma*gamma)/gamma;
n=sqrt(cPlot.^2.*sigma_A(3).^2+q.^2);

%loglog(cPlot/(2*Na*L*pi*R^2),exp(sigma_A(2))*n,'--','linewidth',1.5,'color',clr);
semilogx(cPlot/(2*Na*L*pi*R^2),exp(sigma_A(2))*n,'--','linewidth',1.5,'color',clr);

e=1.6e-19;
mCSquare2elecPerAngSquare= 6.24e-5;
surfaceCharge=sigma_A(1)/mCSquare2elecPerAngSquare; 
mu=exp(sigma_A(:,2))*(L*1e-10)^2/e*1e-9; %G[nS]=10^9*G[s]=ne\mu*A/L=N/L^2*e\mu
Eb=log(sigma_A(3)); %in units of kT

function err=myFitErr(c,j,gamma, sigma,A)

q=findNanoPoreCharge(c*gamma, sigma*gamma)/gamma;
n=sqrt(c.^2+q.^2);
err=norm(A+log(n)-log(j));






