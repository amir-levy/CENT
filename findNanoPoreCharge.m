function [Q,Q0]=findNanoPoreCharge(A_Array,Qext_Array)

if length(A_Array)==1
    A_Array=Qext_Array*0+A_Array;
end

if length(Qext_Array)==1
    Qext_Array=A_Array*0+Qext_Array;
end

for k=1:length(A_Array)
    Qext=Qext_Array(k);
    A=A_Array(k);
    fun=@(Q) Q-A*sinh(Qext-Q);
    initGuess=[A*Qext/(A+1), sign(Qext)*lambertw(A*exp(abs(Qext))/2)];
    %[minFun,initGuessInd]=min(abs(fun(initGuess)));
    %initGuess=initGuess(initGuessInd);
    if isinf(initGuess(2))
        if  abs(fun(initGuess(1)))<Qext
            initGuess=initGuess(1);
        else
            initGuess=Qext;
        end
    else
        funInitGuess=1./abs(fun(initGuess)); 
        funInitGuess(funInitGuess>1e5)=1e5;

        %p=exp(-(abs(fun(initGuess))))/sum(exp(-(abs(fun(initGuess)))));
        p=funInitGuess/sum(funInitGuess);

        initGuess = p*initGuess';
    end
    options = optimoptions('fsolve','Display','none');
    
   %Q(k)=fsolve(fun,initGuess,options);
   Q(k)=initGuess;
    %Q0(k)=initGuess;
%     if mod(k,100)
%         k/length(A_Array)
%     end
end