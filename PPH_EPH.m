function [PPH,PESS,EPH,Es,Eu] = PPH_EPH(x,time)
Es0 = 0;
t = length(time);
%% PPH
    for i = 1:length(x)
        if max(x) > 0
           PPH(i) = 1-(mean(x)/max(x));
        else 
           PPH(i) = 1;
        end
    end

%% EPH
    for i = 1:length(x)
        PESS(i) = x(i) - mean(x);
    end
    Es = Es0 + cumtrapz(PESS);
    Eu = max(Es) - min(Es);

    for i = 1:length(x)
        if max(x)>0
            EPH = max(x)/Eu;
        else
            EPH = inf;
        end
    end
end