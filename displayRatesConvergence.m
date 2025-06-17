function displayRatesConvergence(h,E,pMin,pMax)
logE = log10(E);
logH = log10(h);
    for i=pMin:pMax
       if ( i == pMin)
            fprintf('\n & %.2e & \\bf{---} & ',E(i,1));
       else 
           if(i~=pMax)
           fprintf('%.2e & \\bf{---} & ',E(i,1));
           else
               fprintf('%.2e & \\bf{---} \\\\ \n',E(i,1));
           end
       end
    end
    for j=2:length(h)
        for i=pMin:pMax
            rate = abs((logE(i,j)-logE(i,j-1))/(logH(j)-logH(j-1)));
            if ( i == pMin)
                fprintf('\n & %.2e & \\bf{%.2f} &',abs(E(i,j)),rate);
            else
                if(i~= pMax)
                    fprintf('%.2e & \\bf{%.2f} & ',abs(E(i,j)),rate);
                else
                    fprintf('%.2e & \\bf{%.2f} \\\\ \n',abs(E(i,j)),rate);
                end
            end
        end
    end
    fprintf('\n');
    rates = zeros(1,pMax-pMin+1);
    for k = pMin:pMax
    	[r,~] = polyfit(log10(h),log10(E(k,:)),1);
    	rates(k) = r(1);
    end
    fprintf('\t\t Linear Regresion  of rates\n')
    for i=pMin:pMax
        fprintf('p=%d ( %.2f )\t',i,rates(i));
    end
    fprintf('\n');
end