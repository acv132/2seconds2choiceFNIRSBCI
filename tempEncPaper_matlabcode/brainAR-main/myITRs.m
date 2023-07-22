function cur_ITR = myITRs(res)
    cur_ITR = zeros(10,10);
    for subject = 1:10
        cur_res = res(subject,:);
        for trial_rep = 1:10
            cur_ITR(subject,trial_rep) = calcITR(2, cur_res(trial_rep)/100, trial_rep*10);
        end
    end
end

function ITR = calcITR(N, P, tau)
ITR = (log2(N) +  P*log2(P) + (1-P)*log2((1-P)/(N-1)))*(60/tau);
end
