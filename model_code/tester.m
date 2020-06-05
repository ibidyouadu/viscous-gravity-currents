n = 50; %linspace(25,100,4);
lambdavec = linspace(7,8,2);
bta = 1e-4; %10.^linspace(-6,-4,3)
tFinalvec = 10.^linspace(1,3,2);
avec = [0.1, 0.5, 1, 2, 10];
% increase domain of lambda, and startng from large betas, incorporate a
% beta vector

% the following will be length(nvec) x length(lambdavec) x length(betavec)
xN = []; 
err = [];
runTimes = [];
params = [];
total_runs = num2str(length(avec)*length(lambdavec)*length(tFinalvec));
done = 0;
for i=1:length(avec)
    a = avec(i);
    for j=1:length(lambdavec)
        lambda = lambdavec(j);
        for k=1:length(tFinalvec)
            tFinal = tFinalvec(k);
            tic
            gc_molND_nonuniform_s; % plotting is off for these purposes
            runTime = toc;
            TH_sol;
            UA_sol = [y(end,1:n-1) 0];
            astr = num2str(a);
            lambdastr=num2str(lambda);
            tfstr=num2str(tFinal);
            display(['Solving a=',astr,', lambda=',lambdastr,', tFinal=',tfstr]);
            sol_err = norm(h_TH_f-UA_sol,1)/norm(h_TH_f,1);
            
            xN(i,j,k) = y(end,end);
            err(i,j,k) = sol_err;
            runTimes(i,j,k) = runTime;
            done = done + 1;
            done_str = num2str(done);
            display(['Completed ', done_str, ' out of ', total_runs]);
%             dlmwrite('xN.txt', xN);
            %type xN.txt;
%             dlmwrite('err.txt.', err);
            %type err.txt;
%             dlmwrite('runTimes.txt', runTimes);
            %type runTimes.txt;
            save('tester_out_05_22.mat', 'xN', 'err', 'runTimes','params');
        end
    end
end
%save('xN.m', '-ascii', 'xN')
