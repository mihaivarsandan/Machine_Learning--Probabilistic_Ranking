load tennis_data

randn('seed',27); % set the pseudo-random number generator seed

M = size(W,1);            % number of players
N = size(G,1);            % number of games in 2011 season 

pv = 0.5*ones(M,1);           % prior skill variance 
iter = 10000;
w = zeros(M,1);
nr= 1:M;% set skills to prior mean

Novak_Djokovic = 16;
Rafael_Nadal = 1;
Roger_Federer = 5;
Andy_Murray = 11;

Novak_Djokovic_skill =[];
Rafael_Nadal_skill = [];
Roger_Federer_skill =[];
Andy_Murray_skill = [];

total_samples = zeros(M, int16(iter/15)-1);
for i = 1:iter
    i;
  % First, sample performance differences given the skills and outcomes
  
  t = nan(N,1); % contains a t_g variable for each game
  for g = 1:N   % loop over games
    s = w(G(g,1))-w(G(g,2));  % difference in skills
    t(g) = randn()+s;         % performace difference sample
    while t(g) < 0  % rejection sampling: only positive perf diffs accepted
      t(g) = randn()+s; % if rejected, sample again
    end
  end 
 
  
  % Second, jointly sample skills given the performance differences
  
  m = nan(M,1);  % container for the mean of the conditional 
                 % skill distribution given the t_g samples
  for p = 1:M
    m(p) = sum(t.* ((p==G(:, 1)) - (p==G(:, 2))));
  end
  
  iS = zeros(M,M); % container for the sum of precision matrices contributed
                   % by all the games (likelihood terms)
  for g = 1:N
      
      iS(G(g, 1), G(g, 1)) = iS(G(g, 1), G(g, 1)) + 1;
      iS(G(g, 2), G(g, 2)) = iS(G(g, 2), G(g, 2)) + 1;
      iS(G(g, 1), G(g, 2)) = iS(G(g, 1), G(g, 2)) - 1;
      iS(G(g, 2), G(g, 1)) = iS(G(g, 2), G(g, 1)) - 1;% (***TO DO***) build the iS matrix
  end

  iSS = diag(1./pv) + iS; % posterior precision matrix
  % prepare to sample from a multivariate Gaussian
  % Note: inv(M)*z = R\(R'\z) where R = chol(M);
  iR = chol(iSS);  % Cholesky decomposition of the posterior precision matrix
  mu = iR\(iR'\m); % equivalent to inv(iSS)*m but more efficient
    
  % sample from N(mu, inv(iSS))
  w = mu + iR\randn(M,1);
  
  if rem(i,15)==0
      Novak_Djokovic_skill=[Novak_Djokovic_skill;w(Novak_Djokovic)];
      Rafael_Nadal_skill=[Rafael_Nadal_skill;w(Rafael_Nadal)];
      Roger_Federer_skill=[Roger_Federer_skill;w(Roger_Federer)];
      Andy_Murray_skill=[Andy_Murray_skill;w(Andy_Murray)];
      total_samples(:,i/15)=w;
   
  end
  
  %if rem(i,100)==0
%   if i==500
%         [c,lags]=xcov(w,100,'coeff');
%         %plot(lags,c)
%         scatter(nr,w)
%         xlabel('Players')
%         ylabel('Skills')
%end
 drawnow;   
end


%Question D
% Mu_1 = mean(Rafael_Nadal_skill)
% Mu_2 = mean(Andy_Murray_skill)
% Var_1= var(Rafael_Nadal_skill)
% Var_2= var(Andy_Murray_skill)
% 
% 
% 1 - normcdf(0,Mu_1-Mu_2,sqrt(Var_1+Var_2))

%Question E

prob_win = zeros(M,size(total_samples,2));
for P1=1:M
    for P2=1:M
        for iter=1:size(total_samples,2)
            prob_win(P1, iter) =normcdf(total_samples(P1,iter)-total_samples(P2,iter));
        end
    end
end

mean_prob= mean(prob_win,2);

[kk_em, ii_em] = sort(mean_prob, 'descend');

np = 107;
barh(kk_em(np:-1:1), 'r')
set(gca,'YTickLabel',W(ii_em(np:-1:1)),'YTick',1:np,'FontSize',5)
axis([0 1 0.5 np+0.5])
xlabel('Average Probability of Winning a Match', 'FontSize', 14);
ylabel('Player', 'FontSize', 14);

            








