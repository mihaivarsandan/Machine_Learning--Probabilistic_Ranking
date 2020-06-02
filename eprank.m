load tennis_data

M = size(W,1);            % number of players
N = size(G,1);            % number of games in 2011 season 

psi = inline('normpdf(x)./normcdf(x)');
lambda = inline('(normpdf(x)./normcdf(x)).*( (normpdf(x)./normcdf(x)) + x)');

pv = 0.5;            % prior skill variance (prior mean is always 0)

% initialize matrices of skill marginals - means and precisions
Ms = nan(M,1); 
Ps = nan(M,1);

% initialize matrices of game to skill messages - means and precisions
Mgs = zeros(N,2); 
Pgs = zeros(N,2);
nr = 1:M;

% allocate matrices of skill to game messages - means and precisions
Msg = nan(N,2); 
Psg = nan(N,2);


Novak_Djokovic = 16;
Rafael_Nadal = 1;
Roger_Federer = 5;
Andy_Murray = 11; 


Novak_Djokovic_mean = [];
Rafael_Nadal_mean = [];
Roger_Federer_mean =[];
Andy_Murray_mean = [];

Novak_Djokovic_var = [];
Rafael_Nadal_var = [];
Roger_Federer_var =[];
Andy_Murray_var = [];

Mean_Matrix=zeros(M,30);
Variance_Matrix=zeros(M,30);


for iter=1:30
  % (1) compute marginal skills
  Ms_old = Ms;
  Ps_old = Ps;
   
  for p=1:M
    % precision first because it is needed for the mean update
    Ps(p) = 1/pv + sum(Pgs(G==p)); 
    Ms(p) = sum(Pgs(G==p).*Mgs(G==p))./Ps(p);
  end
  
  if  (norm(Ms)-norm(Ms_old) < 10^(-5)) &(norm(Ps)-norm(Ps_old)<10^(-5));
      iter
      break
  end

  % (2) compute skill to game messages
  % precision first because it is needed for the mean update
  Psg = Ps(G) - Pgs;
  Msg = (Ps(G).*Ms(G) - Pgs.*Mgs)./Psg;
    
  % (3) compute game to performance messages
  vgt = 1 + sum(1./Psg, 2);
  mgt = Msg(:,1) - Msg(:,2); % player 1 always wins the way we store data
   
  % (4) approximate the marginal on performance differences
  Mt = mgt + sqrt(vgt).*psi(mgt./sqrt(vgt));
  Pt = 1./( vgt.*( 1-lambda(mgt./sqrt(vgt)) ) );
    
  % (5) compute performance to game messages
  ptg = Pt - 1./vgt;
  mtg = (Mt.*Pt - mgt./vgt)./ptg;   
    
  % (6) compute game to skills messages
  Pgs = 1./(1 + repmat(1./ptg,1,2) + 1./Psg(:,[2 1]));
  Mgs = [mtg, -mtg] + Msg(:,[2 1]);
  
  Mean_Matrix(:,iter) = Ms;
  Variance_Matrix(:,iter) = 1./Ps;
  
end
% i=Roger_Federer
% j=Andy_Murray
% 
% normcdf((Ms(i)-Ms(j))/sqrt((1+1/Ps(i)+1/Ps(j))))


% Question E
Mu = Ms;
Var = 1./Ps;
prob_win=zeros(M,1);

for P1=1:M
    for P2=1:M
        if (not(P1==P2))
            prob_win(P1)= prob_win(P1) + normcdf( (Mu(P1)-Mu(P2)) / sqrt(1+Var(P1)+Var(P2)) );
        end
    end
end
av_prob_win = prob_win./(M);

[kk_em, ii_em] = sort(av_prob_win, 'descend');

np = 107;
barh(kk_em(np:-1:1), 'm')
set(gca,'YTickLabel',W(ii_em(np:-1:1)),'YTick',1:np,'FontSize',5)
axis([0 1 0.5 np+0.5])
xlabel('Average Probability of Winning a Match', 'FontSize', 14);
ylabel('Player', 'FontSize', 14);





