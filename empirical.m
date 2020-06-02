load('tennis_data.mat')
% [kk, ii] = sort(G, 'descend');

np = 107;

count_won = zeros(np,1);
count_lost = zeros(np,1);
ratio = zeros(np,1);
for i = 1:107
   count_won(i) = sum(G(:,1)==i);
   count_lost(i) = sum(G(:,2)==i);    
   ratio(i) = count_won(i)/ ( count_won(i) + count_lost(i) );
end

[kk_em, ii_em] = sort(ratio, 'descend');

np = 107;
barh(kk_em(np:-1:1), 'b')
set(gca,'YTickLabel',W(ii_em(np:-1:1)),'YTick',1:np,'FontSize',5)
axis([0 1 0.5 np+0.5])
xlabel('Games Won /Total Number of Games Played', 'FontSize', 14);
ylabel('Player', 'FontSize', 14);

int16(10000/15)