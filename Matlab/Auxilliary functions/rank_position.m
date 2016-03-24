function out = rank_position(vec)
[~,presort]= sort(vec);
[~, out] = sort(presort);
end