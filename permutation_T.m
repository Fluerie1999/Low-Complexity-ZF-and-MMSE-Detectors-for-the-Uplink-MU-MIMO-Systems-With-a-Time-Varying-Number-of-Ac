function T= permutation_T(K, I, inactive_user)
T = eye(K);
index = sort(inactive_user, 'descend');
for k = 1:I
    T([index(k) end-k+1],:) = T([end-k+1 index(k)],:);
end
end

