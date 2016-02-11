function dotProductNorm = dotProductNormalized(P1,P2,P3)
%calculates the normalized dot Product for consecutive points P1, P2, P3

%building normalized vectors
AB = P2-P1;
BC = P3-P2;
ABn = AB/norm(AB);
BCn = BC/norm(BC);

dotProductNorm = dot(ABn,BCn);