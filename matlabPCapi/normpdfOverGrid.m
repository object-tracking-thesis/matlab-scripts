function pdf = normpdfOverGrid(mu,P,X,Y) 
    XY = cat(3,X,Y);
    % subtract mu
    XYmmu = bsxfun(@minus,XY,shiftdim(mu(:),-2));

    isigXY = squeeze(sum(bsxfun(@times,shiftdim(inv(P),-2),XYmmu),3));
    XYisXY = sum(isigXY .* XYmmu,3);

    pdf = (1/(2*pi*sqrt(det(P)))) * exp(-0.5 * XYisXY);
end