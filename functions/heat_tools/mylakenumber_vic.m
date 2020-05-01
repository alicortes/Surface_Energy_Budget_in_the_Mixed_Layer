function [Ln,z_t,St,z_g,drhodz] = mylakenumber_vic(z,A,rho,tau)

% function Ln = lakenumber(depth,area,rho,ustar);  
%
% DESCRIPTION: calculates Lakenumber
%
%       Ln = (
%
% INPUT:
%	z   - depth, positive down, zero at the surface (m)
%   A   - lake area as a function of depth z (m^2) (max. at surface).
%   rho - density of the water column at depths z 
%   tau - wind stress (N/m^2)
%
% OUTPUT:
%	Ln - lakenumber (dimensionless)
%
% REFERENCES:
%	Imberger and Patterson, Adv. Applied Mech., 27:303, 1989
%
% NOTE:  Won't work for non-monotonic hypsographic curves, 
%        e.g. if lake has undersea caves or overhangs.

g = 9.8;

z = z(:)'; A = A(:)'; tau = tau(:);

% Make sure dimensions are right and that z starts at the bottom
si = size(rho);
% Replace Nan's with zeros
for j = 1:si(1)
    for k = 1:si(2)
        if isnan(rho(j,k))
            rho(j,k) = 0;
        end
    end
end
if length(z) ~= si(2) 
    rho = rho';
    si = size(rho);
    if length(z) ~= si(2),
        error('Dimension Mismatch');
    end
end
[z,in] = sort(z); rho = rho(:,in);
%- (b.l.) with this reference, area is maximal at the surface.
A = sort(A,'descend');

% Center of volume of lake
z_g = trapz(z,z.*A) ./ trapz(z,A);

% Center of 'density gradient'
%if size(rho,1)>1,
%[drhodz,x] = gradient(rho,z,1);
%else
drhodz = gradient(rho,z,1);
%end
z_t = trapz(z,z(ones(si(1),1),:).*drhodz,2) ./ trapz(z,drhodz,2);

% Overall Stability
St = trapz(z,(z(ones(si(1),1),:)-z_g).*A(ones(si(1),1),:).*rho,2);

% Lake Number
Ln = g*St.*z_t ./ (tau*max(A).^(3/2)*z_g);

dum = find(Ln < 1E-4);
Ln(dum) = 1E-4;

return

