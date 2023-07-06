function crvs = modified_nrbextract(srf, side)

% This is a version of nrbextract from the nurbs-package, which has been modified 
% to be used in extracting the boundaries of multipatch 
% domains.

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

if (~iscell (srf.knots))
  crvs(1).knots = srf.knots(1);
  crvs(1).coefs = srf.coefs(:,1);
  crvs(2).knots = srf.knots(end);
  crvs(2).coefs = srf.coefs(:,end);
  return
end


for dimension = 1:numel(srf.knots)
  ord = srf.order(dimension);
  if (srf.knots{dimension}(1) ~= srf.knots{dimension}(ord) || ...
      srf.knots{dimension}(end) ~= srf.knots{dimension}(end-ord+1))
    error ('nrbextract: only working for open knot vectors')
  end
end


if (numel (srf.knots) == 2)
  for i = 1:2
    if (i == 1)
      coefs1 = squeeze (srf.coefs(:,1,:));
      coefs2 = squeeze (srf.coefs(:,end,:));
    elseif (i == 2)
      coefs1 = squeeze (srf.coefs(:,:,1));
      coefs2 = squeeze (srf.coefs(:,:,end));
    end
    crvs((i - 1) * 2 + 1) = nrbmak (coefs1, srf.knots{mod (i, 2) + 1});
    crvs((i - 1) * 2 + 2) = nrbmak (coefs2, srf.knots{mod (i, 2) + 1});
  end
elseif (numel (srf.knots) == 3)
  for i = 1:3
    inds = setdiff (1:3, i);
    if (i == 1)
      coefs1 = squeeze (srf.coefs(:,1,:,:));
      coefs2 = squeeze (srf.coefs(:,end,:,:));
    elseif (i == 2)
      coefs1 = squeeze (srf.coefs(:,:,1,:));
      coefs2 = squeeze (srf.coefs(:,:,end,:));
    elseif (i == 3)
      coefs1 = squeeze (srf.coefs(:,:,:,1));
      coefs2 = squeeze (srf.coefs(:,:,:,end));
    end
    crvs((i - 1) * 2 + 1) = nrbmak (coefs1, {srf.knots{inds(1)} srf.knots{inds(2)}});
    crvs((i - 1) * 2 + 2) = nrbmak (coefs2, {srf.knots{inds(1)} srf.knots{inds(2)}});
  end
else
  error ('The entity is not a surface nor a volume')
end

if nargin > 1
  crvs = crvs (side);
end

end
