% This example shows how Octave's nurbs package can be used to create geometries. 
% It requires the nurbs package to be installed. For boundaries of multipatch
% domains (in the volume sense), as modified script is required. 

clear all
close all

% Here we load the nurbs-package
pkg load nurbs;

% This constructs the geometry. This geometry, although having a cube topology, cannot be represented exactly by polynomial elements
crv = nrbline([1 0 0],[2 0 0]);
srf = nrbrevolve(crv,[0 0 0],[0 0 1],.5*pi);
vol1 = nrbrevolve(srf,[3 0 0],[ 0 -1 0],.5*pi);
vol2 = nrbrevolve(srf,[-1 0 0],[ 0 -1 0],.5*pi);

% Here, we consider a geometry induced by a single volume patch
vol = vol1;

% This extracts the surface
srfs = nrbextract(vol);

% This extracts topology information. This is not needed by Bembel, but prevents the export routine from throwing warnings
[intrfc, bnd] = nrbmultipatch (srfs);

% Here, we transpose patches, since all normals need to be outward directed
srfs(1) = nrbtransp(srfs(1));
srfs(4) = nrbtransp(srfs(4));
srfs(5) = nrbtransp(srfs(5));

% We visualize the output
figure 1;
hold on
for i = 1:6
	nrbkntplot(srfs(i))
end

% This writes the geometry to file, readable by bembel
nrbexport(srfs,intrfc,bnd,'example_single_volume.dat')


%	The first part of the example ends here. The remainder of the script deals with the case of considering the joint boundary of multiple volume patches.
vol = [vol1,vol2];

% Here, we extract the boundary information of the multiple volume patches
[intrfc, bnd] = nrbmultipatch (vol);

for k = 1:numel(bnd)
	% The modified nrbextract extracts the surfaces of the multipatch volume
  	srfs(k) = modified_nrbextract (vol(bnd(k).patches), bnd(k).faces);
end

% Here, we transpose patches, since all normals need to be outward directed
srfs(3) = nrbtransp(srfs(3));
srfs(4) = nrbtransp(srfs(4));
srfs(6) = nrbtransp(srfs(6));
srfs(7) = nrbtransp(srfs(7));
srfs(10) = nrbtransp(srfs(10));


% We visualize the output
figure 2;
hold on
for i = 1:10
	nrbkntplot(srfs(i))
end

% This extracts topology information. This is not needed by Bembel, but prevents the export routine from throwing warnings
[intrfc, bnd] = nrbmultipatch (srfs);

% This writes the geometry to file, readable by bembel
nrbexport(srfs,intrfc,bnd,'example_multi_volume.dat')





