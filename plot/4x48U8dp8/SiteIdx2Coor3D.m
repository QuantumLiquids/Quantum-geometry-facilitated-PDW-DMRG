function [x,y,z] = SiteIdx2Coor3D(Ly,site_idx)
z = mod(site_idx,2);
site_idx_2d = floor(site_idx/2);
y = mod(site_idx_2d,Ly);
x = floor(site_idx_2d/Ly);
end

