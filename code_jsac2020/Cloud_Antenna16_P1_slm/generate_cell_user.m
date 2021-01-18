function d_ap_user = generate_cell_user(r_max, r_min, layer, K, fig_en)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate hexagonal cells and locations of users and access points (APs)
% Input: 
%   r_max: radius of cell 
%   r_min: minimum user-AP distance
%   layer: layer of cells
%   K: number of users in each cell
%   fig_en: fingure enable
% Output:
%   d_ap_user: distance between APs and users
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% generate a central cell
n = 6;              % vertices of hexagonal cell
a = 0:2*pi/n:2*pi;  % angles of hexagonal cell
x0 = r_max*cos(a);  % x-axis of central cell
y0 = r_max*sin(a);  % y-axis of central cell


%% generate cells
cell_num = 1 + 3*layer*(layer-1); % number of cells
x = zeros(n+1, cell_num);         % x-axis 
y = zeros(n+1, cell_num);         % y-axis
m = layer - 1;                    % column of cells
layer_count = layer;              % number of cells in each column
s = 1;                            % column index
bs_x = zeros(1, cell_num);        % x location of cell center
bs_y = zeros(1, cell_num);        % y location of cell center
c = 1;                            % index of cell center

for i = -m : m                    % column index
    if mod(layer_count, 2) == 1   
       p = 0; 
       for j = - fix(layer_count/2) : fix(layer_count/2)  
           bs_x(c) = i*1.5*r_max; 
           bs_y(c) = j*sqrt(3)*r_max; 
           c = c + 1;
           x(:, p+s) = x0 + i*1.5*r_max;
           y(:, p+s) = y0 + j*sqrt(3)*r_max;
           p = p + 1;
       end
    else % mod(layer_count, 2) == 0
       p = 0;
       for j = [- fix(layer_count/2) : -1 , 1 : fix(layer_count/2)]
           if j > 0
              t = 0.5;
           else
              t = -0.5;
           end
           bs_x(c) = i*1.5*r_max; 
           bs_y(c) = (j-t)*sqrt(3)*r_max;
           c = c + 1;
           x(:, p+s) = x0 + i*1.5*r_max;
           y(:, p+s) = y0 + (j-t)*sqrt(3)*r_max;
           p = p + 1;
        end
   end
   s = s + layer_count;
   if i < 0 
      layer_count = layer_count + 1;
   else 
      layer_count = layer_count - 1;
   end
end


%% generate the locations of users
user_x = zeros(K, cell_num);
user_y = zeros(K, cell_num);
for c = 1:cell_num
    r_k = sqrt(r_min^2 + (r_max^2-r_min^2)*rand(K,1));
    a_k = 2*pi*rand(K,1);
    user_x(:,c) = bs_x(c) + r_k.*cos(a_k);
    user_y(:,c) = bs_y(c) + r_k.*sin(a_k);
end


%% distance between APs and users, distance between APs and fog function nodes
d_ap_user = zeros(cell_num, K*cell_num);
% AP-user
for c = 1:cell_num
    d_ap_user(c,:) = sqrt((bs_x(c) - reshape(user_x,[K*cell_num,1])).^2 +...
                          (bs_y(c) - reshape(user_y,[K*cell_num,1])).^2);
end



%% figure
if fig_en == 1
   figure;
   plot(x, y, 'k');  % cells
   axis equal;
   hold on;
   scatter(bs_x, bs_y, 'b<', 'linewidth', 1.2)  % AP 
   scatter(reshape(user_x,[K*cell_num,1]), reshape(user_y,[K*cell_num,1]), 'g.');  % users
   set(gcf,'Color', [1 1 1], 'Units', 'Normalized', 'Position', [0.25 0.25 0.5 0.5]);
end
end


