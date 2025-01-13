%%变换光学设计三维平面龙伯透镜
%
%
% programmed by Dr. Yan Zhang, yanzhang_ise@seu.edu.cn
%  Nov. 2, 2022

close all;

% Original Luneburg Lens： R - raidus，eps - permittivity distribution in LL
R = 50; % Unit:mm半径
r = 0 : R; %intial of coordinate r原始坐标r
eps = 2-(r./R).^2;%如果是矩阵a./b就是对应元素想除，如果是两个数就是普通的除法

cratio = 0.2;
b = R * cratio;

x = -R :1: R;
y = (-R :1: R)';
z = (-R :1: R)';

r =  (( x .^ 2 + y .^ 2) .^ 0.5) .* (( x .^ 2 + y .^ 2) .^ 0.5 <= R);
r((( x .^ 2 + y .^ 2) .^ 0.5) > R) = NaN;%球内坐标

eps = 2-(r./R).^2 .* (( x .^ 2 + y .^ 2) .^ 0.5 <= R);
mu = ones( size ( eps ) );%返回eps的数组
mu( isnan ( r ) ) = NaN;%判断数组元素是否是NaN

% F1 = figure;图
% [xx , yy] = meshgrid(x , y);生成网络采样点，列数x行数y
% surface(xx, yy, eps);
% colormap(jet);
% h = colorbar;
% shading interp;
% axis equal;
% axis off;
% set(h, 'fontsize', 15,'Fontname', 'Times New Roman', 'FontWeight', 'bold');
% title('Spherical Luneburg Lens', 'Fontname', 'Times New Roman', 'Fontsize', 22);

% Ellipsoidal transformation

xe = x;
ye = y;
ze = b / R .* z;

TK1 = [1 , 0 , 0 ; 0 , 1, 0; 0, 0, b/R];
T1 = TK1 * TK1.'/ det(TK1);

eps_xx_e = eps .* T1(1, 1);
eps_yy_e = eps .* T1(2, 2);
eps_zz_e = eps .* T1(3, 3);

mu_xx_e = mu .* T1(1, 1);
mu_yy_e = mu .* T1(2, 2);
mu_zz_e = mu .* T1(3, 3);

% re =  (( ye .^ 2 + ze .^ 2) .^ 0.5) .* (( ze .^ 2 ./ (R^2) + ze .^ 2 ./ (b^2)) .^ 0.5 <= 1);
% % re(( ze .^ 2 ./ (R^2) + ze .^ 2 ./ (b^2)) .^ 0.5 > 1) = NaN;

F2 = figure;
[yy , zz] = meshgrid(ye , ze);

subplot(2, 2, 1);
surface(yy, zz, eps_xx_e);
colormap(jet);
h1 = colorbar('position',[0.47 0.68 0.02 0.1],'Ticks', [min(min(eps_xx_e)), max(max(eps_xx_e))]);
shading interp;
axis equal;
axis off;
set(h1, 'fontsize', 15,'Fontname', 'Times New Roman', 'FontWeight', 'bold');
text(-30, b+25, '\epsilon_x_x & \epsilon_y_y', 'Fontname', 'Times New Roman', 'Fontsize', 15);

subplot(2, 2, 2);
surface(yy, zz, eps_zz_e);
colormap(jet);
h2 = colorbar('position',[0.909 0.68 0.02 0.1],'Ticks', [min(min(eps_zz_e)), max(max(eps_zz_e))]);
shading interp;
axis equal;
axis off;
set(h2, 'fontsize', 15,'Fontname', 'Times New Roman', 'FontWeight', 'bold');
text(-30, b+25, '\epsilon_z_z', 'Fontname', 'Times New Roman', 'Fontsize', 15);

subplot(2, 2, 3);
surface(yy, zz, mu_xx_e);
colormap(jet);
h3 = colorbar('position',[0.47 0.17 0.02 0.2]);
caxis([min(min(mu_xx_e))-0.001 max(max(mu_xx_e))+0.001])
shading interp;
axis equal;
axis off;
set(h3, 'fontsize', 15,'Fontname', 'Times New Roman', 'FontWeight', 'bold');
text(-30, b+25, '\mu_x_x & \mu_y_y', 'Fontname', 'Times New Roman', 'Fontsize', 15);

subplot(2, 2, 4);
surface(yy, zz, mu_zz_e);
colormap(jet);
h3 = colorbar('position',[0.909 0.17 0.02 0.2]);
caxis([min(min(mu_zz_e))-0.001 max(max(mu_zz_e))+0.001])
shading interp;
axis equal;
axis off;
set(h3, 'fontsize', 15,'Fontname', 'Times New Roman', 'FontWeight', 'bold');
text(-30, b+25, '\mu_z_z', 'Fontname', 'Times New Roman', 'Fontsize', 15);

sgtitle('Ellipse Luneburg Lens', 'Fontname', 'Times New Roman', 'Fontsize', 22);

close(F2);

% Planar transformantion
xp = x;
yp = y;
zp = b .* z ./ ((R^2 - y.^2).^0.5);
zp(abs(zp)>b) = NaN; 

[yy , zz] = meshgrid(y , z);
[zz] = b .* zz ./ ((R^2 - yy.^2).^0.5);
zz(abs(zz)>b) = NaN; 

%TK2 = [1 , 0 , 0 ; 0 , 1, 0; 0, zp .* yp ./ (R ^ 2 - y.^2), b/R];
%T2 = TK2 * TK2'/ det(TK2) = [(R ^ 2 - yp ^ 2) ^ (1/2) / b, 0, 0; 0, (R ^ 2 - yp ^ 2) ^ (1/2) / b, (yp * zp) / (b * (R ^ 2 - yp ^ 2) ^ (1/2)); 0, (yp * zp) / (b * (R ^ 2 - yp ^ 2) ^ (1/2)), ((R ^ 2- yp ^ 2) ^ (1/2) ) * (b ^ 2 /(R ^2 - yp ^2) + (yp ^ 2 * zp ^ 2) / (R ^ 2 - yp ^ 2) ^ 2 )/b ];
T2_11 = ((R ^ 2 - y .^ 2).^0.5) ./ b;
T2_22 = T2_11;
T2_33 = ((R ^ 2- yy.' .^ 2) .^ (1/2) .* (b ^ 2 ./ (R ^ 2 - yy.' .^ 2) + (yy.' .^ 2) .* (zz.' .^ 2) ./ (R ^ 2 - yy.' .^ 2) .^ 2) )/b; 
T2_23 = ((yy.') .* (zz.')) ./ (b * (R ^ 2 - (yy.') .^ 2) .^ (1/2)); 

eps_xx_p = eps .* T2_11;
eps_yy_p = eps .* T2_22;
eps_zz_p = eps .* T2_33;
eps_yz_p = eps .* T2_23;

eps_zz_p(eps_zz_p>1) = 1;
eps_zz_p(eps_zz_p<=0.2) = 0;
%eps_yz_p(eps_yz_p==inf)=0;

mu_xx_p = mu .* T2_11;
mu_yy_p = mu .* T2_22;
mu_zz_p = mu .* T2_33;
mu_yz_p = mu .* T2_23;
mu_zz_p(mu_zz_p>1) = 1;
%mu_zz_p(mu_zz_p<=0.2) = 0;
%eps_yz_p(eps_yz_p==inf)=0;

r =  (( ye .^ 2 + ze .^ 2) .^ 0.5) .* (( ze .^ 2 ./ (R^2) + ze .^ 2 ./ (b^2)) .^ 0.5 <= 1);
r(( ze .^ 2 ./ (R^2) + ze .^ 2 ./ (b^2)) .^ 0.5 > 1) = NaN;

F3 = figure;


subplot(2, 3, 1);
surface(yy, zz, eps_xx_p.');
colormap(jet);
h = colorbar('position',[0.345 0.71 0.02 0.1],'Ticks', [min(min(eps_xx_p)), max(max(eps_xx_p))]);
shading interp;
axis equal;
axis off;
set(h, 'fontsize', 15,'Fontname', 'Times New Roman', 'FontWeight', 'bold');
text(-58, b+28, '\epsilon_x_x & \epsilon_y_y', 'Fontname', 'Times New Roman', 'Fontsize', 15);

subplot(2, 3, 2);
surface(yy, zz, eps_zz_p.');
colormap(jet);
h = colorbar('position',[0.635 0.71 0.02 0.1],'Ticks', [min(min(eps_zz_p)), max(max(eps_zz_p))]);
shading interp;
axis equal;
axis off;
set(h, 'fontsize', 15,'Fontname', 'Times New Roman', 'FontWeight', 'bold');
text(-28, b+28, '\epsilon_z_z', 'Fontname', 'Times New Roman', 'Fontsize', 15);

subplot(2, 3, 3);
surface(yy, zz, eps_yz_p.');
colormap(jet);
h = colorbar('position',[0.92 0.71 0.02 0.1],'Ticks', [min(min(eps_yz_p)), max(max(eps_yz_p))]);
shading interp;
axis equal;
axis off;
set(h, 'fontsize', 15,'Fontname', 'Times New Roman', 'FontWeight', 'bold');
text(-28, b+28, '\epsilon_y_z', 'Fontname', 'Times New Roman', 'Fontsize', 15);

subplot(2, 3, 4);
surface(yy, zz, mu_xx_p.');
colormap(jet);
h = colorbar('position',[0.345 0.23 0.02 0.1],'Ticks', [min(min(mu_xx_p)), max(max(mu_xx_p))]);
shading interp;
axis equal;
axis off;
set(h, 'fontsize', 15,'Fontname', 'Times New Roman', 'FontWeight', 'bold');
text(-58, b+28, '\mu_x_x & \mu_y_y', 'Fontname', 'Times New Roman', 'Fontsize', 15);

subplot(2, 3, 5);
surface(yy, zz, mu_zz_p.');
colormap(jet);
h = colorbar('position',[0.635 0.23 0.02 0.1],'Ticks', [min(min(mu_zz_p)), max(max(mu_zz_p))]);
shading interp;
axis equal;
axis off;
set(h, 'fontsize', 15,'Fontname', 'Times New Roman', 'FontWeight', 'bold');
text(-28, b+28, '\mu_z_z', 'Fontname', 'Times New Roman', 'Fontsize', 15);

subplot(2, 3, 6);
surface(yy, zz, mu_yz_p.');
colormap(jet);
h = colorbar('position',[0.92 0.23 0.02 0.1],'Ticks', [min(min(mu_yz_p)), max(max(mu_yz_p))]);
shading interp;
axis equal;
axis off;
set(h, 'fontsize', 15,'Fontname', 'Times New Roman', 'FontWeight', 'bold');
text(-28, b+28, '\mu_y_z', 'Fontname', 'Times New Roman', 'Fontsize', 15);

% %------------------------------
% syms z y R b;
% A = [1, 0, 0; 0, 1, 0; 0, z*y/(R^2-y^2), b/sqrt(R^2-y^2)];
% 1 ./ det(A)
% B = A*A.'/det(A)
% %------------------------------