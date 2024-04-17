clc, clear; format compact; format longG;
%
% Performs PSO. Solves example 7.8 (page 319) of the
% book of Martins and Ning (2021)
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 

% define PSO parameters
ub = [3,3];
lb = [-3,-3];
n_particle = 40;
max_iter = 1000;
options.alfa = [0.8,1.2];
options.beta_max = 2;
options.gamma_max = 2;
options.max_vel = 0.8*min(abs([ub, lb]) );

% define parameters of objective function and its derivative
params = {};
func.params = {};
func.fobj = @bean;
% func.fcon

[xbest, fbest, hist] = PSO(n_particle, max_iter,lb,ub, func, options)


%%%Plot
x1 = linspace(-2.5,2.5);
x2 = linspace(-1,3);
[X1,X2] = meshgrid(x1,x2);

for ii=1:size(X1,1)
    for jj=1:size(X2,2)
        Z(ii,jj) = bean([X1(ii,jj),X2(ii,jj) ],params);
    end
end

f1 = figure;
set(f1, 'units','inches','position',[1,1,9,5]);


k_set = [1,2,3,10,100,500];

for ik = 1:6

    k = k_set(ik); 
    subplot(2,3,ik)
    
    contour(X1,X2,Z); hold on;
    
    scatter(hist(k).pbest(:,1),hist(k).pbest(:,2), ...
        Marker="o",Color='r',MarkerFaceColor='g',MarkerEdgeColor='g'); hold on
    
    scatter(xbest(1),xbest(2),Marker="pentagram", ...
        Color='r',MarkerFaceColor='r',MarkerEdgeColor='r'); hold on
    
    xlabel('x_1'); ylabel('x_2');
    
    xlim([-2.5,2.5]);
    ylim([-1,3]);
    
    title(strcat('k =',num2str(k)) )

end

%%%
