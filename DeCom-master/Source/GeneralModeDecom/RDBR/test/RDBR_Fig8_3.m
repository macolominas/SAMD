% This code generate Figure 8_3 in the paper "Recursive Diffeomorphism-Based
% Regression for Shape Functions".
%
% By Haizhao Yang

if (1)
    
    close all;
    clear all;
    
    
    opt.maxiter = 200;
    opt.eps_error = 1e-13; % sufficiently small tolerance
    opt.show = 0;
    nknots = 20; % sufficiently large number of knots
    opt.knotremoval_factor= 1.01;
    order = 3;
    opt.eps_diff = opt.eps_error;
    options = struct('animation', 0, 'knotremoval_factor',opt.knotremoval_factor);
    
    NN = 2.^(7:16); % sufficiently large sampling rate
    err = zeros(length(NN),4);
    x_reg = 0:1/1000:(1-1/1000);
    for i = 1:length(NN)
        i
        N = NN(i);
        X = (0:1/N:(1-1/N))';
        sh1 = @(x) gen_shape(x,5);
        sh2 = @(x) gen_shape(x,2);
        YTrue1 = sh1(X);
        YTrue2 = sh2(X);
        ep = (rand(N,1)-0.5);
        Y1 = YTrue1 + ep;
        Y2 = YTrue2 + ep;
        [pp] = BSFK(X',Y1,order,nknots,[],options);
        s = ppval(pp,x_reg);
        strue = sh1(x_reg);
        err(i,1) =   (norm(s-strue)/sqrt(length(s)));
        
%         [pp] = BSFK(X',Y1,order,nknots,1/20:1/20:(1-1/20),options);
%         s = ppval(pp,x_reg);
%         err(i,3) =   (norm(s-strue)/sqrt(length(s)));
        
        [pp] = BSFK(X',Y2,order,nknots,[],options);
        s = ppval(pp,x_reg);
        strue = sh2(x_reg);
        err(i,2) =   (norm(s-strue)/sqrt(length(s)));
        
%         [pp] = BSFK(X',Y2,order,nknots,1/20:1/20:(1-1/20),options);
%         s = ppval(pp,x_reg);
%         err(i,4) =   (norm(s-strue)/sqrt(length(s)));
    end
    
   save ./results/RDBR_fig8_3.mat;
end

if (1)
   load ./results/RDBR_fig8_3.mat;
    mm = 10;
    pic = figure;
    hold on;
    h = zeros(1,2);
    xx=log2(NN(1:mm)');yy=log2(err((1:mm),1)+err((1:mm),2));
    h(1) = plot(xx,yy,'o','LineWidth',2);
    A = [xx'*xx,sum(xx);sum(xx),length(xx)]; rhs = [xx'*yy;sum(yy)];
    p = inv(A)*rhs;
    h(2) = plot(xx,polyval(p,xx),'LineWidth',2);
   % h(2) = plot(log2(NN),log2(err(:,3)+err(:,4)),'LineWidth',2);
   temp = sprintf('linear fitting\n slope = %3.2f',p(1));
    legend(h,'raw data',temp,'Location','northeast');
    axis square;
    hold off;
    xlabel('log_2(L)'); ylabel('log_2(regression error)');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    tit = sprintf('./results/RDBR_fig8_%d.fig',3);
    saveas(pic,tit);
    str = sprintf('./results/RDBR_fig8_%d',3);
    print(gcf, '-depsc', str);
end
