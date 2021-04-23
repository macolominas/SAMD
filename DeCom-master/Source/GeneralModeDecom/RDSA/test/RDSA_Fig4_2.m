% This code generate Figure 4 in the paper "A fast algorithm for 
% multiresolution mode decomposition".
%
% By Haizhao Yang

if (1)
    
    close all;
    clear all;
    
    
    opt.maxiter = 200;
    opt.eps_error = 1e-13; % sufficiently small tolerance
    opt.show = 0;
    opt.eps_diff = opt.eps_error;
    opt.iterStyle = 'GS';
    opt.shapeMethod = 1;
    opt.ampErrBandWidth = 20;
    opt.numSweep = 10;
    
    opt.para.Ls=2000;
    opt.para.bandWidth = 50;
    opt.para.diffeoMethod = 'nufft';
    
    NN = 2.^(9:18); % sufficiently large sampling rate
    FFrange = 100;
    errorRec = cell(length(NN),length(FFrange));
    SL2Rec = cell(length(NN),length(FFrange));
    flag = cell(length(NN),length(FFrange));
    numGroup = 2;
    for i = 1:length(NN)
        for j = 1:length(FFrange)
            N = NN(i);
            
            x = [0:N-1]/N;
            amp = 0.006;%0.01
            F1 = FFrange(j);
            F2 = FFrange(j);
            
            sh1 = @(x) gen_shape(x,5);
            sh2 = @(x) gen_shape(x,2);
            
            num_group = 2;
            ins_freq = zeros(num_group,N);
            ins_amplt = ones(num_group,N);
            ins_pre_phase = zeros(num_group,N);
            
            xx = x + amp*sin(2*pi*x);
            ins_freq(1,:) = (1+amp*2*pi*cos(2*pi*x))*F1;
            f1 = ins_amplt(1,:).*sh1(F1*xx); %%%%%ECG gen_shape(FF*xx,2)
            
            yy = x + amp*cos(2*pi*x);
            ins_freq(2,:) = (1-amp*2*pi*sin(2*pi*x))*F2;
            f2 = ins_amplt(2,:).*sh2(F2*yy); %%%%%%
            
            ins_pre_phase(1,:) = (xx)*F1;
            ins_pre_phase(2,:) = (yy)*F2;
            
            disp(['  i = ',num2str(i),'    j = ',num2str(j)]);
            
            fTrue = cell(1,2);
            fTrue{1} = f1;
            fTrue{2} = f2;
            fff = f1+ f2;
            shapeTrue = cell(1,2);
            shapeTrue{1} = @(x) sh1(x);
            shapeTrue{2} = @(x) sh2(x);
            [~,~,~,errorRec{i,j},SL2Rec{i,j},~,~] = shapeDiffusion(fff,numGroup,ins_amplt,ins_freq,ins_pre_phase,opt,shapeTrue);
            %  [~,~,errorRec{i,j},SL2Rec{i,j},~,flag{i,j}] = srcIterRegJC(fff,N,numGroup,ins_amplt,ins_pre_phase,opt,fTrue,shapeTrue);
            
        end
    end
    
    save ./results/RDSA_fig4_2.mat;
end

if (1)
    load ./results/RDSA_fig4_2.mat;
    mm = numel(NN);
    pic = figure;
    error = zeros(1,length(NN));
    for i = 1:length(NN)
        error(i) = errorRec{i,1}(end);
    end
    hold on;
    xx=log2(NN(1:mm)');yy=log2(error(1:mm)');
    h(1) = plot(xx,yy,'o','LineWidth',2);
    A = [xx'*xx,sum(xx);sum(xx),length(xx)]; rhs = [xx'*yy;sum(yy)];
    p = inv(A)*rhs;
    h(2) = plot(xx,polyval(p,xx),'LineWidth',2);
    temp = sprintf('linear fitting\n slope = %3.2f',p(1));
    legend(h,'raw data',temp,'Location','northeast');
    axis square;
    hold off;
    xlabel('log_2(L)');ylabel('log_2(\epsilon_1^{(j)})');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    tit = sprintf('./results/RDSA_fig4_%d.fig',2);
    saveas(pic,tit);
    str = sprintf('./results/RDSA_fig4_%d',2);
    print(gcf, '-depsc', str);
end
