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
    opt.iterStyle = 'JC';
    opt.shapeMethod = 1;
    opt.ampErrBandWidth = 20;
    opt.numSweep = 10;
    
    opt.para.Ls=2000;
    opt.para.bandWidth = 50;
    opt.para.diffeoMethod = 'nufft';
    
    
    NN = 2.^19; % sufficiently large sampling rate
    FFrange = 30:5:50;
    errorRec = cell(length(NN),length(FFrange));
    SL2Rec = cell(length(NN),length(FFrange));
    numGroup = 2;
    for i = 1:length(NN)
        for j = 1:length(FFrange)
            N = NN(i);
            
            x = [0:N-1]/N;
            amp = 0.006;%0.01
            F1 = 2*FFrange(j);
            F2 = 2*FFrange(j);
            
            sh1 = @(x) gen_shape(x,5);
            sh2 = @(x) gen_shape(x,2);
            
            num_group = 2;
            ins_freq = zeros(num_group,N);
            ins_amplt = zeros(num_group,N);
            ins_pre_phase = zeros(num_group,N);
            
            xx = x + amp*sin(2*pi*x);
            ins_freq(1,:) = (1+amp*2*pi*cos(2*pi*x))*F1;
            f1 = zeros(1,N);
            ins_amplt(1,:) = 1;%+0.05*sin(2*pi*x);
            f1 = sh1(F1*xx); %%%%%ECG gen_shape(FF*xx,2)
            
            yy = x + amp*cos(2*pi*x);
            ins_freq(2,:) = (1-amp*2*pi*sin(2*pi*x))*F2;
            f2 = zeros(1,N);
            ins_amplt(2,:) = 1;%+0.05*cos(2*pi*x);
            f2 = sh2(F2*yy); %%%%%%
            
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
            %  [~,~,errorRec{i,j},SL2Rec{i,j}] = srcIterRegJC(fff,N,numGroup,ins_amplt,ins_pre_phase,opt,fTrue,shapeTrue);
            
        end
    end
    
    save ./results/RDSA_fig4_1.mat;
end

if (1)
    load ./results/RDSA_fig4_1.mat;
    pic = figure;
    mx = 9;
    hold on;
    h = zeros(1,length(FFrange));
    for j = 1:length(FFrange)
        diff = errorRec{1,j}(1:end-1)-errorRec{1,j}(2:end);
        pos = find(abs(diff)<1e-4);
        if numel(pos)>0
            pos = pos(1);
        end
        tt = log(abs(errorRec{1,j}(2:mx)-errorRec{1,j}(1:mx-1)));
        ss = tt(1:end-1)-tt(2:end);
        if numel(pos)>0 & pos<=length(ss)
            ss = ss(1:pos);
        end
        h(j) = plot(1:length(ss),ss,'LineWidth',2);
        if numel(pos)>0 & pos<=length(ss)
            plot(pos,ss(pos),'k.','LineWidth',10);
            [j pos]
        end
    end
    legend(h,'N=30','N=35','N=40','N=45','N=50','Location','northeast');
    axis square;
    hold off;
    xlabel('j');ylabel('\eta_j');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    tit = sprintf('./results/RDSA_fig4_%d.fig',1);
    saveas(pic,tit);
    str = sprintf('./results/RDSA_fig4_%d',1);
    print(gcf, '-depsc', str);
end

