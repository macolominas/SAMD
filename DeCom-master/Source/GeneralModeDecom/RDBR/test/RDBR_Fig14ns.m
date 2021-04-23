% This code generate Figure 13, 14, and 15 in the paper "Recursive Diffeomorphism-Based
% Regression for Shape Functions".
%
% By Haizhao Yang

if (1)
    close all;
    clear all;
    
    %% generate signal
    N = 2^16;      %% sampling size
    NM = 1.5;       %% noise level
    ns = NM*randn(1,N);
    
    x = [0:N-1]/N;
    t = x;
    amp = 0.001;%0.01
    F1 = 32;
    F2 = 48;
    
    sh1 = @(x) gen_shape2(x,3);
    sh2 = @(x) gen_shape2(x,2);
    
    num_group = 2;
    ins_freqTrue = zeros(num_group,N);
    ins_ampltTrue = zeros(num_group,N);
    ins_pre_phaseTrue = zeros(num_group,N);
    
    xx = x + amp*sin(2*pi*x);
    ins_freqTrue(1,:) = (1+amp*2*pi*cos(2*pi*x))*F1;
    f1 = zeros(1,N);
    ins_ampltTrue(1,:) = 1+0.05*sin(2*pi*x);
    f1 = ins_ampltTrue(1,:).*sh1(F1*xx); %%%%%ECG gen_shape(FF*xx,2)
    
    yy = x + amp*cos(2*pi*x);
    ins_freqTrue(2,:) = (1-amp*2*pi*sin(2*pi*x))*F2;
    f2 = zeros(1,N);
    ins_ampltTrue(2,:) = 1+0.05*cos(2*pi*x);
    f2 = ins_ampltTrue(2,:).*sh2(F2*yy); %%%%%%
    
    ins_pre_phaseTrue(1,:) = (xx)*F1;
    ins_pre_phaseTrue(2,:) = (yy)*F2;
    
    fff = f1 + f2  + ns;
    
    if (1)
        numGroup = 2;
        opt.eps = 1e-3;
        opt.res = 0.05;
        opt.freq_range = [0 N/512];
        opt.NG = N;
        opt.dt = 1/N;
        opt.t_sc = 0.7;
        opt.NM = NM;
        opt.st = round([400 700]/16/opt.res);
        opt.ed = round([700 950]/16/opt.res);
        opt.num_select = 2;
        opt.red = 8;
        opt.C = 100;
        opt.rad = 1;
        opt.show = 0;
        opt.ex = 14; opt.idx = 6;
        [freq,ins_amplt,ins_pre_phase,comp_select] = insInfo(fff,opt);
        % find peaks
        peaks = zeros(size(comp_select));
        for cnt = 1:numGroup
            figure;plot(comp_select(cnt,:));
            peaks(cnt,:) = peakDetection(comp_select(cnt,:),freq(cnt,:));
            hold on;plot(peaks(cnt,:));
            title('selected one component and its peaks');
        end
        
        % corect phases
        ins_pre_phase = phaseShift(ins_pre_phase,peaks);
        
        if (0) % check results
            for cnt = 1:2
                figure;subplot(1,2,1);plot(ins_pre_phase(cnt,:),'r');hold on; plot(ins_pre_phaseTrue(cnt,:),'b');
                subplot(1,2,2);plot(ins_pre_phase(cnt,:)-ins_pre_phaseTrue(cnt,:));
                figure;subplot(1,2,1);plot(ins_amplt(cnt,:),'r');hold on; plot(ins_ampltTrue(cnt,:),'b');
                subplot(1,2,2);plot(ins_amplt(cnt,:)-ins_ampltTrue(cnt,:));
                figure;subplot(1,2,1);plot(freq(cnt,:),'r');hold on; plot(ins_freqTrue(cnt,:),'b');
                subplot(1,2,2);plot(freq(cnt,:)-ins_freqTrue(cnt,:));
            end
        else
            for cnt = 1:2
                pic = figure;plot(x,freq(cnt,:),'LineWidth',2);hold on; plot(x,ins_freqTrue(cnt,:),'LineWidth',2);
                axis square;legend('Est','True');
                title('remove');
                set(gca, 'FontSize', 16);
                b=get(gca);
                set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
                tit = sprintf('fig14_%d.fig',12+cnt-1);
                saveas(pic,tit);
                str = sprintf('./results/RDBR_fig14_%d',12+cnt-1);
                print(gcf, '-depsc', str);
                pic = figure;plot(x,ins_amplt(cnt,:),'LineWidth',2);hold on; plot(x,ins_ampltTrue(cnt,:),'LineWidth',2);
                axis square;legend('Est','True');
                title('remove');
                set(gca, 'FontSize', 16);
                b=get(gca);
                set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
                tit = sprintf('fig14_%d.fig',14+cnt-1);
                saveas(pic,tit);
                str = sprintf('./results/RDBR_fig14_%d',14+cnt-1);
                print(gcf, '-depsc', str);
            end
        end
    end
    
    opt.maxiter = 200;
    opt.eps_error = 1e-10;
    opt.nknots = 20;
    opt.knotremoval_factor= 1.0001;
    opt.order = 3;
    opt.eps_diff = opt.eps_error;
    
    if (1) % test example: two components
        fTrue = cell(1,2);
        fTrue{1} = f1;
        fTrue{2} = f2;
        numGroup = 2;
        shapeTrue = cell(1,2);
        shapeTrue{1} = @(x) sh1(x);
        shapeTrue{2} = @(x) sh2(x);
        [shapeInLoop,comp,errorRec,SL2Rec,iter,flag] = srcIterRegJC(fff,N,numGroup,ins_amplt,ins_pre_phase,opt,fTrue,shapeTrue);
        % make shape unitary
        for cnt1 = 1:numGroup
            shapeInLoop{cnt1} = shapeInLoop{cnt1}/(norm(shapeInLoop{cnt1})/sqrt(length(shapeInLoop{cnt1})));
        end
    end
    
    %close all;
    
    save ./results/RDBR_fig14ns.mat;
end

if (1)
    load ./results/RDBR_fig14ns.mat;
    close all;
    error = zeros(1,num_group);
    for i = 1:num_group
        [trans_est_shape,min_error]=shape_phase_trans(shapeInLoop{i}.',shapeTrue{i}(linspace(0,1,1000)));
        L = length(trans_est_shape);
        gd = 0:1/L:(1-1/L);
        pic = figure;plot(gd,trans_est_shape,'LineWidth',2);hold on;plot(gd,shapeTrue{i}(linspace(0,1,1000)),'LineWidth',2); hold off;
        legend('Est','True'); title([num2str(i),'th shape with error=',num2str(min_error)]);axis square;
        title('remove');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        tit = sprintf('./results/RDBR_fig14_%d.fig',i+2);
        saveas(pic,tit);
        str = sprintf('./results/RDBR_fig14_%d',i+2);
        print(gcf, '-depsc', str);
    end
end

