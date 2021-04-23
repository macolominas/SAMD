% This code generate Figure 10 in the paper "Recursive Diffeomorphism-Based
% Regression for Shape Functions".
%
% By Haizhao Yang

if (1)
    close all;
    clear all;
    
    %% generate signal
    N = 2^12;      %% sampling size
    NM = 0;       %% noise level
    ns = NM*randn(1,N);
    
    x = [0:N-1]/N;
    t = x;
    amp = 0.006;%0.01
    F1 = 120;
    F2 = 185;
    
    sh1 = @(x) gen_shape(x,3);
    sh2 = @(x) gen_shape(x,7);
    
    num_group = 2;
    ins_freq = zeros(num_group,N);
    ins_pre_phase = zeros(num_group,N);
    
    xx = x + amp*sin(2*pi*x);
    ins_freq(1,:) = (1+amp*2*pi*cos(2*pi*x))*F1;
    f1 = zeros(1,N);
    ins_amplt(1,:) = 1+0.4*sin(4*pi*x);
    f1 = ins_amplt(1,:).*sh1(F1*xx); %%%%%ECG gen_shape(FF*xx,2)
    
    yy = x + 2*amp*cos(4*pi*x);
    ins_freq(2,:) = (1-amp*8*pi*sin(4*pi*x))*F2;
    f2 = zeros(1,N);
    ins_amplt(2,:) = 1-0.3*sin(2*pi*x);
    f2 = ins_amplt(2,:).*sh2(F2*yy); %%%%%%
    
    ins_pre_phase(1,:) = (xx)*F1;
    ins_pre_phase(2,:) = (yy)*F2;
    
    
    opt.maxiter = 200;
    opt.eps_error = 1e-6;
    opt.show = 1;
    opt.nknots = 20;% use large nkots for discontinous shapes
    opt.knotremoval_factor= 1.01;
    opt.order = 2;
    opt.eps_diff = opt.eps_error;
    
    if (1) % test example: two components
        fTrue = cell(1,2);
        fTrue{1} = f1;
        fTrue{2} = f2;
        numGroup = 2;
        fff = f1 + f2  + ns;
        shapeTrue = cell(1,2);
        shapeTrue{1} = @(x) sh1(x);
        shapeTrue{2} = @(x) sh2(x);
        [shapeInLoop,comp] = srcIterRegJC(fff,N,numGroup,ins_amplt,ins_pre_phase,opt,fTrue,shapeTrue);
    end
    
    %% post processing to get better components
    %case 2 + true amplitude
    for cntGroup=1:numGroup
        shapepp{cntGroup} = spline(linspace(0,1,length(shapeInLoop{cntGroup}))',shapeInLoop{cntGroup});
        comp{cntGroup} = ppval(shapepp{cntGroup},mod(ins_pre_phase(cntGroup,:),1));
        cnt = cntGroup; cnty = cnt- sign(cnt-1.5  );
        if cnt ==1; fff = f2; else fff=f1; end
        comp{cntGroup} = comp{cntGroup}.* ins_amplt(cntGroup,:)...
            *(fff/(comp{cntGroup}.* ins_amplt(cntGroup,:)));
    end
   
    
    close all;
    
    save ./results/RDBR_fig10.mat;
end

if (1)
    load ./results/RDBR_fig10.mat;
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
        tit = sprintf('./results/RDBR_fig10_%d.fig',i);
        saveas(pic,tit);
        str = sprintf('./results/RDBR_fig10_%d',i);
        print(gcf, '-depsc', str);
    end
end  

