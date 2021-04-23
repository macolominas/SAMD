% This code generate Figure 11 and 12 in the paper "Recursive Diffeomorphism-Based
% Regression for Shape Functions".
%
% By Haizhao Yang

if (1)
    close all;
    clear all;
    
    %% generate signal
    change_amplitude = 1;
    N = 2^12;      %% sampling size
    dif = 0.05;
    NM = 0;       %% noise level
    ns = NM*randn(1,N);
    
    x = [0:N-1]/N;
    t = x;
    amp = 0.01;%0.01
    FF = 160;%160
    
    num_group = 4;
    ins_freq = zeros(num_group,N);
    ins_amplt = zeros(num_group,N);
    ins_pre_phase = zeros(num_group,N);
    
    xx = x + amp*sin(2*pi*x);
    ins_freq(1,:) = (1+amp*2*pi*cos(2*pi*x))*FF;
    f1 = zeros(1,N);
    ins_amplt(1,:) = 1+0.05*sin(4*pi*x);
    f1 = ins_amplt(1,:).*gen_shape2(FF*xx,3); %%%%%ECG gen_shape(FF*xx,2)
    
    yy = x + dif + amp*sin(2*pi*(x+dif));
    ins_freq(2,:) = (1+amp*2*pi*cos(2*pi*(x+dif)))*FF;
    f2 = zeros(1,N);
    ins_amplt(2,:) = 1+0.05*sin(2*pi*x);
    f2 = ins_amplt(2,:).*gen_shape(FF*yy,2); %%%%%%
    
    
    zz = x + dif*2 + amp*sin(2*pi*(x+dif*2));
    ins_freq(3,:) = (1+amp*2*pi*cos(2*pi*x+dif*2))*FF;
    f3 = zeros(1,N);
    ins_amplt(3,:) = 1+0.05*sin(2*pi*x);
    f3 = ins_amplt(3,:).*gen_shape(FF*zz,1);%%%%%%%%%%%%
    
    ww = x + dif*3 + amp*sin(2*pi*(x+dif*3));
    ins_freq(4,:) = (1+amp*2*pi*cos(2*pi*x+dif*3))*FF;
    f4 = zeros(1,N);
    ins_amplt(4,:) = 1+0.05*sin(2*pi*x);
    f4 = ins_amplt(4,:).*gen_shape(FF*ww,3);%%%%%%%%%%
    
    ins_pre_phase(1,:) = (xx)*FF;
    ins_pre_phase(2,:) = (yy)*FF;
    ins_pre_phase(3,:) = (zz)*FF;
    ins_pre_phase(4,:) = (ww)*FF;
    
    
    opt.maxiter = 200;
    opt.eps_error = 1e-6;
    opt.show = 1;
    opt.nknots = 20;
    opt.knotremoval_factor= 1.01;% 1.0001
    opt.order = 3;
    opt.eps_diff = opt.eps_error;
    
    if (1) % test example: four components
        fTrue = cell(1,4);
        fTrue{1} = f1;
        fTrue{2} = f2;
        fTrue{3} = f3;
        fTrue{4} = f4;
        numGroup = 4;
        fff = f1 + f2 +f3 + f4 + ns; ffff = fff;
        shapeTrue = cell(1,4);
        shapeTrue{1} = @(x) gen_shape2(x,3);
        shapeTrue{2} = @(x) gen_shape(x,2);
        shapeTrue{3} = @(x) gen_shape(x,1);
        shapeTrue{4} = @(x) gen_shape(x,3);
        [shapeInLoop,component] = srcIterRegJC(fff,N,numGroup,ins_amplt,ins_pre_phase,opt,fTrue,shapeTrue);
    end
    
    %% post processing to get better components
    for cntGroup=1:numGroup
        shapepp{cntGroup} = spline(linspace(0,1,length(shapeInLoop{cntGroup}))',shapeInLoop{cntGroup});
        comp{cntGroup} = ppval(shapepp{cntGroup},mod(ins_pre_phase(cntGroup,:),1));
        cnt = cntGroup; cnty = cnt- sign(cnt-1.5  );
        if cnt ==1; fff = f2; else fff=f1; end
        comp{cntGroup} = comp{cntGroup}.* ins_amplt(cntGroup,:)...
            *(fff/(comp{cntGroup}.* ins_amplt(cntGroup,:)));
    end
    
    close all;
    
    save ./results/RDBR_fig12.mat;
end

if (1)
    load ./results/RDBR_fig12.mat;
    
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
        tit = sprintf('./results/RDBR_fig12_%d.fig',i);
        saveas(pic,tit);
        str = sprintf('./results/RDBR_fig12_%d',i);
        print(gcf, '-depsc', str);
    end
end


