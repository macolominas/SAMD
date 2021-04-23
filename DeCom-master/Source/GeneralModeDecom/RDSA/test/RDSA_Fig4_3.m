% This code generate Figure 4 in the paper 
% "A Fast Algorithm for Multiresolution Mode Decomposition"
% by Haizhao Yang.

if (1)
    
    close all;
    clear all;
    
    
    opt.maxiter = 200;
    opt.eps_error = 1e-13; % sufficiently small tolerance
    opt.show = 0;
    opt.nknots = 20; % sufficiently large number of knots
    opt.knotremoval_factor= 1.01;
    opt.order = 3;
    opt.eps_diff = opt.eps_error;
    opt.Ls=200;
    opt.bandWidth = 50;
    
    NN = 2.^(10:19); % sufficiently large sampling rate
    FFrange = 50:20:110;
    errorRec = cell(length(NN),length(FFrange));
    SL2Rec = cell(length(NN),length(FFrange));
    t1 = zeros(length(NN),length(FFrange));
    t2 = zeros(length(NN),length(FFrange));
    flag = cell(length(NN),length(FFrange));
    numGroup = 2;
    
    
    fidFile = sprintf('./results/time.log');
    fid = fopen(fidFile,'a+');
    fprintf(fid,'\\toprule\n');
    for j = 1:length(FFrange)
        for i = 1:length(NN)
            N = NN(i);
            
            x = [0:N-1]/N;
            amp = 0.006;%0.01
            F1 = FFrange(j);
            F2 = FFrange(j);
            
            sh1 = @(x) gen_shape(x,5);
            sh2 = @(x) gen_shape(x,2);
            
            num_group = 2;
            ins_freq = zeros(num_group,N);
            ins_amplt = zeros(num_group,N);
            ins_pre_phase = zeros(num_group,N);
            
            xx = x + amp*sin(2*pi*x);
            ins_freq(1,:) = (1+amp*2*pi*cos(2*pi*x))*F1;
            f1 = zeros(1,N);
            ins_amplt(1,:) = 1+0.05*sin(2*pi*x);
            f1 = ins_amplt(1,:).*sh1(F1*xx); %%%%%ECG gen_shape(FF*xx,2)
            
            yy = x + amp*cos(2*pi*x);
            ins_freq(2,:) = (1-amp*2*pi*sin(2*pi*x))*F2;
            f2 = zeros(1,N);
            ins_amplt(2,:) = 1+0.05*cos(2*pi*x);
            f2 = ins_amplt(2,:).*sh2(F2*yy); %%%%%%
            
            ins_pre_phase(1,:) = (xx)*F1;
            ins_pre_phase(2,:) = (yy)*F2;
            
            disp(['  i = ',num2str(i),'    j = ',num2str(j)]);
            
            fTrue = cell(1,2);
            fTrue{1} = f1;
            fTrue{2} = f2;
            fff = f1+ f2;
            tic;
            shapeRegBSFK(fff,ins_amplt(1,:),ins_pre_phase(1,:),opt);
            t1(i,j) = toc;
            tic;
            shapeRegDSA(fff,ins_amplt(1,:),ins_freq(1,:),ins_pre_phase(1,:),opt);
            t2(i,j) = toc;
            
            if i < length(NN)
                fprintf(fid,'%3.2f & ',t1(i,j)/t2(i,j));
            else
                fprintf(fid,'%3.2f  \\\\ \n',t1(i,j)/t2(i,j));
            end
            
        end
    end
    spd = t2./t1;
    save ./results/RDSA_fig4_3.mat;
end

if (1)
    load ./results/RDSA_fig4_3.mat;
    
    pic = figure;
    plot(log2(NN),log2(t1(:,1))','r-','LineWidth',2);
    hold on;
    plot(log2(NN),log2(t2(:,1))','b-','LineWidth',2);
    legend('RDBR','RDSA');
    axis square; xlabel('log_2(L)');ylabel('log_2(time)');
    str = sprintf('./results/RDSA_fig4_3');
    saveas(pic,[str,'.fig']);
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end
