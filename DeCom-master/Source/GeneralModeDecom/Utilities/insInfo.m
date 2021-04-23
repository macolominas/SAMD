function [freq,ins_amplt,ins_pre_phase,comp_select] = insInfo(sig,opt)
% This code computes several instantaneous properties of the given signal
% sig
%
% By Haizhao Yang

N = numel(sig);
NM = opt.NM;
st = opt.st;
ed = opt.ed;

%synchrosqueezed wave packet transform
eps = opt.eps;
res = opt.res;
freq_range = opt.freq_range;
NG = opt.NG;
dt = opt.dt;
t_sc = opt.t_sc;
x = 0:1/N:(1-1/N);
red = opt.red;
rad = opt.rad;
[T_f coef kk] = ss_wp1_fwd(sig,1,1,1,x,NG,freq_range(2),freq_range(1),rad,1,t_sc,red,eps,res);
loc = find(T_f<NM^2);
T_f(loc) = 0;
if (opt.show)
   % pic = figure;imagesc([0 1],[1 max(ed)]*res,log(1+real(T_f(1:max(ed),:))));title('log(SST) in selected range');
    pic = figure;imagesc(log(1+T_f));title('log(SST) in selected range');
    xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
    axis square; colormap (1-gray);
    if any(strcmp('figName',fieldnames(opt)))
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        str = sprintf('./results/%s',opt.figName);
        tit = [str,'.fig'];
        saveas(pic,tit);
        print(gcf, '-depsc', str);    command = sprintf('epstopdf %s.eps',str);      system(command);
    end
end

%extract isolated components in given ranges
thre = 0;
C = opt.C;
max_num = 1;
pct = 0.01;
num_select = opt.num_select;
cluster = cell(num_select,1);
mode = cell(num_select,1);
T= cell(num_select,1);
ins_amplt = zeros(num_select,N);
freq = zeros(num_select,N);
comp_select = zeros(num_select,N);
ins_pre_phase = zeros(num_select,N);

for cnt = 1:num_select
    T_temp = cell(1,1);
    T_temp{1} = zeros(size(T_f));
    T_temp{1}(st(cnt):ed(cnt),:) = T_f(st(cnt):ed(cnt),:);
    [mode{cnt,1}, ~] = ss_wp1_invT(T_temp, coef, kk, 1, N, freq_range(2), freq_range(1), rad, 1, t_sc, res);
    ins_amplt(cnt,:) = amplt_est(1, mode{cnt,1});
    if (opt.show)
        figure;hold on;plot(mode{cnt,1},'b'); plot(ins_amplt(cnt,:) ,'r');plot(-ins_amplt(cnt,:) ,'r');hold off;head = sprintf('%d recovered mode',cnt);title(head);
    end
    [T{cnt,1},~,~] = ss_wp1_fwd(mode{cnt,1},1,1,1,x,NG,freq_range(2),freq_range(1),rad,1,t_sc,red,eps,res);
    [cluster{cnt,1}, freqTemp, ~, ~] = freq_selection(T{cnt,1}, 1, eps*20, C, max_num, thre, res, pct, freq_range);
    freqTemp = filterFreq(freqTemp,1);
    freq(cnt,:) = DeCom_fft_interpolation(freqTemp,N);
    ins_amplt(cnt,:) = filterAmplt(ins_amplt(cnt,:),max(mean(freq(cnt,:))/10,3));
    comp_select(cnt,:) = mode{cnt,1};
    ins_pre_phase(cnt,:) = pre_phase_est(freq(cnt,:),dt);
end

