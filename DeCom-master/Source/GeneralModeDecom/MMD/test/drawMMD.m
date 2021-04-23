if (1)
    name = {'org_sst','MMD_residual_sst','MMD_residual_sst2','GMD_residual_sst'};
    for cnte = 1:length(name)
        figName = sprintf('%s.fig',name{cnte});
        open(figName);
        title('remove');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        str = sprintf('./results/MMD_fig1_%s',name{cnte});
        print(gcf, '-depsc', str);
      %  command = sprintf('epstopdf %s.eps',str);      system(command);
    end
end
