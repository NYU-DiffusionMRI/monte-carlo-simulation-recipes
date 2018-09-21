clear all

root = '/Users/magda/Desktop/diffusion_simulation/generate_packing/axon_radius_histogram/aboitiz_human_cc/';

fnames = {'genu.mat','ant_body.mat','mid_body.mat','post_body.mat','splenium.mat'};

cc = struct;

figure;

for i = 1:numel(fnames)
    fname = fnames{i};
    cc(i).name = fname;
    load(fullfile(root,fname));
    di = dist(:,1);
    cc(i).di = round(di*10)/10;
    hi = dist(:,2);
    cc(i).hi = hi/sum(hi);
    subplot(5,1,i)
    bar(cc(i).di,cc(i).hi)
    title(cc(i).name)
    xlim([0 9])
    ylim([0 0.3])
    pbaspect([2 1 1])
end

save(fullfile(root,'cc_hist.mat'),'cc')