% manually finding the limit cycle
% we will compare this output to the output from the bif. diagram

v1_vec = [0:0.5:250];
ttotal = 3000;
mmc = zeros(length(v1_vec),2);

parfor k = 1:length(v1_vec)
    [~,~,mmc(k,:)] = calciumsystem('ttotal',ttotal,'v1',v1_vec(k));
end

%% figure
figure
plot(v1_vec,mmc(:,1),'b','linewidth',2); hold on;
plot(v1_vec,mmc(:,2),'b','linewidth',2,'HandleVisibility','off');

%% find bif results to compare.
bi_LC = load('MatCont7p2\MatCont7p2\bi_LC.mat');
plot(bi_LC.biLC_v1p,bi_LC.biLC_min,'r','linewidth',2);
plot(bi_LC.biLC_v1p,bi_LC.biLC_max,'r','linewidth',2,'HandleVisibility','off');
legend('poor mans','cont. software')


