%% fitting example

clear all; close all; clc; % clear workspace
load csi_az.mat; % load data
a = a.*1e4; % rescale to correct units

m = .0281; % or scan m values -- different conditions yield m of .0117, .0332, .0281
H = .4; % or scan H values -- different conditions yield m of .464, .382
for j = 1:length(m);
for i = 1:length(H);
    y = z./(sqrt(2).*sqrt(a.*m(j)).^H(i)); % define rescaled max depth
    yp = linspace(min(y),prctile(y,99.9),1000); % define quantiles
    fHy = exp(-yp.^2./2).*yp.^(1/H(i)-2).*exp((H(i)-1/2).*(4.*log(yp)+gfbm(yp))); % evaluate theoretical pdf at percentiles
    fHy = cumsum(fHy,'omitnan')./nansum(fHy); % get theoretical cdf
    [yc xc] = ecdf(y); % get empirical cdf
    yc = interp1(xc(2:end),yc(2:end),yp); % make comparable
    D(i,j) = max(abs(yc-fHy)); % get goodness of fit
    [j./length(m) i./length(H)]
end
end

figure % plot
h = histogram(y,100);
yy = h.Values;
xx = h.BinEdges; xx = (xx(2:end)+xx(1:end-1))./2;
scatter(xx,yy,100,(1/200)*[110 23 78],'filled')
hold on; 
plot(yp(1:end-1),(xx(2)-xx(1))./(yp(2)-yp(1))*8164./sum(diff(fHy(yp<4))).*diff(fHy),'k','linewidth',3)
axis([0 4 0 800])
box on
set(gca,'ticklabelinterpreter','latex','fontsize',16)
ylabel('Counts','interpreter','latex','fontsize',18)
lgnd = legend('Histogram of Lakes','Theoretical Prediction');
set(lgnd,'interpreter','latex','fontsize',18)
xlabel('$y~ \big(= z/\sqrt{2} (\ell \sqrt{a})^H\big)$','interpreter','latex','fontsize',18)
text(3.9,625,'$H = 0.46$','interpreter','latex','fontsize',18,'HorizontalAlignment','right')
text(3.9,575,'$\ell = 0.11$','interpreter','latex','fontsize',18,'HorizontalAlignment','right')
text(3.9,525,'$D_{alt} = 0.012$','interpreter','latex','fontsize',18,'HorizontalAlignment','right')

%% upscaling example

clear all; close all; clc; % clear workspace
load qx.mat; load lagos_area.mat; % load lake data and theoretical distribution
a = a'; % transpose data

for j = 1:1; % for each percentile
r = rand(1,length(a)); % randomly generate data following distribution
for i = 1:length(r);
    [~,prout] = min(abs(r(i)-q));
    y(i) = x(prout);
end
z = y.*0.696137.*a.^.2; % rescale to correct units
f = 10.^(.42+.097.*log10(a)-1.08.*log10(z)); % flux as function of area and depth
f = f./1000.*24.*365.*a; % individual flux
F(j) = sum(f); % total flux
j
end