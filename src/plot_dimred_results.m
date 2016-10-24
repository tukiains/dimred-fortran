
clear all
close all

fpath = '/home/tukiains/Koulujutut/special_assignment/src/output/';

% read files
air = importdata([fpath,'air.dat']);
prof = importdata([fpath,'profile.dat']);
prior = importdata([fpath,'prior.dat']);
dens = importdata([fpath,'dens.dat']);
t = importdata([fpath,'t.dat']);
wn = importdata([fpath,'wn.dat']);
r = importdata([fpath,'r.dat']);
alt = importdata([fpath,'alt.dat']);
d = importdata([fpath,'d.dat']);
P = importdata([fpath,'P.dat']);

% residual and esimated values
ourtheta = r(end-(d-1):end);
r = r(1:length(wn));

% example random draws from the prior
for n=1:500
    theta = randn(d,1);
    draw(n,:) = prior + bsxfun(@times,P*theta,air)/1e9;
end

% show prior, random draws and retrieved profile
figure(1)
hold on
plot(bsxfun(@rdivide,draw,air'),alt,'-','color',[.7 .7 .7],'handlevisibility','off')
plot(prior./air,alt,'r-','linewidth',3)
plot(dens./air,alt,'g-','linewidth',3)
plot(prof./air,alt,'b-','linewidth',3)
h = legend('Prior','Truth','Solution')
set(h,'location','southwest')
legend boxoff
set(gca,'ylim',[0 40])
    
% show simulated transmission and residual of the fit
figure(2)
subplot(2,1,1)
plot(wn,t,'-','color',[.5 .5 .5], 'linewidth',2)
set(gca,'xlim',[min(wn),max(wn)])
ylabel('Transmission')
subplot(2,1,2)
plot(wn,r,'-','color',[.5 .5 .5])
hold on
set(gca,'xlim',[min(wn),max(wn)])
xlabel('Wavenumber')
ylabel('Residual')

