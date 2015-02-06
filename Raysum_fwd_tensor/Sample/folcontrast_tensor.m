% read synth output and get into same format as rfs
% this version with moveout correction to vertical incidence time

% call a different version of raysum with hardwired eta changed (0.4-0.9
% crustal range cf. Godfrey et al. 2000? instead of original 1.03 for mantle according to Farra et
% al. 1991)
sizelabel = 12;

%clear all
format compact
format short g

!rm trace.out


%!cp mohoonly.mod sample.mod
%!cp moho_with_slab.mod sample.mod
%!cp amplcompslow_plunge.mod sample.mod
%!cp ampl_p_s_test.mod sample.mod
%!cp moho_Vp.mod sample.mod
%!cp moho_Vs.mod sample.mod
%!cp intracrustal.mod sample.mod
%!cp epsl_intraflat.mod sample.mod
%!cp epsl_intradip_moho.mod sample.mod
%!cp epsl_anihor.mod sample.mod
%!cp folcontrast_hor.mod sample.mod
%!cp folcontrast2int_dip.mod sample.mod
%!cp test.mod sample.mod


% ray parameter
%!cp ray04.geom sample.geom
%!cp ray08.gls -lteom sample.geom
%!cp ray06.geom sample.geom
slow = 0.06; % synth slowness in s/km for correcting moveout - should match geom file slowness

! ./raysum_grease.cmd
%! ./raysum_tatham.cmd
%! ./raysum_tensor.cmd
%! ./raysum_tathamsingletens.cmd
%! ./raysum_greasesingletens.cmd
%! ./raysum_tathamsingletenshor.cmd
%! ./raysum_greasesingletenshor.cmd
filename='trace.out';
fid=fopen(filename,'r');


% Read header:
line=fgetl(fid);
while line(1) == '#'
  line=fgetl(fid);
end
dum=str2num(line);
ntr=dum(1); nsamp=dum(2); dt=dum(3); align=dum(4); shift=dum(5);
[ntr,nsamp,dt,align,shift]

% Read each trace
line=fgetl(fid);
for itr=1:ntr
  itr;
  while line(1) == '#'
    line=fgetl(fid);
  end
  for isamp=1:nsamp
    traces(:,isamp,itr)=str2num(line);
    line=fgetl(fid);
  end
end

tvec=-shift:dt:(dt*(nsamp-1)-shift);
size(traces)

fclose(fid)

display('make sure correction slowness in this script matches slowness in geom file!')

tmax =6;

allmax = max(max(max(abs(traces))));

baz = 0:15:345;
baz = baz';

rseis = squeeze(traces(1,:,:));
ntrar = size(rseis,2);
tseis = squeeze(traces(2,:,:));
tmaxamp = max(max(abs(tseis)));
rmaxamp = max(max(abs(rseis)));
rmaxampconv = max(max(abs(rseis(6:201,:))));
tmaxampconv = max(max(abs(tseis(6:201,:))));
display('max abs conversion radial, max tangential, skipping first half second:')
[rmaxampconv tmaxampconv]
ntrat = size(tseis,2);

% moveout correction as for data
dzi = .5         ;% depth increment
Rv = 1.73               ;% Vp/Vs value
dzmax =100.0   ;% max depth
dZ = 0.0:dzi:dzmax;    ;% depth vector
ndz = size(dZ,2); 
zthk = ones(1,ndz)*dzi;
pvel = ones(ndz,1)*6.4; 
svel = pvel./Rv;
sv2  = (svel).^(-2); 
pv2 = (svel*Rv).^(-2); 
vtt = cumsum( (sqrt(sv2) - sqrt(pv2))*dzi );
nrseis = zeros(ndz,ntrar);
ntseis = zeros(ndz,ntrat);
p2 = ones(ndz,1).*slow*slow; 
mtt = cumsum( (sqrt(sv2 - p2) - sqrt(pv2-p2))*dzi );
ntt = round(real(mtt)/dt); 
ntt(1) = 1;
nrseis = rseis(ntt',:);
clear rseis;
rseis = nrseis;
ntseis = tseis(ntt',:);
clear tseis;
tseis = ntseis;
tt = vtt';

vp = 6.3;
rv = 1.73;
vs = vp/rv;
zvec = tt./(1/vs - 1/vp);

ampscale = 70;

% corrected radial and tangential:
avgrseis = mean(rseis,2);
size(rseis)
size(avgrseis)
% subtract median R RF from all R RFs
rseiscorr = rseis;
for itra = 1:ntrar
 rseiscorr(:,itra) = rseis(:,itra) - avgrseis;
end

% shift corrected T by +90 deg to match R - nodes are in strike orientation
clear tbaz
tbaz = baz+90;
clear ibazneg
ibazneg = find(tbaz>360);
tbaz(ibazneg) = tbaz(ibazneg)-360;


figure

subplot(141)
%subplot(243)
for iaz = 1:length(baz),
 fill([0 tt max(tt)],[baz(iaz) ampscale.*rseis(:,iaz)'+baz(iaz) baz(iaz)],'r')
 hold on
 ipos = find(rseis(:,iaz)>0);
 negseis =rseis;
 negseis(ipos,iaz) = 0;
 fill([0 tt max(tt)],[baz(iaz) ampscale.*negseis(:,iaz)'+baz(iaz) baz(iaz)],'b')
end
ypos = 345;
grid on
axis([0 tmax 0-baz(2) ypos+baz(2)+10])
ylabel('back-azimuth (deg)','fontsize',sizelabel,'FontName','Times')
xlabel('time(s)','fontsize',sizelabel,'FontName','Times')
h = text(tmax-4,ypos+10,'RAD','fontweight','bold','FontName','Times')
set(gca,'ytick',[0 90 180 270 360],'fontsize',sizelabel,'FontName','Times')
set(gca,'xtick',0:1:tmax,'fontsize',sizelabel,'FontName','Times')
%h=title('synthetic receiver functions ')
%set(h,'fontweight','bold','fontsize',10,'fontname','times')
%plot([0.6 0.6],[-15 385],col1)
%plot([time2 time2],[-15 385],col2)

clear negseis;


subplot(142)
%subplot(244)
for iaz = 1:length(baz),
 fill([0 tt max(tt)],[baz(iaz) ampscale.*tseis(:,iaz)'+baz(iaz) baz(iaz)],'r')
 hold on
 ipos = find(tseis(:,iaz)>0);
 negseis =tseis;
 negseis(ipos,iaz) = 0;
 fill([0 tt max(tt)],[baz(iaz) ampscale.*negseis(:,iaz)'+baz(iaz) baz(iaz)],'b')
end
ypos = 345;
grid on
axis([0 tmax 0-baz(2) ypos+baz(2)+10])
%ylabel('back-azimuth (deg)','fontsize',sizelabel')
xlabel('time(s)','fontsize',sizelabel,'FontName','Times')
set(gca,'ytick',[0 90 180 270 360],'fontsize',sizelabel,'FontName','Times')
%set(gca,'ytick',[0 90 180 270 360])
set(gca,'xtick',0:1:tmax,'fontsize',sizelabel,'FontName','Times')
%plot([0.6 0.6],[-15 385],col1)
%plot([time2 time2],[-15 385],col2)
text(tmax-3.5,ypos+10,'TAN','fontweight','bold','FontName','Times')

