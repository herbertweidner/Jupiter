%Jupiter, Carson-BW
%Autor = Herbert Weidner, 31-Jul-2026
%Kernrotation = 9 h 55 min 29,7 s = 9.925 h
%k=9.925*3600; fGW0=2/k; %5.597537083683180e-05 astro
%Jeweils 25..50 Barometer (randperm)
%-----------Zeitbereich 2000..2020 ----------------
clear %[y, Ts]=hol_DWD(3); %1=2000..2009; 2=2010..2019
k=fopen('q2000.txt'); J2000 = textscan(k,'%s %s %s %s'); fclose(k); 
nJ2000=numel(J2000{1,1}); %Anzahl der Barometer in Jahren 2000..2009
k=fopen('p2010.txt'); J2010 = textscan(k,'%s %s %s %s'); fclose(k); 
nJ2010=numel(J2010{1,1}); %Anzahl der Barometer in Jahren 2010..2019

%gobale Konstanten und Dateien definieren
Ts=3600; k=87660*2; %=length(y) im Vorgriff
tJ=365.25636042; fJ=1/tJ/24/3600; [b,a] = cheby1(6,0.01,0.005);
fZF=1e-6; %willkürlicher Wert, danach weiter reduzieren
kdec=64; %decimations-Faktor
fEnd=1/(20*kdec*Ts); %willkürlicher, glatter Wert
t=(1:Ts:21/fJ)'; tt=(1:(k))'; dftt=1e-24*t(tt); dftt2=1e-30*t(tt).^2;
    t=2*pi*t(tt); %21 Jahre FM darstellen
    t2=(1:k/kdec+1)'*2*pi*(kdec*Ts); %Zeitskala für 2.Stufe wegen decimate(y,k3)
    F1=exp(1j*fZF*t); %zum hochmixen
    F2=exp(-1j*t2*fZF); %runtermixen
    F3=exp(1j*t2*fEnd); %raufmixen
    LP=window_sinc_filter(1.5e5,1e-9,1,1/(kdec*Ts),'low','blackman');
fGW=55.97005e-6; % optimiert
drift=7.7e5; drift2=-780;
axJ=2.458; pxJ=2.698; fxJ=29.0243e-9; %29 nHz???????
%ax=0.134; px=3.6; fx=694.97e-9; %Kallisto, T=16,689 Tage
%ax2=3.035; px2=4.75; fx2=1.610062e-6; %Ganymed, T=7,155 Tage (1.627)
%ax3=1.80; px3=-0.18; fx3=3.25815e-6; %Europa, T=3,551 Tage=5114 Minuten (3.254)
%ax4=6.942; px4=-0.425; fx4=6.54522e-6; %Io, T=1,769 Tage=2547.6 Minuten (6.508)

Ph=axJ*sin(t*fxJ+pxJ);
%Ph=Ph+ax*sin(px+t*fx)+ax2*sin(px2+t*fx2)+ax3*sin(px3+t*fx3)+ax4*sin(px4+t*fx4);
fGWdrift=fGW+drift*dftt+drift2*dftt2; %gemeinsam in allen Funktionen
    %ys=exp(-1j*(t.*fGWdrift+Ph)); %lokaler Oszi mit Jahres-PM
    ys=sin(t.*fGWdrift+Ph); %lokaler Oszi mit Jahres-PM (Dauer=20 Jahre)
    x=T_zero(ys,Ts); %Zeitpunkte der Nulldurchgänge
    erg=diff(x); %Abstände der Nulldurchgänge
    %2*8929 s -> 56 µHz
    [b,a] = cheby1(6,0.02,0.03);
    ys=flipud(filtfilt(b,a,flipud(erg-mean(erg))))+mean(erg); 
    yy = spline(x(1:end-1),ys,t/2/pi);
    plot(1./(2*yy)), ylabel('redshift ---- blueshift Hz')
    %bis hierher gilt Ts=3600
    plot(t/(2*pi*3600*24*365.256),5e5./yy), xlabel('Zeit (Jahre)')
    ylabel('(red) RX-Frequenz (µHz) (blue)'), title('GW from Jupiter')
return

R(1,1:20)=mean(Res8(1:661,1:20)); %Mittelwerte
R(2,1:20)=std(Res8(1:661,1:20)); %Fehler

%PSD berechnen
h6=window_sinc_filter(2e3,fZF-20e-9,1,1/Ts,'low','blackman');
yi=t.*F0+z2; ys=y.*sin(yi); yc=y.*cos(yi); 
y1=conv(ys,h6,'same').*stt+conv(yc,h6,'same').*ctt;
yi=conv(ys,h6,'same').*stt-conv(yc,h6,'same').*ctt; %invers
[sp,~]=zeig_sp2(y1,Ts,3000.4); [spi,f]=zeig_sp2(yi,Ts,3000.4);
plot(1e9*f,sp,1e9*f,spi)
k=16*1024;[sp,f]= pwelch(cat(1,y1,zeros(3e5,1)),rectwin(k),k/8,16*k,1/Ts);
plot(1e9*f,sp/1e5), xlabel('Frequency (nHz)')
return

%a und b sind die Frequenzgrenzen für Rekonstr. von y
%peak ist das bin, dessen peak gemessen wird
function [yss]=FFT_LP(y,a) %a<numel(y)/2 entspricht max Frequenz
%a und b begrenzen den Bereich des Maximums
%L=(1:numel(y))'; %nur, falls f-Skala benötigt wird:
Y = fft(y); %plot(abs(Y)) 
    %Ts=3600; NFFT = 2^(nextpow2(length(y))+0);
    %f=1/Ts/2*linspace(0,1,NFFT/2+1); plot(1e6*f,abs(Y(1:NFFT/2+1))), xlabel('Frequency (µHz)')
Y(a:numel(y)+2-a)=0; %Tiefpass, Summe=8194!
%Y(b:numel(y)-b)=0; %Tiefpass, Summe=8192+2! 
%Y(1:a)=0; Y(L-a:end)=0; %Hochpass
yss=ifft(Y); 
%yss=yss(L);
end


%----------------------------------------
%
%A,F,P sind Kenndaten für 29 nHz.
function Richtung(A,F,P,t,fGWdrift,Ph,Ts)
%----------y ist IQ -------------
    %ys=exp(-1j*(t.*fGWdrift+Ph)); %lokaler Oszi mit Jahres-PM
    ys=sin(t.*fGWdrift+Ph); %lokaler Oszi mit Jahres-PM (Dauer=20 Jahre)
    x=T_zero(ys,Ts); %Zeitpunkte der Nulldurchgänge
    erg=diff(x); %Abstände der Nulldurchgänge
    %2*8929 s -> 56 µHz
    [b,a] = cheby1(6,0.02,0.03);
    ys=flipud(filtfilt(b,a,flipud(erg-mean(erg))))+mean(erg); 
    yy = spline(x(1:end-1),ys,t/2/pi);
    plot(1./(2*yy)), ylabel('redshift ---- blueshift Hz')
    %bis hierher gilt Ts=3600
    plot(t/(2*pi*3600*24*365.256),5e5./yy), xlabel('Zeit (Jahre)')
    ylabel('RX-Frequenz (µHz)'), title('GW from Jupiter')
end

function erg=T_zero(y,Ts) %ohne Zeitzuordnung
yc=y.*cat(1,0,y(1:end-1)); erg=zeros(1e5,1); %1=Zeitpunkt der Nulldurchgänge
ys=(yc<0); %Vorzeichenwechsel
nn=1; j=1; while ys(j)==0, j=j+1; end %1. VZ-Wechsel
b=j-1+y(j-1)/(y(j-1)-y(j)); %genauer 0-Durchgang 
while j<length(ys)-20
    j=j+1; while ys(j)==0, j=j+1; end %nächster Wechsel
    a=j-1+y(j-1)/(y(j-1)-y(j)); %genauer 0-Durchgang 
    erg(nn,1)=Ts*(a+b)/2; nn=nn+1; b=a;
end 
erg=erg(1:nn-1,:);
end

function erg=T_zaehl_T(y,Ts) %mit Zeitzuordnung
yc=y.*cat(1,0,y(1:end-1)); erg=zeros(1e5,2); %1=Zeit, 2=Zeitdiff
ys=(yc<0); %Vorzeichenwechsel
nn=1; j=1; while ys(j)==0, j=j+1; end %1. VZ-Wechsel
b=j-1+y(j-1)/(y(j-1)-y(j)); %genauer 0-Durchgang 
while j<length(ys)-30
    j=j+1; while ys(j)==0, j=j+1; end %nächster Wechsel
    a=j-1+y(j-1)/(y(j-1)-y(j)); %genauer 0-Durchgang 
    erg(nn,1)=Ts*(a+b)/2;
    erg(nn,2)=Ts*(a-b); nn=nn+1; b=a;
end 
erg=erg(1:nn-1,:);
end


function y=sync(yGW, f, Ts)
a=2*pi*Ts*f; b=sin(a); a=cos(a); L=length(yGW);
y=zeros(L,1); y1=y;
if yGW(1)>0, y(1)=1; %willkürlicher Startwert
else, y(1)=-1;
end
for j=2:L %vorwärts
y(j)=y(j-1)*a+y1(j-1)*b+yGW(j); %Add-Theorem
y1(j)=y1(j-1)*a-y(j-1)*b;
end%, figure(1),plot(y), title('Coherent Demodulation')
end

%frequenztransformierendes Bandfilter um f_xx
%mit wählbarer Flankensteilheit (n)
function ssb = filtTrans_n_Hz(y,Ts,n,f_ein,f_aus,B)
%f_xx und B(andbreite) in Hz
%n ist die Länge von BF
%wenn B<0 wird cheby-Filter verwendet
L=length(y); j=(1:L)'*2*pi*f_ein*Ts;
ys=y.*sin(j); yc=y.*cos(j);
if B>0
    BF=window_sinc_filter(n,B,1,1/Ts,'low','blackman');
    ys=conv(ys,BF,'same'); yc=conv(yc,BF,'same');
else
    [b,a] = cheby1(6,0.01,-B); %0.006
    yc=filtfilt(b,a,yc); ys=filtfilt(b,a,ys);
end
if f_ein~=f_aus
    j=(1:L)'*2*pi*f_aus*Ts; %verschobene MittenFrequenz
end
ys=ys.*sin(j); yc=yc.*cos(j);
ssb=2*(yc+ys).*wei(100,L); %Ts ist unverändert, B wie gewählt
%[erg f]=zeig_sp2(ssb,-Ts,'log');
end
%

%y=Signal (variable Frequenz), gesampled mit Ts
%stt und ctt, gesampled mit 1/fZF erzeugen aus y die ZF (feste Frequenz)
%

function [A,F,P,Ph,yi]=ITER(A,F,P,t,fGWdrift,y,Ph,Ts,b,a,LP,F1,F2,F3)
if A==0, yi=0; return, end 
erg=zeros(2000,2); dd=0.001;
Ph=Ph-A*sin(P+t*F); %nach 3 Iterationen wieder einsetzen
for k=1:size(erg,1), z=A*sin(P+t*F)+Ph; 
    ys=y.*exp(-1j*(t.*fGWdrift+z)); %runtermixen mit variabler Frequenz
    ys=filtfilt(b,a,ys).*F1; %filtern, hochmixen auf 1 µHz
    %bis hierher gilt Ts=3600, fZF=1000 nHz
    yi=decimate(decimate(ys,8,'fir'),8,'fir'); %Ts=64*Ts;
    %ab jetzt: Ts2=230400, Signalfreq= fEnd~193 nHz, schmalbandig
    ys=yi.*F2; %runtermixen
    yi=conv(ys,LP,'same').*F3; %filtern, raufmixen auf fEnd~193 nHz
    erg(k,:)=[P,sum(abs(yi))]; %max suchen
    if k==2 && (erg(k,2))<(erg(k-1,2))
        dd=-dd; P=P+dd; %Suchrichtung umkehren
        erg(5,:)=erg(1,:); erg(1,:)=erg(2,:); erg(2,:)=erg(5,:);        
    else, if k>1 && erg(k,2)<erg(k-1,2), break, end
    end, P=P+dd;
end, p=polyfit(erg(k-2:k,1)-P,erg(k-2:k,2),2); %max suchen
P=P-p(2)/2/p(1);

df=F*1e-6; 
for k=1:size(erg,1), z=A*sin(P+t*F)+Ph; 
    ys=y.*exp(-1j*(t.*fGWdrift+z)); %runtermixen mit variabler Frequenz
    ys=filtfilt(b,a,ys).*F1; %filtern, hochmixen auf 1 µHz
    %bis hierher gilt Ts=3600, fZF=1000 nHz
    yi=decimate(decimate(ys,8,'fir'),8,'fir'); %Ts=64*Ts;
    %ab jetzt: Ts2=230400, Signalfreq= fEnd~193 nHz, schmalbandig
    ys=yi.*F2; %runtermixen
    yi=conv(ys,LP,'same').*F3; %filtern, raufmixen auf fEnd~193 nHz
    erg(k,:)=[F,sum(abs(yi))]; %max suchen
    if k==2 && (erg(k,2))<(erg(k-1,2))
        df=-df; F=F+df; %Suchrichtung umkehren
        erg(5,:)=erg(1,:); erg(1,:)=erg(2,:); erg(2,:)=erg(5,:);        
    else, if k>1 && erg(k,2)<erg(k-1,2), break, end
    end, F=F+df;
end, p=polyfit(1e9*(erg(k-2:k,1)-F),erg(k-2:k,2),2); %max suchen
F=F-p(2)/2e9/p(1);

for k=1:size(erg,1), z=A*sin(P+t*F)+Ph; 
    ys=y.*exp(-1j*(t.*fGWdrift+z)); %runtermixen mit variabler Frequenz
    ys=filtfilt(b,a,ys).*F1; %filtern, hochmixen auf 1 µHz
    %bis hierher gilt Ts=3600, fZF=1000 nHz
    yi=decimate(decimate(ys,8,'fir'),8,'fir'); %Ts=64*Ts;
    %ab jetzt: Ts2=230400, Signalfreq= fEnd~193 nHz, schmalbandig
    ys=yi.*F2; %runtermixen
    yi=conv(ys,LP,'same').*F3; %filtern, raufmixen auf fEnd~193 nHz
    erg(k,:)=[A,sum(abs(yi))]; %max suchen
    if k==2 && (erg(k,2))<(erg(k-1,2))
        dd=-dd; A=A+dd; %Suchrichtung umkehren
        erg(5,:)=erg(1,:); erg(1,:)=erg(2,:); erg(2,:)=erg(5,:);        
    else, if k>2 && erg(k,2)<erg(k-1,2), break, end
    end, A=A+dd;
end, p=polyfit(erg(k-2:k,1)-A,erg(k-2:k,2),2); %max suchen
A=A-p(2)/2/p(1); 
Ph=Ph+A*sin(P+t*F); %nach 3 Iterationen wieder einsetzen
end

function glocke=wei(k,L) %k~2000
%Flanke des Weidnerfensters wählen
glocke=ones(L,1); 
if k>0, j=hann(2*k); else k=-k; j=hamming(2*k); end
glocke(1:k)=glocke(1:k).*j(1:k); %Enden ignorieren
glocke(L-k:L)=glocke(L-k:L).*j(k:2*k);

%Spektrum der Funktion y zeigen
function [erg, f]=zeig_sp2(y,Ts,w) %Ts in s
%w ist die Anzahl der angezeigten Punkte
%die angefügte Dezimale wird zur Erhöhung
%der Auflösung zu NFFT addiert.
%Die effektive Bildbreite bleibt konstant
    %figure(1)
%[b,a] = cheby1(4,0.1,0.95); y=filtfilt(b,a,y);
[z,p,k] = cheby1(6,0.1,0.95); [sos,g] = zp2sos(z,p,k);
y=g(1,1)*sosfilt(sos,y);
%if Ts<0, Ts=-Ts; y=y.*blackman(length(y));
if Ts<0, Ts=-Ts; y=y.*hamming(length(y));
end
if isnumeric(w)
  if w<0, w=-w; a=1; else, a=0; end %kein Graph
    verg=round(10*(w-floor(w)));
    NFFT = 2^(nextpow2(length(y))+verg);% z.B. +3
    Fs=1/Ts; %Ein Wert pro Minute?
    f = Fs/2*linspace(0,1,NFFT/2+1);
    Y = fft(y,NFFT);
    erg=abs(Y(1:NFFT/2+1));
    w=floor(w);
    w=min(2^verg*floor(w),length(f));
    erg=erg(1:w); f=f(1:w);
    if a==0
      plot(1e6*f(1:w),erg(1:w)) %Ausschnitt vergrößern
      %semilogy(1e6*f,erg) %Ausschnitt vergrößern
      title('Spektrum'), xlabel('Frequenz in µHz')
      ylabel('relative Amplitude')
    end
else
    NFFT = 2^(nextpow2(length(y))+2);%Verfeinerung
    Fs=1/Ts; %Ein Wert pro Minute?
    f = Fs/2*linspace(0,1,NFFT/2+1);
    Y = fft(y,NFFT);
    erg=abs(Y(1:NFFT/2+1));
    if strcmp(w,'lin'),plot(f,erg),end %Ausschnitt vergrößern
    if strcmp(w,'log'),semilogy(f,erg),end
    title('Spektrum'), xlabel('Frequenz in Hz')
    ylabel('relative Amplitude')
end



