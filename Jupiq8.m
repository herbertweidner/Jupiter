%Jupiter, Carson-BW
%Autor = Herbert Weidner, 13-Jul-2026
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
Res=zeros(1000,20); nBaro=40; %<nJ2000
for Respoi=1:400 %Jeweils nBaro Barometer zusammenfassen, 400 Versuche
    tt=zeros(2,nBaro); Ts=3600; 
    y=0; %Summe von nBaro Barometern, aneinandergefügt 2000..2010..2019
    %Generate a random permutation of nBaro unique integers (without repeating elements) selected randomly from the integers 1 to nJ2000.
    tt(1,:)= randperm(nJ2000,nBaro); %Zufallszahlen 1..nJ2000 erzeugen
    tt(2,:)= randperm(nJ2010,nBaro); %Zufallszahlen 1..nJ2010 erzeugen
    for k=1:nBaro
        zeile1=J2000{1,2}{tt(1,k)}; zeile1 = strip(zeile1, ''''); y1={load(zeile1)}; 
        zeile2=J2010{1,2}{tt(2,k)}; zeile2 = strip(zeile2, ''''); y2={load(zeile2)}; 
        y=y+cat(1,y1{1,1}.y,y2{1,1}.y); %plot(y)
    end
ys=FFT_LP(y,7400); %fGrenz ~ 11.5 µHz 
y=(y-ys).*wei(40,numel(y)); %alle tieferen Frequenzen löschen

y=hilbert(y); %komplexes Basisband
fGW=55.97005e-6; % optimiert
drift=8.3e5; drift2=-820;
axJ=2.46; pxJ=2.7; fxJ=29.025e-9; %29 nHz???????
ax=0.134; px=3.6; fx=694.97e-9; %Kallisto, T=16,689 Tage
ax2=3.035; px2=4.75; fx2=1.610062e-6; %Ganymed, T=7,155 Tage (1.627)
ax3=1.80; px3=-0.18; fx3=3.25815e-6; %Europa, T=3,551 Tage=5114 Minuten (3.254)
ax4=6.942; px4=-0.425; fx4=6.54522e-6; %Io, T=1,769 Tage=2547.6 Minuten (6.508)

R=zeros(398,21);
%fakF=-4.8e-13; fakD1=-3.4e5; fakD2=-27000;
fakF=-4.8e-13; fakD1=-3.4e5; fakD2=-27000;
for k3=1:size(R,1)
z2=axJ*sin(t*fxJ+pxJ);
z2=z2+ax*sin(px+t*fx)+ax2*sin(px2+t*fx2)+ax3*sin(px3+t*fx3)+ax4*sin(px4+t*fx4);
fGWdrift=fGW+drift*dftt+drift2*dftt2; %gemeinsam in allen Funktionen

%jetzt die Modulationen anpassen
[axJ,fxJ,pxJ,z2,yi]=ORBIT(axJ,fxJ,pxJ,t,fGWdrift,y,z2,Ts,fZF,b,a,LP,F1,F2,F3);
[ax,fx,px,z2,~]=ITER(ax,fx,px,t,fGWdrift,y,z2,Ts,b,a,LP,F1,F2,F3);
[ax2,fx2,px2,z2,~]=ITER(ax2,fx2,px2,t,fGWdrift,y,z2,Ts,b,a,LP,F1,F2,F3);
[ax3,fx3,px3,z2,~]=ITER(ax3,fx3,px3,t,fGWdrift,y,z2,Ts,b,a,LP,F1,F2,F3);
[ax4,fx4,px4,z2,yEnd]=ITER(ax4,fx4,px4,t,fGWdrift,y,z2,Ts,b,a,LP,F1,F2,F3);

z0=T_zero(real(yEnd),kdec*Ts); z00=diff(z0(:,1)); %z00=z00(3:end);
if mod(k3,30)==0, plot(z00), drawnow, end
k=mean(z00)-0.5/fEnd; %disp(k) %Abweichung von fEnd(Iter) (2.268e6)
fGW=fGW+k*fakF; %Frequenz korrigieren
p=polyfit(1:size(z00),z00,2); drift=drift+p(2)*fakD1; drift2=drift2+p(1)*fakD2;
R(k3,:)=[fGW,drift,drift2,axJ,pxJ,fxJ,ax,px,fx,ax2,px2,fx2,ax3,px3,fx3,ax4,px4,fx4,k,p(1),p(2)];
end, R=R(1:k3-1,:); %return
[sp,f]=zeig_sp2(yi,-kdec*Ts,-100000.5); 
plot(f*1e9,sp), xlabel('Frequency (nHz)')
Res(Respoi,1:20)=mean(R(end-50:end,1:20));
end, return
R(1,1:20)=mean(Res8(1:661,1:20)); %Mittelwerte
R(2,1:20)=std(Res8(1:661,1:20)); %Fehler

%PSD berechnen
k=16*1024;[sp,f]= pwelch(cat(1,yi,zeros(3e5,1)),rectwin(k),k/8,16*k,1/(64*Ts));
plot(1e9*f,sp/1e6), xlabel('Frequency (nHz)')
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
function [A,F,P,Ph,yi]=ORBIT(A,F,P,t,fGWdrift,y,Ph,Ts,fZF,b,a,LP,F1,F2,F3)
if A==0, return, end 
%----------y ist IQ -------------
erg=zeros(2000,2); dd=0.001;
Ph=Ph-A*sin(t*F+P); %nach 3 Iterationen wieder einsetzen
%zuerst A,F,P optimieren
for k=1:size(erg,1), z=Ph+A*sin(t*F+P); 
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

df=F*5e-8; %1e-7
for k=1:size(erg,1), z=Ph+A*sin(t*F+P); 
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
end, p=polyfit(1e11*(erg(k-2:k,1)-F),erg(k-2:k,2),2); %max suchen
F=F-p(2)/2e11/p(1);

for k=1:size(erg,1), z=Ph+A*sin(t*F+P); 
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
Ph=Ph+A*sin(t*F+P); %nach 3 Iterationen wieder einsetzen
end

function erg=T_zero(y,Ts) %ohne Zeitzuordnung
yc=y.*cat(1,0,y(1:end-1)); erg=zeros(3000,1); %1=Zeitpunkt der Nulldurchgänge
ys=(yc<0); %Vorzeichenwechsel
nn=1; j=1; while ys(j)==0, j=j+1; end %1. VZ-Wechsel
b=j-1+y(j-1)/(y(j-1)-y(j)); %genauer 0-Durchgang 
while j<length(ys)-200
    j=j+1; while ys(j)==0, j=j+1; end %nächster Wechsel
    a=j-1+y(j-1)/(y(j-1)-y(j)); %genauer 0-Durchgang 
    erg(nn,1)=Ts*(a+b)/2; nn=nn+1; b=a;
end 
erg=erg(1:nn-1,:);
end

function erg=T_zaehl_T(y,Ts) %mit Zeitzuordnung
yc=y.*cat(1,0,y(1:end-1)); erg=zeros(300,2); %1=Zeit, 2=Zeitdiff
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

function [output] = window_sinc_filter(M, fc1, fc2, fs, filter_type, window_type) 
% Die Funktion berechnet die Filterkoeffizienten eines Window-Sinc Filters. 
% Hierzu werden die sinc Funktion sin(x)/x und ein wählbares Fenster verwendet. 
% 
% Inputparameter: 
% order        = Filterordnung ->  Anzahl Filterkoeff. = order + 1 
%                Bedingung: order = gerade ganze Zahl 
% 
% fc1 [Hz]     = cutoff frequency 1 = Grenzfrequenz der ersten Filterstufe 
%                -> für Hochpass, Tiefpass, Bandpass und Bandsperre 
%                Bedingung: 0 < fc1 < fs/2 
% 
% fc2 [Hz]     = cutoff frequency 2 = Grenzfrequenz der zweiten Filterstufe 
%                -> neben fc1 zusätzlich nur für Bandpass und Bandsperre 
%                Bedingung: 0 < fc1 < fc2 < fs/2 
% 
% fs [Hz]      = Sampling-frequency = Abtastfrequenz 
% 
% filter_type  = 'low' = Tiefpass, 'high' = Hochpass,, 'bandpass' = Bandpass 
%                und 'bandreject' = Bandsperre 
% 
% window_type  = Fenstertyp: 'Hamming' oder 'Blackman' 
% 
% analyse_plot = Bodediagramm, Frequenzantwort linear und Sprungfunktion(nur Tiefpass) 
%                des Filters darstellen: 'y' oder 'n' 
% 
% Outputparameter: 
% output     = Filterkoeffizienten = Impulsantwort des Filters 
% 
% Quelle des Filter-Algorithmus: 
% The Scientist and Engineer's Guide to Digital Signal Processing 
% By Steven W. Smith, Ph.D. , Chapter 16: Window-Sinc Filter 
% http://www.dspguide.com/ch16.htm 
% ------------------------------------------------------------------------- 

% Überprüfung der Inputparameter: 
% Grenzfrequenz 1 überprüfen 
Fc1 = fc1 / fs ; % Normierung auf fn (Fc1 = 0...0.5 fs) 
if (Fc1 <= 0) || (Fc1 >= 0.5) 
    error(['Wrong cutoff frequency fc1 [Hz]: ', sprintf('%8.1f',fc1),' ...choose: 0 < fc1 < ',sprintf('%8.1f',fs/2),' = fs/2']); 
end 
% Grenzfrequenz 2 überprüfen 
% wird nur für Bandpass- und Bandsperrfilter benötigt 
if (strcmp(filter_type, 'bandpass')) || (strcmp(filter_type, 'bandreject')) 
    Fc2 = fc2 / fs ; % Normierung auf fn (Fc2 = 0...0.5 fs) 
    if ((Fc2 <= 0) || (Fc2 >= 0.5) || (Fc2 <= Fc1)) 
        error(['Wrong cutoff frequency fc2 [Hz]: ', sprintf('%8.1f',fc2),' ...choose: ',sprintf('%8.1f',fc1) ' < fc2 < ',sprintf('%8.1f',fs/2),' = fs/2']); 
    end 
end 

% Auswahl Filtertyp 
switch lower(filter_type) 
    case 'high' % Hochpass - Filter 
        B = zeros(M+1, 1); % Init 
        % Fensterfunktion erstellen 
        window = Fenster(M+1, window_type); 
        for i = 0:M 
            if 2 * i == M 
                B(i+1) = 2*pi*Fc1; 
            else 
                B(i+1) = sin(2*pi*Fc1 * (i-(M/2))) / (i-(M/2)); 
            end 
            B(i+1) = B(i+1) * window(i+1); 
        end                 

        % Verstärkungsfaktor des Filters auf 1 normieren 
        B = B./sum(B); 
        % Tiefpass in Hochpass durch Inversion des Spektrums wandeln 
        output      = - B; 
        output((M/2)+1) = output((M/2)+1) + 1; 
        
    case 'low' % Tiefpass - Filter 
        B = zeros(M+1, 1); % Init 
        window = Fenster(M+1, window_type); 
        for i = 0:M 
            if 2 * i == M   % Multiplikation ist schneller als Division 
                B(i+1) = 2*pi*Fc1; 
            else 
                B(i+1) = sin(2*pi*Fc1 * (i-(M/2))) / (i-(M/2)); 
            end 
            B(i+1) = B(i+1) * window(i+1); 
        end 

        % Verstärkungsfaktor des Filters auf 1 normieren 
        B      = B./sum(B); 
        output = B; 

    case 'bandpass' % Bandpass - Filter 
        A = zeros(M+1, 1); B=A; % Init 
        % Fensterfunktion erstellen 
        window = Fenster(M+1, window_type); 
        
        % 1. Filterstufe: Berechnung der Filterkoeffizienten 
        for i = 0:M 
            if 2 * i == M   % Multiplikation ist schneller als Division 
                A(i+1) = 2*pi*Fc1; 
            else 
                A(i+1) = sin(2*pi*Fc1 * (i-(M/2))) / (i-(M/2)); 
            end 
            A(i+1) = A(i+1) * window(i+1); 
        end 

        % 2. Filterstufe: Berechnung der Filterkoeffizienten 
        for i = 0:M 
            if 2 * i == M 
                B(i+1) = 2*pi*Fc2; 
            else 
                B(i+1) = sin(2*pi*Fc2 * (i-(M/2))) / (i-(M/2)); 
            end 
            B(i+1) = B(i+1) * window(i+1); 
        end 

        % Verstärkungsfaktor des Filters auf 1 normieren 
        A = A./sum(A); B = B./sum(B); 
        % Tiefpass in Hochpass durch Inversion des Spektrums wandeln 
        B          = - B; 
        B((M/2)+1) = B((M/2)+1) + 1; 
        output     = A + B; 
        % Bandsperre in Bandpass durch Inversion des Spektrums wandeln 
        output          = - output; 
        output((M/2)+1) = output((M/2)+1) + 1; 

    case 'bandreject' % Bandsperr - Filter 
        A = zeros(M+1, 1); B=A; % Init 
        % Fensterfunktion erstellen 
        window = Fenster(M+1, window_type); 

        % 1. Filterstufe: Berechnung der Filterkoeffizienten 
        for i = 0:M 
            if 2 * i == M   % Multiplikation ist schneller als Division 
                A(i+1) = 2*pi*Fc1; 
            else 
                A(i+1) = sin(2*pi*Fc1 * (i-(M/2))) / (i-(M/2)); 
            end 
            A(i+1) = A(i+1) * window(i+1); 
        end 

        % 2. Filterstufe: Berechnung der Filterkoeffizienten 
        for i = 0:M 
            if 2 * i == M   
                B(i+1) = 2*pi*Fc2; 
            else 
                B(i+1) = sin(2*pi*Fc2 * (i-(M/2))) / (i-(M/2)); 
            end 
            B(i+1) = B(i+1) * window(i+1); 
        end 

        % Verstärkungsfaktor des Filters auf 1 normieren 
        A = A./sum(A); B = B./sum(B); 
        % Tiefpass in Hochpass durch Inversion des Spektrums wandeln 
        B          = - B; 
        B((M/2)+1) = B((M/2)+1) + 1; 
        output = A + B; 
        
    otherwise % ungültiger Filtertyp 
        error(['Unknown filter type: ', filter_type ' ...choose: high, low, bandpass or bandreject']); 
end % switch 

%-----------------------------------
function output = Fenster(window_size, window_type) 
% Input der Funktion: 
% window_size = Fensterlänge 
% window_typ = 'Hamming', 'Hann', 'Blackman', 'Blackman-Harris' 
N = window_size; 
% Fensterlänge output = zeros(N, 1); 
if mod(N, 2) == 0 % N gerade 
  m = fix(N / 2); n = m; 
else % N ungerade 
  m = fix(N / 2)+1; n = m-1; 
end
switch lower(window_type) 
  case 'hamming' 
    window = 0.54 - 0.46 * cos(2*pi*(0:m) / (N-1)); 
  case 'hann' 
    window = 0.50 - 0.50 * cos(2*pi*(0:m) / (N-1)); 
  case 'blackman' 
    window = 0.42 - 0.50 * cos(2*pi*(0:m) / (N-1)) + 0.08 * cos(4*pi* (0:m) / (N-1)); 
  case 'blackmanharris' 
  window = 0.35875 - 0.48829 * cos(2*pi*(0:m) / (N-1)) + 0.14128 * cos(4*pi* (0:m) / (N-1)) - 0.01168 * cos(6*pi* (0:m) / (N-1)); 
  otherwise, error(['Unknown window type: ', window_type]); 
end % Ergebnisvektor 
output = transpose([window(1:m),window(n:-1:1)]);
end
end

function glocke=wei(k,L) %k~2000
%Flanke des Weidnerfensters wählen
glocke=ones(L,1); 
if k>0, j=hann(2*k); else k=-k; j=hamming(2*k); end
glocke(1:k)=glocke(1:k).*j(1:k); %Enden ignorieren
glocke(L-k:L)=glocke(L-k:L).*j(k:2*k);
end

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
%Yb=abs(Y(1:NFFT/2+1));
%plot(Yb(1:5000,1)) 
%semilogy(Yb(1:5000,1))
end
