%% Skylar Tamke - Final Project
% Will parse the input wave and find the formants after a filtering process
% to smooth the wave.

function [dataset, datachar, Storedbits, liftering_vowels] = parsewaveform(filename, phnfile)
%file that says 'Artificial intelligence is for real.'

% clc 
% clear 

% filename = 'SX295.WAV';
% phnfile = "SX295.PHN";
phntable = import_phn_table(phnfile)
% phntable = readtable("SX29phoneam.txt");

% filename2 = 'SX29.WAV';
% phnfile2 = 'SX29.PHN';
% phntable = import_phn_table(phnfile2)



%code provided on handout by Snider
fid = fopen(filename,'r');
status = fseek(fid, 1024, -1);
[wave,count] = fread(fid,inf,'int16');
fclose(fid);

Fs = 16000;


windowSize = 320;
%to change the wave length into a nice number to divide into
count = 320*(floor(count/windowSize));


numWindows = floor(count/windowSize);

windowOverlap = floor(windowSize/2);

%as of Ross's recommendation this only needs to be around 16
coefCount = 16;                                             %how many coef per section, will probably go up for end result

%length of the fft/ifft output
L = 1024;
kepstrumFilter = round(.0135 * L);
f = Fs/(L)*(1:L);

Clength = windowSize;
tempcoefCount = coefCount;

wloadindows = zeros(numWindows,windowSize);
h = zeros(numWindows,(Clength-coefCount),coefCount);
b = [];                                                                     %if the number of coefficients change this does too, so b and h have equal m vs n lengths
atotal_ls = zeros(numWindows,coefCount);


% Create a vector of overlapping windows
% This part is important since if this is done improperly the output will
% sound choppy.  If the windows are not overlapped the output will only
% show the changes between the phoneams, which is where the choppyness
% comes from.  When the windows are overlapped and combined correctly later
% the windows will blend the changes together keeping some of the
% choppyness out.

phnnumrows = height(phntable);
phnwindows = {phnnumrows};

tempphn = squeeze(phntable{:,1:2});
tempphn(1,1) =  1;  % will error out otherwise

for i = 1:phnnumrows
    phnwindows{i,:} = wave(tempphn(i,1):tempphn(i,2));
end


% for i = 1:numWindows
%     if i-1 == 0
%         windows(i,:) = wave(1:windowSize);
%     else
%         windows(i,:) = wave(((i-1)*windowOverlap+1):(((i-1)*windowOverlap)+windowSize));
%     endHelp Greg to pick a peck of potatoes
% end
% 

numWindows = phnnumrows;

fftWindows = {numWindows}
for i = 1:numWindows
   fftWindows{i,:} = fft(phnwindows{i,:},L);
end


kepstrumWindows = {numWindows};

for i = 1:numWindows
    kepstrumWindows{i,:} = ifft(log(fftWindows{i,:}),L);
end

% liftering the parts we want
liftering_vowels = kepstrumWindows;
for i = 1:numWindows
    temp = liftering_vowels{i,:};
    Storedbits(i,:) = temp(kepstrumFilter+1:end-kepstrumFilter);
    temp(kepstrumFilter+1:end-kepstrumFilter) = 0;
    
    liftering_vowels{i,:} = temp;
end

% Taking the exponential to get back into the frequency domain
delog = {numWindows};
for i = 1:numWindows
   delog{i,:} = exp((fft(liftering_vowels{i,:},L)));
end



% subplot(2,1,1)
% plot(delog)

for i = 1:numWindows
   synth{i,:} = (ifft(delog{i,:},L,"symmetric")); 
end

concat = []
for i = 1:numWindows
   concat = vertcat(concat, (synth{i,:}));
end

Storedbits;
liftering_vowels;

% subplot(2,1,2)
% plot(concat)

% soundsc(concat,16000)
% soundsc(wave,16000)
for i = 1:numWindows
    [pks,indc] = findpeaks(abs(delog{i,:}),f);

    formant1{i,:} = round(indc(1));
    formant2{i,:} = round(indc(2));
end

% %% plot the new data set from timit file
% figure(1)
% for i = 1:numWindows
% %    scatter(formant1{i,:},formant2{i,:})
%    hold on
% end
% 
% 
%% Group new dataset and character array
dataset = [];
datachar = phntable(:,3);
for i = 1:numWindows
    dataset(i,:) = [formant1{i,:} formant2{i,:}];
    
end

end
