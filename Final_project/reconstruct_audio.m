%% Reconstruct_audio
% This function will recreate an audio wave with the inputs., the idea is
% to be able to take different liftered sets of data and apply them to a
% different speaker
function concat = reconstruct_audio(lifteringVowels, storedBits)

lifteredData = lifteringVowels%{1,:}
fillIn = storedBits%{1,2}

Fs = 16000;
windowSize = 320;
%to change the wave length into a nice number to divide into
count = 40650
count = 320*(floor(count/windowSize));
numWindows = floor(count/windowSize);
%need to keep this the same as what is in parsewaveform.m
L = 1024;    
kepstrumFilter = 14

numWindows = 35; %same as number of phoneams

for i = 1:numWindows
    i
    temp = lifteredData{i,:};
    temp(kepstrumFilter+1:end-kepstrumFilter) = fillIn(i,:);
    
    lifteredData{i,:} = temp;
end


% Taking the exponential to get back into the frequency domain
delog = {numWindows};
for i = 1:numWindows
   delog{i,:} = exp((fft(lifteredData{i,:},L)));
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


% subplot(2,1,2)
% plot(concat)

 soundsc(concat,16000)
% soundsc(wave,16000)



end