%% Need to change project description a little bit.
% can do a few different things from here.
% 1. swap vocal cords with vocal tracts
% 2. average and output synthesized


clc
clear


wavfiles = struct2cell(dir("*.WAV"))
phnfiles = struct2cell(dir("*.PHN"))

numberoffiles = length(wavfiles(1,:))



dataset = {numberoffiles}
charset = {numberoffiles}
storedBits = {numberoffiles}
lifteringVowels = {numberoffiles}

% wavcells= =struct2cell(wavfiles)
% phncells= =struct2cell(phnfiles)

%% Formant finding function
% based on the code that we wrote for project 2, this finds the formants
% for each files and outputs corresponding formants and the pheonems that
% go with them.
for i = 1:numberoffiles
    [dataset{:,:,i},charset{:,i},storedBits{:,i},lifteringVowels{:,i}] = parsewaveform(wavfiles{1,i},phnfiles{1,i});
end



% could take one speaker and swap formants with another speaker to see when
% they sound like

%% reconstruct audio wave

testoutput = reconstruct_audio(lifteringVowels{:,1},storedBits{:,2})


% 
% for i = 1:numberoffiles
%    chartable{i} = unique(charset{i});
% end
% 
% testcat = []
% 
% for i = 1:numberoffiles
%    testcat = vertcat(testcat,chartable{i}); 
% end
% 
% %all characters are sorted for splitting all corresponding values into
% %indivisual vectors, 
% charlist = unique(testcat)
% charnum = height(charlist)
% 
% datasorted = {}
% 
% for i = 1:numberoffiles
%    typelength = length(dataset{:,:,i});
%    tempdata = dataset{:,:,i};
%    temptable = charset{:,i};
%    for j = 1:typelength
%        for k = 1:charnum
%             if(charlist{k,:} == temptable{j,:})
%                 sortedData{k,j,:} = tempdata(j,:);
%             end
%        end  
%    end
% end
% 
% 
% %% Averaging PHN values based on location
% empty = 0;
% nonempty = 0;
% totalVal = 0;
% for i = 1:charnum
%    for j = 1:typelength
%       if(isempty(sortedData{i,j,:}))
%           empty = empty + 1;
%       else
%           nonempty = nonempty + 1;
%           totalVal = totalVal + sortedData{i,j,:}
%       end
%    end
%    total(i,:,:) = totalVal
%    totalVal = 0;
%    values(i) = nonempty;
%    nonempty = 0;
% end
% 
% 
% for i = 1:charnum
%    Averagephn(i,1) = total(i,:,1)/values(i) 
%    Averagephn(i,2) = total(i,:,2)/values(i) 
% end
% 
% 









%% Scaling to sound
% as Ross put it in class one day: s = s/max(abs(s))

