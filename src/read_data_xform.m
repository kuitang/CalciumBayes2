% load('../data/EE210_for_Conrad.mat');

n = zeros(12,15000);
nChoices = [21 82 62 81 76 9 34 17 78 15 12 84];

for i = 1:12
    for t = 1:1201
        if(spikeInd(t) == nChoices(i))
            n(i,ceil(spikeT(t)/10)) = ...
                n(i,ceil(spikeT(t)/10)) + 1;
        end
    end
end

save('../data/real_data.mat');

