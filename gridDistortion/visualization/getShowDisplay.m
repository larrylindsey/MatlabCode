function [showDisplay f] = getShowDisplay

showDisplay = false;
ch = sort(get(0, 'Children'));
for ii = 1:numel(ch)
    ud = get(ch(ii), 'UserData');
    if and(~isempty(ch), strcmp(ud, 'showExtractDistortionFigs'));
        showDisplay = true;
        f = ch(ii);
    end
end
