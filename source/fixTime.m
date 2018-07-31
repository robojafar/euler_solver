%% CALCULATE TIME REMAINING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function time = fixTime(s)
% Input: time in seconds

min = fix(s/60);
sec = rem(s,60);
sec = fix(sec);

time = sprintf('%d:%02d',min,sec);