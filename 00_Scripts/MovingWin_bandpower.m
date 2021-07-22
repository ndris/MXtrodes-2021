function output=MovingWin_bandpower(x,fs,band,winLen,winDisp)

% winLen and winDisp are in ms
NumWindows=floor(((length(x)/(fs/1000))-((winLen-winDisp)*1000))/(winDisp*1000));
output=zeros(1,NumWindows); % initialize output vector for number of windows 
for i=1:NumWindows
    output(i)= bandpower(x(floor(1+((i-1)*winDisp*fs)):floor((winLen*fs)+((i-1)*winDisp*fs))),fs,band);  
end
end

