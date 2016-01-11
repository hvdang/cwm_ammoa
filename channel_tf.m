% channel_tf.m
% To generate the cross-channel parameters
% Written by Hieu V. Dang

function [Hqq,Hrq] = channel_tf(Q,Nc)

for i=1:Q
    for j=1:Nc
        Hqq(i,j) = 0.5*exp(-0.005*i*j);
        Hrq(i,j) = 0.1*exp(-0.001*i*j);
    end
end

