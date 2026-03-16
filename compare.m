% load("D:\과제\삼성 과제\Data\2차년도\Processed\BMHR21030714D_2025-09-19_15-25-44_processed.mat")
load("D:\과제\삼성 과제\Data\2차년도\Processed\receiver_opensky_0514_processed.mat")


idx = 108;

true_val = [-1.853262e+07  3.089512e+07 -2.199198e+07 -157.046868  1279.862607  1902.538736  85472.79625];
current_val = [SVpos_x(1, idx) SVpos_y(1, idx) SVpos_z(1,idx) SVvel_x(1, idx), SVvel_y(1, idx), SVvel_z(1, idx), sv_clock_bias(1,idx)];

current_val - true_val

for i=1:5
    if idx < constellation_idx(i+1) && idx >= constellation_idx(i)
        constellation_name(i)
        idx - constellation_idx(i) + 1
        break
    end
end