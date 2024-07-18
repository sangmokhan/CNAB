% SangMok Han, Turbulence Lab, Yonsei University, August 2023
function [cmax] = converge(u,u_prev,v,v_prev)
    cmax1 = max(max(abs(u-u_prev)));
    cmax2 = max(max(abs(v-v_prev)));
    cmax = max(cmax1,cmax2);
end