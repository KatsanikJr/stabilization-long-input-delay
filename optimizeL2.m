function [M_bar] = optimizeL2(M1,M2,M3)
    cvx_begin sdp quiet
        variable M_bar(size(M1,1), size(M1,2))  % same dimension as each A_i
        variable R nonnegative

        minimize( R )
        subject to
            norm( M_bar - M1, 2 ) <= R;
            norm( M_bar - M2, 2 ) <= R;
            norm( M_bar - M3, 2 ) <= R;
    cvx_end
   
end