function [rate_MRC,bound_MRC,rate_ZF,bound_ZF,rate_MMSE,bound_MMSE] = PCSI(Eu,scale)

    ITER = 200;
    K = 10; %No. of users
    Mv = 20:30:500; %No. of BS antennas
    rate_MRC = zeros(1,length(Mv));
    bound_MRC = zeros(1,length(Mv));
    rate_ZF = zeros(1,length(Mv));
    bound_ZF = zeros(1,length(Mv));
    rate_MMSE = zeros(1,length(Mv));
    bound_MMSE = zeros(1,length(Mv));

    for it=1:ITER
    
    it
    D = Dmatrix(K);
    beta = diag(D);
    
    for mx = 1:length(Mv)
        M = Mv(mx);
        if scale==1
            pu = Eu/M;
        else
            pu = Eu/sqrt(M);
        end
        H = sqrt(1/2)*(randn(M,K) + j*randn(M,K));
        G = H*sqrt(D);
        
        for k = 1:K
            gk = G(:,k);
            nr_MRC = pu*norm(gk)^4;
            dr_MRC = norm(gk)^2;
            nr_bound_MRC = pu*(M-1)*beta(k);
            dr_bound_MRC = 1;
            
            for iu=1:K
                
                if(iu~=k)
                    dr_MRC = dr_MRC + pu*abs(gk'*G(:,iu))^2;
                    dr_bound_MRC = dr_bound_MRC + pu*beta(iu);
                end
            end
            
            rate_MRC(mx) = rate_MRC(mx)+log2(1+nr_MRC/dr_MRC);
            bound_MRC(mx) = bound_MRC(mx)+log2(1+nr_bound_MRC/dr_bound_MRC);
            
            nr_ZF = pu;
            invGG = inv(G'*G);
            dr_ZF = invGG(k,k);
            rate_ZF(mx) = rate_ZF(mx) + log2(1+nr_ZF/dr_ZF);
            
            nr_ZF_bound = pu*beta(k)*(M-K);
            bound_ZF(mx) = bound_ZF(mx)+log2(1+nr_ZF_bound);
            
            nr_MMSE = 1;
            dr_MMSE_temp = inv(eye(K) + pu*G'*G);
            dr_MMSE = dr_MMSE_temp(k,k);
            rate_MMSE(mx) = rate_MMSE(mx) + log2(nr_MMSE/dr_MMSE); 
        end
    end
    end

    rate_MRC = rate_MRC/ITER;
    bound_MRC = bound_MRC/ITER;
    rate_ZF = rate_ZF/ITER;
    bound_ZF = bound_ZF/ITER;
    rate_MMSE = rate_MMSE/ITER;
    bound_MMSE = bound_MMSE/ITER;
