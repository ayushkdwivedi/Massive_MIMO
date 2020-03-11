function [rate_MRC,bound_MRC,rate_ZF,bound_ZF,rate_MMSE,bound_MMSE] = IPCSI(Eu,scale)
    ITER = 50;
    K = 10;
    Mv = 20:30:500;
    nTsym = K;
    rate_MRC = zeros(1,length(Mv));
    bound_MRC = zeros(1,length(Mv));
    rate_ZF = zeros(1,length(Mv));
    bound_ZF = zeros(1,length(Mv));
    rate_MMSE = zeros(1,length(Mv));
    bound_MMSE = zeros(1,length(Mv));
    for it=1:ITER

        it
        D = Dmatrix(K);
        ubeta = diag(D);

        for mx = 1:length(Mv)
            M = Mv(mx);
            if scale==1
                pu = Eu/M;
            else
                pu = Eu/sqrt(M);
            end
            Pp = nTsym*pu;
            Phi = sqrt(1/nTsym)*dftmtx(nTsym);
            pilotMtx = sqrt(Pp)*Phi;
            H = sqrt(1/2)*(randn(Mv(mx),K) + j*randn(Mv(mx),K));
            G = H*sqrt(D);
            N = sqrt(1/2)*(randn(M,nTsym)+j*randn(M,nTsym));
            RxBlk = G*(pilotMtx.')+N;
            Ghat = sqrt(1/Pp)*RxBlk*conj(Phi)*inv(1/Pp*inv(D)+eye(K));
            invGG = inv(Ghat'*Ghat);        

            for k = 1:K
                gkHat = Ghat(:,k);
                MUlmrc = 0;
                boundMUlmrc = 0;
                errZF = 0;
                errMRC = 0;
                errMMSE = 0;
                DrMRC = norm(gkHat)^2;

                for iu=1:K

                    DrErr = (nTsym*pu*D(iu,iu))+1;
                    DrMRC = DrMRC + (ubeta(iu)/DrErr)*pu*norm(gkHat)^2;
                    errZF = errZF + pu*D(iu,iu)/DrErr;
                    errMMSE = errMMSE + D(iu,iu)/DrErr;

                    if(iu~=k)
                        DrMRC = DrMRC + pu*abs(gkHat'*Ghat(:,iu))^2;
                    end
                end

                NrMRC = pu*norm(gkHat)^4;
                rate_MRC(mx) = rate_MRC(mx)+log2(1+NrMRC/DrMRC);
                Nrbound_MRC = nTsym*pu*(M-1)*ubeta(k)^2;
                Drbound_MRC = (nTsym*pu*ubeta(k)+1)*(sum(ubeta)-ubeta(k))+(nTsym+1)*ubeta(k)+(1/pu);
                bound_MRC(mx) = bound_MRC(mx)+log2(1+Nrbound_MRC/Drbound_MRC);

                Nr_ZF = pu;
                Dr_ZF = (errZF+1)*invGG(k,k);
                rate_ZF(mx) = rate_ZF(mx) + log2(1+Nr_ZF/Dr_ZF);

                Nrbound_ZF = nTsym*pu*pu*ubeta(k)^2*(M-K);
                Drbound_ZF = (nTsym*pu*ubeta(k)+1)*errZF+(nTsym*pu*ubeta(k))+1;            
                bound_ZF(mx) = bound_ZF(mx)+log2(1+Nrbound_ZF/Drbound_ZF);

                Nr_MMSE_temp = eye(K) + (errMMSE+pu^-1)^-1*Ghat'*Ghat;
                Nr_MMSE = Nr_MMSE_temp(k,k);
                rate_MMSE(mx) = rate_MMSE(mx) + log2(Nr_MMSE);

            end
        end
    end

    rate_MRC = rate_MRC/ITER;
    bound_MRC = bound_MRC/ITER;
    rate_ZF = rate_ZF/ITER;
    bound_ZF = bound_ZF/ITER;
    rate_MMSE = rate_MMSE/ITER;
    bound_MMSE = bound_MMSE/ITER;