function y = smscma_modu(bit_x,b1,b2,V,K,CB,H)
    
    u = bit_x';
    u1 = u(1:b1,:);
    u2 = u(b1+1:b1+b2,:);
    ante = bi2de_pzp(u1')+1; % column vector 00 --> 1 01-->2, etc.
    scma_codeword_index = bi2de_pzp(u2')+1; % column vector
    trans_sym = zeros(K,V)+1i*zeros(K,V);
    for ii = 1:V
        trans_sym(:,ii) = CB(:,scma_codeword_index(ii),ii);
    end
    
    y = zeros(K,1);

    for vv = 1:V
        y = y + H(:,vv,ante(vv)).*trans_sym(:,vv);
    end

        
end


%bi2de_pzp
function a = bi2de_pzp(b)

    [m,n] = size(b);
    a = zeros(m,1);
    for jj = 1:m
        for ii = 1:n
            a(jj) = a(jj) + b(jj,ii)*2^(n-ii);
        end
    end
    
end