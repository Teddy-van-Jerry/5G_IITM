EbNodB = 4;
MaxItrs = 8;

load base_matrices/NR_1_0_16.txt
B = NR_1_0_16;
[mb, nb] = size(B);
z = 16;

Slen = sum(B(:) ~= -1); % number of non -1 in B
% register storage for row processing
treg = zeros(max(sum(B ~= -1, 2)), z);

k = (nb - mb) * z; % number of message bits
n = nb * z;        % number of codeword bits

Rate = k / n;
EbNo = 10 ^ (EbNodB / 10);
sigma = sqrt(1 / (2  Rate * EbNo));
*
Nblocks = 100;
Nbiterrs = 0;
Nblkerrs = 0;
for i = 1 : Nblocks
    msg = zeros(1, k);   % all-zero message
    % Encoding
    cword = zeros(1, n); % all-zero codeword
    s = 1 - 2 * cword; % BPSK bit to symbol conversion
    r = s + sigma * randn(1, n); % AWGN channel
    
    % Soft-Decision, Iterative Message-Passing Layered Decodingp
    L = r; % total belief (LLR)
    itr = 0; % iterator number
    R = zeros(Slen, z); % storage for row processing
    while itr < MaxItrs
        Ri = 0;
        for lyr = 1 : mb
            ti = 0; % number of non -1 in row = lyr
            for col = 1 : nb
                if B(lyr, col) ~= -1
                    ti = ti + 1;
                    Ri = Ri + 1;
                    % Subtraction and row alignment
                    L((col - 1) * z + 1 : col * z) = L((col - 1) * z + 1 : col * z) - R(Ri, :);
                    % Row alignment and store in treg
                    treg(ti, :) = mul_sh(L((col - 1) * z + 1 : col * z), B(lyr, col));
                end
            end
            % minsum on treg: ti x z
            for i1 = 1 : z % treg(1 : ti, i)
                [min1, pos] = min(abs(treg(:, i1))); % first minimum
                min2 = min(abs(treg([1 : pos - 1 pos + 1 : ti], i1)));
                S = sign(treg(1 : ti, i1));
                parity = prod(S);
                treg(1 : ti, i1) = min1; % absolute value for all
                treg(pos, i1) = min2; % absolute value for min1 position
                treg(1 : ti, i1) = parity * S .* treg(1 : ti, i1); % assign signs
            end
            % column alignment, addition and sotre in R
            Ri = Ri - ti; % reset the storage counter
            ti = 0;
            for col = find(B(lyr,:) ~= -1)
                Ri = Ri + 1;
                ti = ti + 1;
                R(Ri, :) = mul_sh(treg(ti, :), z - B(lyr, col));
                % Addition
                L((col - 1) * z + 1 : col * z) = L((col - 1) * z + 1 : col * z) + R(Ri, :);
            end
        end
        msg_cap = L(1 : k) < 0; % decision
        itr = itr + 1;
    end
    
    % Counting Errors
    Nerrs = sum(msg ~= msg_cap);
    if Nerrs > 0
        Nbiterrs = Nbiterrs + Nerrs;
        Nblkerrs = Nblkerrs + 1;
    end
end

BER_sim = Nbiterrs / k / Nblocks;
FER_sim = Nblkerrs / Nblocks;

disp([EbNodB FER_sim BER_sim Nblkerrs Nbiterrs Nblocks])