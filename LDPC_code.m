load base_matrices/NR_1_2_20.txt
B = NR_1_2_20;
size(B) % 46 * 68
z = 20;
msg = randi([0 1], 1, 22 * z);
c = nrldpc_encode(B, z, msg);
check_cword(B, z, c)