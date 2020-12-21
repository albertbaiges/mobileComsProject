function throughput = throughput(bandwidth, SNR, gap)
    gap_lineal = 10^(gap/10); % Gap to lineal
    throughput = bandwidth*log2(1+SNR/gap_lineal); % compute throughput using the formula and return
end

