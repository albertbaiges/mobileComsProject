function throughput = throughput(bandwidth, SNR, gap)
    gap_lineal = 10^(gap/10);
    throughput = bandwidth*log2(1+SNR/gap_lineal);
end

