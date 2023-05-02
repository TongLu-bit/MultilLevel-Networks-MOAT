# MultilayerNetworks-MOAT  

This is a repository that contains codes for the **M**ultilayer netw**o**rk associ**a**tion me**t**hod (**MOAT**) to systematically investigate a vector-to-matrix association patterns. Specifically, it includes:

1. Matlab functions ("greedy_multilayer.m" and "greedy_innerlayer.m") for implementing MOAT to identify dense sub-networks that containly highly associated edges.

2. A live-script (which is matlab version of Rmarkdown) for MOAT anlaysis on a synthetic data ("toy_example.m")

3. To see how synthetic data is generated, please refer to "simulate_inputdata.m". Otherwise, one can download the input data "W_input.mat" and use it directly in "toy_example.m"



Here is an overview of MOAT. 
\
![Image](/pipeline.png)
