from io import StringIO
from scipy.interpolate import bisplev, bisplrep
from numpy import loadtxt
from cosmoslik import SlikPlugin

class bbn_consistency(SlikPlugin):
    
    def __init__(self):
        super(bbn_consistency,self).__init__()
        self.rep = bisplrep(*loadtxt(StringIO(dat)).T)

    def __call__(self,
                 ombh2,
                 Neff=3.046,
                 **kwargs):
        return bisplev(ombh2,Neff-3.046,self.rep)


dat = \
"""
#BBN prediction of the primordial Helium abundance $Y_p$ as 
#function of the baryon density $\omega_b h^2$ and number of 
#extra radiation degrees of freedom $\Delta N$.
#Calculated with PArthENoPE v1.00 [arXiv:0705.0290] for a 
#neutron lifetime of 885.7 s, using a network of 26 nuclides 
#and 100 reactions.

#ombh2  DeltaN      Yp

 0.005   -3.0    0.178540
 0.007   -3.0    0.183043
 0.009   -3.0    0.186014
 0.011   -3.0    0.188205
 0.013   -3.0    0.189936
 0.015   -3.0    0.191381
 0.016   -3.0    0.192031
 0.017   -3.0    0.192615
 0.018   -3.0    0.193171
 0.019   -3.0    0.193709
 0.020   -3.0    0.194185
 0.021   -3.0    0.194662
 0.022   -3.0    0.195089
 0.023   -3.0    0.195527
 0.024   -3.0    0.195921
 0.025   -3.0    0.196300
 0.026   -3.0    0.196655
 0.027   -3.0    0.197018
 0.028   -3.0    0.197358
 0.029   -3.0    0.197664
 0.030   -3.0    0.197991
 0.032   -3.0    0.198568
 0.034   -3.0    0.199121
 0.036   -3.0    0.199655
 0.038   -3.0    0.200148
 0.040   -3.0    0.200603
  
 0.005   -2.0    0.198733
 0.007   -2.0    0.203575
 0.009   -2.0    0.206688
 0.011   -2.0    0.208970
 0.013   -2.0    0.210757
 0.015   -2.0    0.212239
 0.016   -2.0    0.212910
 0.017   -2.0    0.213515
 0.018   -2.0    0.214066
 0.019   -2.0    0.214609
 0.020   -2.0    0.215112
 0.021   -2.0    0.215569
 0.022   -2.0    0.216027
 0.023   -2.0    0.216451
 0.024   -2.0    0.216840
 0.025   -2.0    0.217226
 0.026   -2.0    0.217595
 0.027   -2.0    0.217962
 0.028   -2.0    0.218300
 0.029   -2.0    0.218626
 0.030   -2.0    0.218929
 0.032   -2.0    0.219538
 0.034   -2.0    0.220096
 0.036   -2.0    0.220601
 0.038   -2.0    0.221106
 0.040   -2.0    0.221611
  
 0.005   -1.0    0.215291
 0.007   -1.0    0.220412
 0.009   -1.0    0.223655
 0.011   -1.0    0.226034
 0.013   -1.0    0.227863
 0.015   -1.0    0.229362
 0.016   -1.0    0.230012
 0.017   -1.0    0.230626
 0.018   -1.0    0.231199
 0.019   -1.0    0.231733
 0.020   -1.0    0.232237
 0.021   -1.0    0.232711
 0.022   -1.0    0.233179
 0.023   -1.0    0.233586
 0.024   -1.0    0.234004
 0.025   -1.0    0.234387
 0.026   -1.0    0.234757
 0.027   -1.0    0.235100
 0.028   -1.0    0.235451
 0.029   -1.0    0.235780
 0.030   -1.0    0.236089
 0.032   -1.0    0.236674
 0.034   -1.0    0.237229
 0.036   -1.0    0.237761
 0.038   -1.0    0.238252
 0.040   -1.0    0.238702
  
 0.005    0.0    0.229294
 0.007    0.0    0.234696
 0.009    0.0    0.238053
 0.011    0.0    0.240463
 0.013    0.0    0.242327
 0.015    0.0    0.243853
 0.016    0.0    0.244504
 0.017    0.0    0.245122
 0.018    0.0    0.245711
 0.019    0.0    0.246237
 0.020    0.0    0.246738
 0.021    0.0    0.247234
 0.022    0.0    0.247681
 0.023    0.0    0.248105
 0.024    0.0    0.248510
 0.025    0.0    0.248892
 0.026    0.0    0.249249
 0.027    0.0    0.249620
 0.028    0.0    0.249941
 0.029    0.0    0.250264
 0.030    0.0    0.250596
 0.032    0.0    0.251168
 0.034    0.0    0.251720
 0.036    0.0    0.252236
 0.038    0.0    0.252725
 0.040    0.0    0.253202
  
 0.005    1.0    0.241373
 0.007    1.0    0.247034
 0.009    1.0    0.250510
 0.011    1.0    0.252954
 0.013    1.0    0.254858
 0.015    1.0    0.256383
 0.016    1.0    0.257056
 0.017    1.0    0.257687
 0.018    1.0    0.258266
 0.019    1.0    0.258809
 0.020    1.0    0.259297
 0.021    1.0    0.259785
 0.022    1.0    0.260220
 0.023    1.0    0.260662
 0.024    1.0    0.261068
 0.025    1.0    0.261433
 0.026    1.0    0.261824
 0.027    1.0    0.262166
 0.028    1.0    0.262508
 0.029    1.0    0.262813
 0.030    1.0    0.263141
 0.032    1.0    0.263724
 0.034    1.0    0.264276
 0.036    1.0    0.264785
 0.038    1.0    0.265266
 0.040    1.0    0.265724
  
 0.005    2.0    0.252003
 0.007    2.0    0.257924
 0.009    2.0    0.261466
 0.011    2.0    0.263970
 0.013    2.0    0.265884
 0.015    2.0    0.267451
 0.016    2.0    0.268125
 0.017    2.0    0.268745
 0.018    2.0    0.269324
 0.019    2.0    0.269865
 0.020    2.0    0.270371
 0.021    2.0    0.270850
 0.022    2.0    0.271290
 0.023    2.0    0.271701
 0.024    2.0    0.272104
 0.025    2.0    0.272489
 0.026    2.0    0.272865
 0.027    2.0    0.273201
 0.028    2.0    0.273550
 0.029    2.0    0.273870
 0.030    2.0    0.274168
 0.032    2.0    0.274749
 0.034    2.0    0.275312
 0.036    2.0    0.275813
 0.038    2.0    0.276297
 0.040    2.0    0.276730
  
 0.005    3.0    0.261458
 0.007    3.0    0.267617
 0.009    3.0    0.271270
 0.011    3.0    0.273816
 0.013    3.0    0.275757
 0.015    3.0    0.277299
 0.016    3.0    0.277980
 0.017    3.0    0.278603
 0.018    3.0    0.279182
 0.019    3.0    0.279721
 0.020    3.0    0.280241
 0.021    3.0    0.280700
 0.022    3.0    0.281168
 0.023    3.0    0.281587
 0.024    3.0    0.281984
 0.025    3.0    0.282350
 0.026    3.0    0.282717
 0.027    3.0    0.283062
 0.028    3.0    0.283395
 0.029    3.0    0.283725
 0.030    3.0    0.284019
 0.032    3.0    0.284615
 0.034    3.0    0.285150
 0.036    3.0    0.285642
 0.038    3.0    0.286136
 0.040    3.0    0.286563
  
 0.005    4.0    0.269959
 0.007    4.0    0.276363
 0.009    4.0    0.280091
 0.011    4.0    0.282694
 0.013    4.0    0.284650
 0.015    4.0    0.286203
 0.016    4.0    0.286884
 0.017    4.0    0.287508
 0.018    4.0    0.288109
 0.019    4.0    0.288626
 0.020    4.0    0.289133
 0.021    4.0    0.289625
 0.022    4.0    0.290069
 0.023    4.0    0.290472
 0.024    4.0    0.290887
 0.025    4.0    0.291263
 0.026    4.0    0.291630
 0.027    4.0    0.291957
 0.028    4.0    0.292286
 0.029    4.0    0.292606
 0.030    4.0    0.292910
 0.032    4.0    0.293507
 0.034    4.0    0.294049
 0.036    4.0    0.294538
 0.038    4.0    0.294988
 0.040    4.0    0.295446
  
 0.005    5.0    0.277687
 0.007    5.0    0.284287
 0.009    5.0    0.288140
 0.011    5.0    0.290752
 0.013    5.0    0.292749
 0.015    5.0    0.294309
 0.016    5.0    0.295010
 0.017    5.0    0.295633
 0.018    5.0    0.296195
 0.019    5.0    0.296753
 0.020    5.0    0.297238
 0.021    5.0    0.297709
 0.022    5.0    0.298155
 0.023    5.0    0.298575
 0.024    5.0    0.298971
 0.025    5.0    0.299366
 0.026    5.0    0.299724
 0.027    5.0    0.300081
 0.028    5.0    0.300381
 0.029    5.0    0.300697
 0.030    5.0    0.301001
 0.032    5.0    0.301571
 0.034    5.0    0.302118
 0.036    5.0    0.302617
 0.038    5.0    0.303076
 0.040    5.0    0.303521
  
 0.005    6.0    0.284707
 0.007    6.0    0.291575
 0.009    6.0    0.295503
 0.011    6.0    0.298173
 0.013    6.0    0.300152
 0.015    6.0    0.301754
 0.016    6.0    0.302424
 0.017    6.0    0.303069
 0.018    6.0    0.303633
 0.019    6.0    0.304189
 0.020    6.0    0.304698
 0.021    6.0    0.305167
 0.022    6.0    0.305586
 0.023    6.0    0.306024
 0.024    6.0    0.306404
 0.025    6.0    0.306800
 0.026    6.0    0.307155
 0.027    6.0    0.307482
 0.028    6.0    0.307807
 0.029    6.0    0.308119
 0.030    6.0    0.308420
 0.032    6.0    0.308988
 0.034    6.0    0.309535
 0.036    6.0    0.310032
 0.038    6.0    0.310467
 0.040    6.0    0.310901
  
 0.005    7.0    0.291188
 0.007    7.0    0.298279
 0.009    7.0    0.302304
 0.011    7.0    0.305005
 0.013    7.0    0.307002
 0.015    7.0    0.308598
 0.016    7.0    0.309287
 0.017    7.0    0.309917
 0.018    7.0    0.310498
 0.019    7.0    0.311034
 0.020    7.0    0.311536
 0.021    7.0    0.312009
 0.022    7.0    0.312469
 0.023    7.0    0.312866
 0.024    7.0    0.313276
 0.025    7.0    0.313654
 0.026    7.0    0.314017
 0.027    7.0    0.314360
 0.028    7.0    0.314658
 0.029    7.0    0.314971
 0.030    7.0    0.315268
 0.032    7.0    0.315853
 0.034    7.0    0.316374
 0.036    7.0    0.316840
 0.038    7.0    0.317298
 0.040    7.0    0.317730
"""
