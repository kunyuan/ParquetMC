## Benchmark for $\Gamma_4$

- **Test Case for the commit**: 7312ff669a415ae7b9544903913b35564a78ab63
    - Parameter:
    - ```
       #Order, Beta, rs,  Mass2,  Lambda, charge2,  Dim, Spin  TotalStep(1e6)
          4,    40,  1.0,  3.683,    3.0,    1.0,     3,    2,      401

        #TauGrid, MomGrid, AngGrid, MaxExtMom(*kF) 
        128,     4,       64,       3.0

        #Print, Save, ReWeight, Message, Collection /time interval (in Sec.)
        10,   10,     30,       10,      10

        #Reweight
        # 0    1    2,   3    4    5     6    7    8    9 
        10.0, 0.8, 0.4, 0.1, 0.4, 0.4, 1.0, 1.0, 1.0, 1.0  #gamma
      ``` 

    - Output from merge.py:
    - ```
        Order 1 
        Sum     Q/kF,    Data,    Error
        As:    0.00,   0.340001,   0.000039
        Aa:    0.00,  -0.118377,   0.000031
        Order 2 
        Sum     Q/kF,    Data,    Error
        As:    0.00,   0.357693,   0.000084
        Aa:    0.00,  -0.124859,   0.000047
        Order 3 
        Sum     Q/kF,    Data,    Error
        As:    0.00,   0.360479,   0.000125
        Aa:    0.00,  -0.126796,   0.000078
        Order 4 
        Sum     Q/kF,    Data,    Error
        As:    0.00,   0.360338,   0.000163
        Aa:    0.00,  -0.127452,   0.000085
      ```