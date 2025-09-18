# heterogeneous-Multi-Agent-Systen-
Data and control synthesis for the heterogeneous MAS

The agents' dynamics data and graph structure are stored in [AgentsData.m](https://github.com/TheBigoranger/heterogeneous-Multi-Agent-Systen-/blob/main/AgentsData.m)

The synthesis of the controller is executed in [hete_synthesis.m](https://github.com/TheBigoranger/heterogeneous-Multi-Agent-Systen-/blob/main/hete_synthesis.m). Please note that in order to run the script,  the following libraries are required:
- [YALMIP](https://yalmip.github.io/)
- [Mosek](https://www.mosek.com/)(optional)

The simulation is run separately in [hete_sim.m](https://github.com/TheBigoranger/heterogeneous-Multi-Agent-Systen-/blob/main/hete_sim.m) from the control synthesis. Therefore, if anyone wishes to change, please make sure to uncomment: 
```
hete_synthesis
``` 
and comment
```
%AgentsData
%load("controlGain.mat")
```
