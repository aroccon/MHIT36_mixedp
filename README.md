# MHIT36

Multi-GPU code for intreface-resolved simulations of multiphase turbulence.
The code relies on direct numerical simulation of Navier-Stokes equations coupled with a phase-field method (ACDI) for interface description.
Tracking of Lagrangian particles (tracers) is also supported.
The code parallerelization relies on the cuDecomp library.
This a spin-off branch of the main code where the Poisson solver works in FP32.
This makes the code more energy-efficient and decrease the communication overhead.


~~~text
███    ███ ██   ██ ██ ████████ ██████   ██████       
████  ████ ██   ██ ██    ██         ██ ██              
██ ████ ██ ███████ ██    ██     █████  ███████   
██  ██  ██ ██   ██ ██    ██         ██ ██    ██     
██      ██ ██   ██ ██    ██    ██████   ██████        
~~~


If you use this code, please cite the following work: 
```bibtex
  @article{roccon2025,
  title   = {MHIT36: A Phase-Field Code for Gpu Simulations of Multiphase Homogeneous Isotropic Turbulence},
  author  = {Roccon, A. and Enzenberger, L. and Zaza, D. and Soldati, A.},
  journal = {Computer Physics Communications},
  year    = {2025},
  volume  = {314},
  issue   = {109804},
  doi     = {https://doi.org/10.1016/j.cpc.2025.109804}
}
```

![Test](val/render2.jpg)


## Reference performance and scaling
Performance (NS only)
* 128 x 128 x 128    |   1 x A100@Leonardo  |   1 ms/timestep
* 256 x 256 x 256    |   1 x A100@Leonardo  |   8 ms/timestep
* 512 x 512 x 512    |   1 x A100@Leonardo  |  65 ms/timestep 
* 128 x 128 x 128    |   4 x A100@Leonardo  |   1 ms/timestep
* 256 x 256 x 256    |   4 x A100@Leonardo  |   3 ms/timestep
* 512 x 512 x 512    |   4 x A100@Leonardo  |  18 ms/timestep 
* 512 x 512 x 512    |   4 x H100@MN5-ACC   |  14 ms/timestep 
* 1024 x 1024 x 1024 |   4 x A100@Leonardo  | 150 ms/timestep 
* 2048 x 2048 x 2048 |  64 x A100@Leonardo  | 330 ms/timestep
* 4096 x 4096 x 4096 | 256 x A100@Leonardo  | 780 ms/timestep

