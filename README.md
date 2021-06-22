# dnsbox
dnsbox is a fortran code for the direct numerical simulation (DNS) 
of the sinusoidally forced Navier-Stokes equations 
(Kolmogorov flow) in a triply periodic domain. 

The version provided here is identical to the one used for the results 
of 

[[YHB2021]](https://link.aps.org/doi/10.1103/PhysRevLett.126.244502)
G. Yalniz, B. Hof, N. B. Budanur, 
*Coarse Graining the State Space of a Turbulent Flow Using Periodic Orbits*. 
Physical Review Letters **126**, 244502 (2021), 
[arXiv:2007.02584](https://arxiv.org/abs/2007.02584).

## Compilation

To compile in a Linux environment, you need the development libraries of:
 - `openmpi`
 - `fftw3`
 - `BLAS`
 - `LAPACK`

and `GFortran` (version 9 or higher).
On Ubuntu,
the packages "`build-essential gfortran libfftw3-dev libopenmpi-dev libblas-dev liblapack-dev`"
should be sufficient.

To compile the simulator, do
```
make
```

## Running

Simulations are started from pair of files:
`state.000000` containing the state data and `parameters.in` containing
the physical, output and debugging parameters.

Such a sample initial condition that leads to turbulence with a lifetime
longer than 10 000 is in `test/`.
To run it, you can create a folder, say, `rundir`,
```
mkdir rundir
```
copy the simulator binary (`dns.x`), the initial condition (`state.000000`) and the parameter input file
(`parameters.in`) there,
```
cp dns.x rundir/
cp test/* rundir/
```
go to `rundir`,
```
cd rundir
```
and run the simulator on `N` cores, detaching from the terminal, and redirecting
`stdout` (`1`) and `stderr` (`2`) to the file `log`,
```
nohup mpirun -np N dns.x > log 2>&1 &
```
Take note of the output process ID, say `123456`, you can use it to kill
the simulation,
```
kill 123456
```
unless it stops by itself due to laminarization, runtime limits (see `parameters.in`) 
or errors.

One can also start simulations from random initial conditions by setting `IC` to `0`.

## Utilities

As the run goes on, it will write state files (`state.123456`) and a file
containing observables and time-stepper data (`stat.gp`).
To visualize `stat.gp`, you can do
```
analysis ./ 0 -1
```

## Periodic orbit database

Along with the code, we provide a database of all periodic orbits that 
we discovered at the parameter values reported in 
[YHB2021](https://arxiv.org/abs/2007.02584) .
These can be found in the subfolders `periodic_orbits/*`. 
In these folders, `parameters.in` 
files are set such that running `dns.x` integrates the periodic orbit 
for one full period. 
Thus, the first line of `guess.in` corresponds to the 
converged period of the respective periodic orbit. 

## Citing

If you use `dnsbox` in your research, please cite

- [[YHB2021]](https://link.aps.org/doi/10.1103/PhysRevLett.126.244502)
G. Yalniz, B. Hof, N. B. Budanur, 
*Coarse Graining the State Space of a Turbulent Flow Using Periodic Orbits*. 
Physical Review Letters **126**, 244502 (2021), 
[arXiv:2007.02584](https://arxiv.org/abs/2007.02584).
