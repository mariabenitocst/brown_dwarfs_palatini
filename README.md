## Brown Dwarfs in Palatini gravity

Code for calculating the evolution of the degeneracy parameter, i.e. $\Psi$, with respect to time in Palatini gravity (Starobinsky model).

The evolution of the degeneracy parameter is given by equation (XX) in Benito & Wojnar 2021. The integration of this equation is done in C++ and the plotting of the corresponding results is done in python3.


### Requirements

C++ requirements:
- [gsl 2.5](https://www.gnu.org/software/gsl/)

python3 requirements:
- jupyter notebook
- normal stuff (matplotlib, numpy)

### How to run
To calculate the time evolution of the degeneracy parameter, type the following two commands in terminal once you are inside the folder `MG/src`
```
$ make
$ ./evolution_degeneracy
```
You will be asked to introduce the value of the parameters $$\alpha$$, $$\gamma$$ and $$\delta$$.
After this, the code will generate the file `evolution_degeneracy_alpha=XXX.dat` that contains the resulted evolution of $$\Psi$$. This can be plotted using the jupyter-notebook named `evolution_degeneracy.ipynb`.
