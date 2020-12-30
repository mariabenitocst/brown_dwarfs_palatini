## Brown Dwarfs in Palatini gravity

Code for calculating the evolution of the degeneracy parameter, i.e. $$\Psi$$, with respect to time in Palatini gravity (Starobinsky model).

The evolution of the degeneracy parameter is given by equation (XX) in Benito & Wojnar 2021. The integration of this equation is done in C++ and the plotting of the corresponding results is done in python3.


### Setup

The following steps have to be taken in order to use the python code in this repository:
- Download/clone this repository in order to have a local copy
- Use virtual environments, e.g. `venv`, in order to create a virtual python3 environment for this repository-
- Once the virtual environment is active use `pip3 install -r requirements.txt` in order to install all necessary packages.

C++ requirements:
- [gsl 2.5](https://www.gnu.org/software/gsl/)


### How to run
To calculate the time evolution of the degeneracy parameter, type the following two commands in terminal once you are inside the folder `MG/src`
```
$ make
$ ./evolution_degeneracy
```
You will be asked to introduce the value of the parameters $$\alpha$$, $$\gamma$$ and $$\delta$$.
After this, the code will generate the file `evolution_degeneracy_alpha=XXX.dat` that contains the resulted evolution of $$\Psi$$. This can be plotted using the jupyter-notebook named `evolution_degeneracy.ipynb`.
