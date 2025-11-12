# Astro 528 [Class Project](https://psuastro528.github.io/Fall2025/project/)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://PsuAstro528.github.io/project-template/stable)

GitHub Actions : [![Build Status](https://github.com/PsuAstro528/project-template/workflows/CI/badge.svg)](https://github.com/PsuAstro528/project-template/actions?query=workflow%3ACI+branch%3Amain)


## Class Project Schedule
- [Project proposal](proposal.md) (due Sept 12)
- Serial version of code (due Oct 8)
- Peer code review (due Oct 15)
- Parallel version of code (multi-core) (due Nov 5)
- Second parallel version of code (distributed-memory or GPU) (due Nov 21)
- Completed code, documentation, tests, packaging (optional) & reflection (due Dec 3)
- Class presentations (Dec 3 - Dec 12, [detailed schedule](https://github.com/PsuAstro528/PresentationsSchedule2025) )

# Code Description
Our code provides a semi-analytic model of mass transfer in a close white dwarf-neutron star binary (henceforth WD-NS, respectively) to determine an time series gravitational waveform. In principal, assuming typical neutron star mass $M_{\rm NS}=1.35\, M_\odot$, the user should be able to run the code with only four parameters: the initial white dwarf mass ($M_{\rm WD}$) and parameters specifying the binary orientation (inclination $i$ and initial phase angle $\varphi$) and luminosity distance ($D_L$ in Mpc). In the current version of the code, the second derivative of the quadrupole moment ($\ddot I^{ij}$) is calculated as a time series from initial $M_{\rm WD}$, and a separate method is given to calculate waveform $h^{ij}$ given $\ddot I^{ij}$ and $i,\,\varphi,\,D_L$. (In practice, the expensive and delicate part of the calculation is calculating $\ddot I^{ij}$ and it is fairly straightforward to then calculate gravitational waveforms for a variety of geometries afterward. To use our code for analysis, we would likely direct the user to output timeseries $\ddot{I}^{ij}$ as permanent files, and then from these template $\ddot{I}^{ij}$ to calculate $h^{ij}$ as temporary objects on which to perform analysis.)

Our code can broadly be divided into two parts. The first implements the mathematics behind mass transfer (also referred to as 'MT' in our documentation) to evolve a set of orbital elements $\vec{\theta}:=\[a, M_{\rm WD}, M_{\rm NS}, R_{\rm WD}, J, \theta\]$: orbital separation, white dwarf mass, neutron star mass, white dwarf radius, total angular momentum, and system phase, respectively. (\textit{Note: Throughout we will use $\vector{\theta}$ to refer to the set of parameters while we use $\theta$ to refer to the single angular parameter.} In principal, this may seem to add an unnecessary notational challenge, but in practice we refer to parameter set $\vec\theta$ often, while  angular $\theta$ only appears in calculating $\ddot I^{ij}$.) We largely follow [the lecture notes of Onno Pols](https://www.astro.ru.nl/~onnop/education/stev_utrecht_notes/). The mathematical underpinnings of this code can be found in the next section. `victor.jl` contains most of this part, including a suite of unit tests for our functions.

The second part of the code involves uniting the effects of mass transfer and '2.5 Post-Newtonian order' gravitational wave backreaction (also referred to as 'radiation reaction' or 'RR' in our documentation; see Chapter 12 of [Poisson & Will](https://www.cambridge.org/core/books/gravity/1F0CF38A59B1E51A63C7C3138268BE5D)) through autodifferentiation of the $\vec{\theta}$ parameters; integration of the system taking both MT and RR effects into account; calculation of timeseries $\ddot I^{ij}$; and conversion to $h^{ij}$. \textbf{This part of the code is significantly more complicated and time intensive than the previous section with serialized runtimes likely on the order of hours, and it is thus crucial to optimize. This is where peer review could be especially beneficial.} `tristan.jl` contains most of this part.

The key producable of our code is a time-series waveform given an initial mass WD orbiting a 1.35 $M_\odot$ NS. The gravitational waveform is calculated with `integrate(M_WD)` in `tristan.jl`. 

# Mathematical Basis
Please see [this note](https://www.overleaf.com/read/zzzprcdvzkzz#d5f227) which serves as the mathematical basis for our project and especially for `src/math.jl`. (I was unable to get these notes in `README.md` in a nice format.)

The units for everything are in solar masses, km, and seconds. Notably, the simulation's output file contains:
- t, time in seconds
- ddI_xx, second time derivative of quadrupole moment's xx component in solar mass * km^2
- ddI_xy, second time derivative of quadrupole moment's xy component in solar mass * km^2
- a, total orbital separation in km
- M_WD, mass of the white dwarf in solar mass
- M_NS, mass of the neutron star in solar mass
- R_WD, radius of the white dwarf in km
- J, total angular momentum in solar mass * km^2 / s
- $\theta$, orbital phase

# Installation and Using Package
To install, you must have Julia installed on your system first. Afterwards, type in the following three terminal commands:
```
git clone https://github.com/PsuAstro528/project-wd_ns_interaction.git
cd project-wd_ns_interaction
julia install.jl
```

To run the code, import Integrate... (type in rest of description).
The output is in ...
Parallelizing gives us...

<!--
The following section documents the mathematical treatment of the 2.5-PN RR and MT forces used to calculate a gravitational waveform in a close WD-NS binary. 
## Mass Transfer
Below are the kinetics used to evolve a circular orbit. The set of elements evolved in this framework $\vec\theta= \{a, M_{\rm WD}, M_{\rm NS}, R_{\rm WD}, J, \theta\}$, the orbital separation, white dwarf mass, neutron star mass, white dwarf radius, and angular velocity. Both radiation reaction and (isotropic) mass transfer tend to be circularizing effects, so we assume that the orbit stays circular. This makes tracking both $a$ and $J$ slightly redundant in theory, but in practice it is helpful since radiation reaction and mass transfer forces become comparable so the instantaneous relationship between velocity and separation is non-trivial.

The simple semi-analytic formalism that we use here is largely based on the lecture notes of [Pols](https://www.astro.ru.nl/~onnop/education/stev_utrecht_notes/). First, we establish the rate of mass transfer due to Roche Lobe overflow of the White Dwarf. We take the Eggleton approximation for the Roche Limit, $R_L$ where mass ratio $q:=\frac{M_{\rm WD}}{M_{\rm NS}}$:
$R_L=a\frac{0.49 q^{2/3}}{0.6q^{2/3}+\ln(1+q^{1/3})}$
Then, approximating orbital period, $P\equiv P(a,M_{\rm WD},M_{\rm NS})$, as keplerian, the mass loss rate of the white dwarf is cubic in the overflow of the white dwarf over the Roche limit:
$\dot M_{\rm WD}=-A\frac{M_{\rm WD}}{P}\left(\frac{R_{\rm WD}-R_L}{R_{\rm WD}}\right)^3$
where $A$ is a numerical constant determined from detailed fluid mechanics which we take to be $A=10$. The natural radius of a white dwarf is a function of its mass approximated as:
$R_{\rm WD,\,0}=10^4\unit{km}\left(\frac{M_{\rm WD}}{0.7\,M_\odot}\right)^{-1/3}\sqrt{1-\left(\frac{M_{\rm WD}}{M_{\rm Ch}}\right)^{4/3}}\left(\frac{\mu_e}{2}\right)^{-5/3}$
where $\mu_e$ is the mean molecular mass per electron which is $2$ for white dwarfs of all but the most exotic chemical compositions (e.g., ONeMg). We assume that the white dwarf radius evolves on the white dwarf dynamical timescale $\tau_{\rm dyn}=\sqrt{\frac{R_{\rm WD}^3}{GM_{\rm WD}}}$. Then:
$\dot R_{\rm WD}=\frac{R_{\rm WD,\,0}-R_{\rm WD}}{\tau_{\rm dyn}}$
We assume that mass transfer is conservative until the mass loss rate of the white dwarf exceeds the Eddington accretion rate of the neutron star:
$\dot M_{\rm edd,\,NS}=\frac{4\pi G m_{ p}}{\sigma_T\eta_{\rm acc}c}M_{\rm NS}$
where we assume that the accretion rate efficiency is $\eta_{\rm acc}=0.1$. This is to say that we take:
$\dot M_{\rm NS}=\min({-\dot M_{\rm WD},\,\dot M_{\rm edd,\,NS}})$
It is helpful to define a parameter to quantify mass loss rate: $\beta:=\frac{ \dot M_{\rm NS}}{-\dot M_{\rm WD}}=\min(1,\frac{ \dot M_{\rm NS,\,edd}}{-\dot M_{\rm WD}})$. Angular momentum loss due to non-conservative mass transfer, $\dot J_{\rm MT}$, is given by:
$\dot J_{\rm MT}=J\gamma(1-\beta)\frac{\dot M_{\rm WD}}{M_{\rm WD}+M_{\rm NS}}$
where $\gamma$ parameterizes the specific angular momentum of the material ejected from the binary. We assume isotropic re-emission (e.g., mass is transferred from the white dwarf to the neutron star and then a fraction $\beta$ of the mass is ejected), in which case $\gamma=\frac{M_{\rm WD}}{M_{\rm NS}}$. Finally, we derive the evolution of the binary separation due to mass loss. In general, changes in angular momentum, separation, mass, and eccentricity are related by:
$2\frac{\dot J}{J}=\frac{\dot a}{a}+2\frac{\dot M_{\rm WD}}{M_{\rm WD}}+2\frac{\dot M_{\rm NS}}{M_{\rm NS}}-\frac{\dot M_{\rm WD}+\dot M_{\rm NS}}{M_{\rm WD}+ M_{\rm NS}}-\frac{2 e\dot e}{1-e^2}\label{general}$
Then the change in binary separation due to mass transfer is:
$\dot a_{\rm MT}=-2a\frac{\dot M_{\rm WD}}{M_{\rm WD}}\bracs{1-\beta\frac{M_{\rm WD}}{M_{\rm NS}}-(1-\beta)\\left(\gamma+\frac{1}{2}\right)\frac{M_{\rm WD}}{M_{\rm WD}+M_{\rm NS}}}$

## Radiation Reaction
At 2.5-PN order, RR contributes to $\dot J,\,\dot a$ as functions of $\ddot I^{ij},\,\dddot I^{ij}$. In an effective one-body system, we have that $I^{ij}=\mu x^i x^j $ 
with $x^1:=x(\vec\theta)=a\cos\theta;\,x^2:=y(\vec\theta)=a\sin\theta$. Then, $\ddot I^{ij}\equiv \ddot I^{ij}\left(\vec\theta,\dot{\vec\theta},\ddot{\vec\theta}\right)$ and $\dddot I^{ij}\equiv \dddot I^{ij}\left(\vec\theta,\dot{\vec\theta},\ddot{\vec\theta},\dddot{\vec\theta}\right)$. This poses a computational challenge as formally we have that $\dot{\vec \theta}\equiv \dot{\vec \theta}\left(\vec\theta,\dot{\vec \theta},\ddot{\vec \theta},\dddot{\vec \theta}\right)$ when taking into account RR effects. We cicumvent this challenge by assuming that radiation reaction effects are small compared to the dominant central keplerian force, and we substitute the dependences on $\vec{\theta}$ and subsequent derivatives in $\ddot I^{ij},\,\dddot I^{ij}$ with $\vec{\theta}_{\rm MT}:=\{a_{\rm MT}, M_{\rm WD}, M_{\rm NS}, R_{\rm WD}, J_{\rm MT}, \theta\}$. We calculate higher derivatives using autodifferentiation.

Angular momentum loss due to gravitational wave emission is calculated as:
$$
\dot{J}^{jk}_{\mathrm{(rr)}} = -\frac{2G}{5c^5} \left( \ddot{I}^{jp} \dddot{I}^k_{\;\;p} - \ddot{I}^{kp} \dddot{I}^j_{\;\;p} \right)
$$

$$
J^j = \epsilon^j_{\;\;pq} J^{pq}
$$  

$$
\implies \dot{J}^z_{\mathrm{(rr)}} = \frac{-2G}{5c^5} \left( \ddot{I}^{xp} \dddot{I}^y_{\;\;p} - \ddot{I}^{yp} \dddot{I}^x_{\;\;p} \right)
$$
where we define axes such that the orbit lies in the $xy$-plane. Then, we have that the derivative in separation due to radiation reaction to 2.5 PN order is:
$\dot a_{\rm RR}=\frac{-4aG}{5Jc^5}(\ddot{I}^{xp}\dddot I^y_{\,\,\,p}-\ddot{I}^{yp}\dddot I^x_{\,\,\,p})$
So the total rate of angular momentum loss is:
$\dot J=\dot J_{\rm MT}+\dot J_{\rm RR}=J\gamma(1-\beta)\frac{\dot M_{\rm WD}}{M_{\rm WD}+M_{\rm NS}}-\frac{2G}{5c^5}(\ddot{I}^{xp}\dddot I^y_{\,\,\,p}-\ddot{I}^{yp}\dddot I^x_{\,\,\,p})$
and the total rate of change in orbital separation is:
$\dot a=\dot a_{\rm MT}+\dot a_{\rm RR}=-a\left\{2\frac{\dot M_{\rm WD}}{M_{\rm WD}}\left[1-\beta\frac{M_{\rm WD}}{M_{\rm NS}}-(1-\beta)\left(\gamma+\frac{1}{2}\right)\frac{M_{\rm WD}}{M_{\rm WD}+M_{\rm NS}}\right]+\frac{4G}{5Jc^5}(\ddot{I}^{xp}\dddot I^y_{\,\,\,p}-\ddot{I}^{yp}\dddot I^x_{\,\,\,p})\right\}$
## Quadrupole Formula and Gravitational Waveform
If we suppose that the event occurs at some inclination $i$ relative to us and define initial phase $\varphi$, we can find the transverse component of $\dddot I^{ij}$ as:
$\ddot I^{ij}_{\rm (T)}=\ddot I^{k\ell}e^{(i)}_{k}e^{(j)}_{\ell}$
for $e_1=(\cos i\cos\varphi,\sin i\sin\varphi,\cos i)$ and $e_2=(-\sin\varphi,\cos\varphi,0)$ are vectors perpendicular to the line of sight. In the Transverse-Traceless gauge, we have:
$\ddot I^{ij}_{\rm (TT)}=\ddot I^{ij}_{\rm (T)}-\delta^{ij}\ddot I^{k}_{\,\,\,k,\,{\rm (T)}}$
Finally, the observed gravitational waveform is proportional via the quadrupole formula:
$h^{ij}_{(TT)}=\frac{2G}{c^4D_L}\ddot I^{ij}_{\rm (TT)}$
for some luminosity distance $D_L$.-->