
# DiseaseOutbreak

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://chelseatrotter.github.io/DiseaseOutbreak.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://chelseatrotter.github.io/DiseaseOutbreak.jl/dev)
[![Build Status](https://travis-ci.com/chelseatrotter/DiseaseOutbreak.jl.svg?branch=master)](https://travis-ci.com/chelseatrotter/DiseaseOutbreak.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/chelseatrotter/DiseaseOutbreak.jl?svg=true)](https://ci.appveyor.com/project/chelseatrotter/DiseaseOutbreak-jl)
[![Codecov](https://codecov.io/gh/chelseatrotter/DiseaseOutbreak.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/chelseatrotter/DiseaseOutbreak.jl)
[![Coveralls](https://coveralls.io/repos/github/chelseatrotter/DiseaseOutbreak.jl/badge.svg?branch=master)](https://coveralls.io/github/chelseatrotter/DiseaseOutbreak.jl?branch=master)
[![Build Status](https://api.cirrus-ci.com/github/chelseatrotter/DiseaseOutbreak.jl.svg)](https://cirrus-ci.com/github/chelseatrotter/DiseaseOutbreak.jl) -->

Tools for modeling disease outbreaks

This Julia package contains functions for simulating the course of an infectious disease outbreak and for fitting observed case count data to estimate parameters of SEIR models.  This is very much a work in progress.  If you plan on using or contributing to this work, please contact the authors or post an issue.

### Simulating disease outbreaks

We use a generative discrete Markov chain to model the progress of an
infectious disease in a population.  It assumes that everyone within a
certain radius of an infected person is susceptible, and there is a
certain probability of being infected.  The input parameters are the
locations of all individuals, radius of infection, probability of
being infected on a given day, how long an infected person is
infectious.  You can seed the simulator with a the number and
positions of currently infected individuals.  The data directory
contains household information for Memphis.

### SEIR model and variants

We use a discrete difference equations to model the course of an infectious disease outbreak.  We have implemented the SEIIIR model of [Hill](https://github.com/alsnhll/SEIR_COVID19) and the SIR-X model of [Maier and Brockmann](https://science.sciencemag.org/content/early/2020/04/07/science.abb4557.abstract.).  Our model fitting strategies differ from the above work and we are aiming for a unified interface to a variety of SEIR model variants.
