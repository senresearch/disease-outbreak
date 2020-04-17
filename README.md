
# DiseaseOutbreak

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://chelseatrotter.github.io/DiseaseOutbreak.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://chelseatrotter.github.io/DiseaseOutbreak.jl/dev)
[![Build Status](https://travis-ci.com/chelseatrotter/DiseaseOutbreak.jl.svg?branch=master)](https://travis-ci.com/chelseatrotter/DiseaseOutbreak.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/chelseatrotter/DiseaseOutbreak.jl?svg=true)](https://ci.appveyor.com/project/chelseatrotter/DiseaseOutbreak-jl)
[![Codecov](https://codecov.io/gh/chelseatrotter/DiseaseOutbreak.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/chelseatrotter/DiseaseOutbreak.jl)
[![Coveralls](https://coveralls.io/repos/github/chelseatrotter/DiseaseOutbreak.jl/badge.svg?branch=master)](https://coveralls.io/github/chelseatrotter/DiseaseOutbreak.jl?branch=master)
[![Build Status](https://api.cirrus-ci.com/github/chelseatrotter/DiseaseOutbreak.jl.svg)](https://cirrus-ci.com/github/chelseatrotter/DiseaseOutbreak.jl) -->

# Modeling disease outbreaks

We use a generative discrete Markov chain to model the progress of an
infectious disease in a population.  It assumes that everyone within a
certain radius of an infected person is susceptible, and there is a
certain probability of being infected.  The input parameters are the
locations of all individuals, radius of infection, probability of
being infected on a given day, how long an infected person is
infectious.  You can seed the simulator with a the number and
positions of currently infected individuals.  The data directory
contains household information for Memphis.

TODO:

- Layer information about severity
- Add testing regime
- Add ability to estimate parameters from observed data
