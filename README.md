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
