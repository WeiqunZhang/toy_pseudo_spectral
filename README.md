# Toy pseudo-spectral solver

In order to build and run the code with `docker`:

- Clone this repository into a directory which also contains picsar and amrex (`development` branch)

- Go into the folder `toy_pseudo_spectral/docker_build` and type
```
docker build -t tps .
```

- Go back in the folder that contains `picsar`, `amrex` and `toy_pseudo_spectral` and type
```
docker run -v $PWD:/test tps
```