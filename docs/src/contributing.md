# Contributing

Additions to LiNCEGE.jl are welcome in the following forms:

1. **Issues** – Report bugs or suggest new features.
2. **Pull Requests** – Submit code changes or documentation updates.

Please ensure all your edits pass the tests located in `test/runtests.jl`.

### Pull Request Process

We use the feature branch workflow for contributions. Please follow these steps:

1. Fork the LiNCEGE.jl repository to your own GitHub account.

3. Create a new branch in your local fork for your specific changes.
   ```bash
   git checkout -b feature/your-feature-name
   ```

4. Implement your changes. Ensure everything works by running the tests, add tests as needed:
   `julia test/runtests.jl`

5. Commit your changes and push the branch to your fork:
   ```bash
   git add .
   git commit -m "Brief description of your changes"
   git push origin feature/your-feature-name
   ```

6. Go to the original Lincege repository on GitHub. Click **New Pull Request** and select your feature branch to submit.

### Running benchmarks

Install the benchmark environment:

```bash
julia --project=benchmark -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
```

To run the full benchmark suite, run the following code in the julia repl with command `julia --project=benchmark`.
```julia
using BenchmarkTools
include("benchmark/benchmarks.jl")
results = run(SUITE; verbose=true)
show(stdout, MIME"text/plain"(), results)
```

Running the full suite can take time, to measure just the square lattice benchmark, use:

```julia
results = run(SUITE["square"]; verbose=true)
show(stdout, MIME"text/plain"(), results)
```

### Building the Documentation for LiNCEGE.jl

Install the documentation environment:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
```

Build the docs locally with:

```bash
julia --project=docs docs/make.jl
```

Open `docs/build/index.html` in a browser to view the result.

### AI Policy
TODO
