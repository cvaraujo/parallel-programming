OmpCluster: Matrix Multiplication
================================================================================

In this lab you should parallelize the matrix multiplication algorithm on
multiple nodes using for example block partitioning (See [Wikipedia][wiki]).
This labs focuses on the speedup obtained when distributing the data and
computation on multiple nodes compared to single node execution. You should use
the OpenMP constructs compatible with OmpCluster.

[wiki]: https://en.wikipedia.org/wiki/Matrix_multiplication

### Input

The program expects two lines as input:

- The first line is the size of the square matrix;
- The second line is the seed to generate random numbers.

### Output

The program will output the result matrix and execution time in seconds.

Tasks & Rules
--------------------------------------------------------------------------------

You should do the following tasks:

- [ ] Understand the single node code in `src/matmul-serial.c`. Note that the
single node implementation is already using `parallel for`.
- [ ] Parallelize the code using OmpCluster in the file `src/matmul-parallel.c`.
- [ ] Run both versions **on Kahuna** and compare them. Did you get any speedup?

You must **not** change the serial implementation, only the parallel one.
Feel free to use C++ if you prefer to implement it using classes.  

Grading
--------------------------------------------------------------------------------

Your assignment will be evaluated in terms of:

- Correctness: your program returns the correct result;
- Performance: your program runs faster than the serial implementation.

In order to test your solution, you can use one of the inputs available inside
the `tests/` directory. Whenever you push your changes to GitHub, the Continuous
Integration (CI) system will compile, run and execute your program using those
tests. Your grade will **not** be computed from the CI runs because other
processes may be running in the server, thus interfering with runtime and
speedup.

Your grade will be computed using an automated routine restricted to the
instructors and TAs. This routine will be run after the assignment deadline,
using the latest commit push to the branch `master` before the deadline. Your
code will be ensured to run in an environment with no competition for resources.

**Note:** Both the CI and automatic grading routine expect your the output of
your program to be formatted correctly. For that reason, you should not add
`printf`s or any other function that writes to `stdout`, otherwise your
assignment will be considered incorrect.

**Note:** Tampering with the serial implementation or the tests is considered
cheating and will result in disqualification of the assignment.

Compiling & Running
--------------------------------------------------------------------------------

After you have accepted this assignment on the course's GitHub Classroom page,
clone it to your machine. First you should use CMake [CMake](https://cmake.org/)
to generate the build system and the OmpCluster fork of Clang to compile the
code. They are both available in the docker image provided for the project
`ompcluster/runtime:latest`. See below for more details on how to use it.

### Local computer

You can use OmpCluster on your own computer using docker and the following
command:

```bash
docker run -v /path/to/OMPC-MatMul/:/root/OMPC-MatMul -it ompcluster/runtime:latest /bin/bash
cd /root/plasma/
```
You can get more information on how to use Docker in the official
[Get Started][https://docs.docker.com/get-started/] guide).

You can also use Singularity (see [here][https://sylabs.io/guides/3.2/user-guide/]
for more informations):

```bash
singularity pull docker://ompcluster/runtime:latest
singularity shell ./runtime_latest.sif
cd /path/to/OMPC-MatMul/
```

Then, run the following commands:

```bash
# Where the build will live
mkdir build && cd build

# Use Clang as compiler
export CC=clang
# Generate the Makefiles
cmake -DCMAKE_BUILD_TYPE=Release ..
# Compile
make
```

Having done that, still inside the `build` directory, run `make` to compile
everything. Finally, from the root directory of the repo you can execute the
serial and parallel versions like so:

```bash
build/serial tests/1.in
mpirun -np 3 build/parallel tests/2.in
```

Remember that the OmpCluster runtime will only be used to schedule the tasks if
`mpirun -np <number_of_process>` is used to execute the application. Not using
`mpirun` means using the x86 device offloading that is mostly useful for
debugging OpenMP target code.

### Kahuna

A PBS script is provided with a functional configuration. Edit the `*.pbs` files
to eventually modify some configurations. Docker is not available on Kahuna so
the scripts are using [Singularity][singularity] to run the container.

You can test your parallelization on your own computer but you MUST run all
scalability tests using multiple nodes on the Kahuna cluster.

**DO NOT EXECUTE ANYTHING ON THE LOGIN NODE, ALWAYS USE PBS !!!**

If you have any doubts or run into problems, please contact the TAs. Happy
coding! :smile: :keyboard:

Contribute
--------------------------------------------------------------------------------

Found a typo? Something is missing or broken? Have ideas for improvement? The
instructor and the TAs would love to hear from you!

About
--------------------------------------------------------------------------------

This repository is one of the assignments handed out to the students in course
"MO644/MC970 - Introduction to Parallel Programming" offered by the Institute of
Computing at Unicamp.
