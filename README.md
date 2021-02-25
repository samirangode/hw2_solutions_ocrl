# HW2: Diving into LQR
## Due date: Tuesday, March 9, 2021

In this homework we'll be delving into LQR and some QP-based MPC. Here's an overview of the problems.
1. Stabilize the quadruped on one leg using LQR.
2. Use TVLQR to stabilize a parallel parking maneuver
3. Use LQR and QP-based MPC to land a rocket

## Running the Autograder
The autograder setup has changed slightly for this assignment. We'll no longer use GitHub actions to run the autograder (since it didn't really work last time, maybe in future assignements...). You will run all your tests locally. To run the full test suite the same way we will when grading your submission follow these instructions:

1. Open a terminal in the root directory of your repo
2. Launch a Julia REPL
3. Enter the package manager using `]` and enter `activate .`
4. Launch the testing suite using `test hw2`

Each notebook now includes a `run_tests()` function at the end, that will run the test suite in your notebook. You can call that test at any point. It will just run the `q1.jl`,
`q2.jl` or `q3.jl` files in the `test/` directory.

## Submitting your homework
Make sure your repo lives under the Class Organization. This will be done automatically when you use the GitHub Classrooms link we send provide.

Follow [these instructions](https://github.com/Optimal-Control-16-745/JuliaIntro/blob/main/docs/Submission%20Instructions.md) for submitting your final time-stamped submission.

## Adding the Upstream Repo
We may release changes to the homework periodically if errors or bugs are found. Follow these instructions for linking your repo to the original template and pulling changes. It's always a good idea to branch your code before pulling from the upstream repo in case something goes wrong or the merge is particularly nasty. 