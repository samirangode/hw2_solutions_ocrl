# HW2: Diving into LQR
## Due date: Friday, February 25, 2022

In this homework we'll be delving into LQR and some QP-based MPC. Here's an overview of the problems.
1. Stabilize the quadruped on one leg using LQR.
2. Use TVLQR to stabilize a parallel parking maneuver
3. Use LQR and QP-based MPC to land a rocket


## Submission
Please follow these instructions. If you do not follow all of these instructions
you will have points deducted from your homework.
1. Fill out your name and AndrewID in `src/HW1.jl`
2. Commit and push all of your changes to GitHub
3. Create a release called "Submission" with a tag `v1.x` where `x` can be any number (start with `v1.0`, if you need to "patch" your submission change it to `v1.1`, then `v1.2`, i.e. follow [semantic versioning](https://semver.org/)). For more details, see [these instructions](https://github.com/Optimal-Control-16-745/JuliaIntro/blob/main/docs/Submission%20Instructions.md).


## Grading Policy
See [this document](https://github.com/Optimal-Control-16-745/JuliaIntro/blob/main/docs/Submission%20Instructions.md) for details on our grading policy. 
Do NOT attempt to hack the autograder. We will be looking at all of your 
solutions and any attempt to modify the autograder in any way will result in
heavy penalties. It's there to help make your lives easier, don't abuse it!

## Getting Started
All the homeworks are distributed as Jupyter notebooks. Follow these instructions to get everything set up.

1. Install Julia. Download v1.5.3 from the [Julia website](https://julialang.org/downloads/). Extract the binaries onto a stable location on your computer, and add the `/bin` folder to your system path.
2. Clone this repo, and open a terminal in the root directory
2. Start a Julia REPL in your terminal using `julia`. This should work as long as the binaries are on your system path.
3. Install the [IJulia](https://github.com/JuliaLang/IJulia.jl) using the Julia package manager. In the REPL, enter the package manager using `]`, then enter `add IJulia` to add it to your system path.
4. In the REPL (hit backspace to exit the package manager), enter `using IJulia`
5. Launch the notebook using `notebook()` or `jupyterlab()`

Check out Kevin's [video walkthrough](https://www.youtube.com/watch?v=I2SC1Mp3Hxs&feature=youtu.be) for HW0 for more details.


## Tips for Success
These homeworks can be pretty challenging, which is why we provide 2 weeks to complete them. Here are some tips for success:
1. Start early (Question 3 takes a lot of time!!!)
2. Commit and push changes often. Don't be afraid to use git branches. This is 
super helpful if you ever need to go back to a previous version because you broke something.
3. Get with a study group!
4. Keep your code clean. Avoid global variables. If you get odd behavior, consider restarting your kernel.
5. Come to office hours


## Questions / Issues
If you have questions as you work through the homework, please post to the 
`hw2` folder on Pizza. 

## Adding the Upstream Repo
We may release changes to the homework periodically if errors or bugs are found. Follow these instructions for linking your repo to the original template and pulling changes. It's always a good idea to branch your code before pulling from the upstream repo in case something goes wrong or the merge is particularly nasty. 
