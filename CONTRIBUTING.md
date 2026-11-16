# Contributing to UPSY-models

## Who can contribute?

You can! There are no rules for who can contribute - but there are rules for *how* to contribute, which are explained here.

## Before you start

- Have a look at the Wiki (HOW TO LINK TO THE ACTUAL WIKI?)
- Get your development environment set up (LINK TO RELEVANT WIKI PAGE)
- Verify that you can check out, compile, and run the code (NEED TO CREATE A SMALL DEMO SIMULATION)
- Ask if anyone is already working on the thing you're planning to work on (LINK TO GITHUB DISCUSSIONS)

## Starting

### Create a branch from `main` for your new work.

*Not sure how to do this? Take a look here: https://docs.github.com/en/desktop/overview/getting-started-with-github-desktop*

`Branches` should be short-lived, ideally about a week or so. Most of the time, this is nowhere near anough to finish your work; that's not the point. This way, you're never away from main for very long, and the risk of merge conflicts between your work and that of other contributors is minimised. Also, it encourages you to chop up your big project into small, manageable chunks, which are much easier to test and review.

We encourage you to use imperative branch naming. E.g., if you're planning to add a new climate model option, your branch will be called something like `add_new_climate_model_option`.


### Start your day by fetching and pulling
(SHOULD THIS BE IN HERE? OR JUST REFER TO CI/CD BEST PRACTICES?)

You'll often find yourself switching back to `main` to look at the code, re-run test cases, run production simulations, etc. Since main is continuously developed, stay up-to-date by fetching from origin and pulling the resulting changes into your local main.

Create a test (forces you to define what your new code should do)

Compile as often as possible (fail early)

## Merging


### Check that you're happy

Does it compile? Does it pass whatever new tests you created?

### Create a pull request

(TEMPLATE?)
Describe what you did and why
Small PR: single change (one single new feature, refactor of one module)
New feature: off by default, includes test
In case of changing output: clearly explain why, update references

### Check the test suite

If all tests pass: great! Ask for reviews.
If not: oh dear! REFER TO WIKI PAGE ABOUT TEST SUITE AND HOW TO INTERPRET IT

### Ask for reviews

Oh dear, reviewer #2 has posted a new comment!

Don't worry, code reviews are not like peer reviews on scientific literature. Yes, it's partly to check whether what you did makes sense - after all, the fact that your code compiles and runs doesn't mean it's good code. But it's also a great way for the rest of us to stay up to date on what's being added to the codebase.

## From errors to bug reports

(TEMPLATE?)

Even automated testing and code reviews are no guaranteed protection against plain old human error, and bugs *do* sneak into our codebase. If you run into one and don't feel up to the task of fixng it yourself, you can submit a bug report for one of our core developers to look into.

The key thing here is *reproducability*. This includes the commit hash of your source code, the versions of the external packages (all of which are printed to the terminal at the start of every simulation), and whatever input parameters you used (the config file and probably some NetCDF input files).

If possible, try to reduce this to the smallest possible set-up. If your simulation read 10 Gigabytes of input and then crashed after an hour of running, reproducing the error is not really feasible for us. Ideally, we want something that takes less than a minute to crash.

## Getting help

LINK TO GITHUB DISCUSSIONS

## Credit

Rules for acknowledgements/coauthorship