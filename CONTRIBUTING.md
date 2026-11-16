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

This inevitably means that, most of the time, you'll be merging unfinished work - which is fine, *as long as it doesn't break existing stuff*. Start with empty subroutines that are not yet called anywhere; as long as they compile, you're good. Then, start adding your new functionality. Want to see how it behaves when it's called from the actual model? Make sure it's **off by default** (best done by using config parameters). Feel like it's good to go? Add a test. Each of these steps is fine as a separate pull request.

We encourage you to use imperative branch naming. E.g., if you're planning to add a new climate model option, your branch will be called something like `add_new_climate_model_option`.

### Fetch/pull/merge every day

After the even more basic stuff that makes up the start of your day (I'd say getting out of bed, putting on clothes, and eating breakfast, but many programmers seem to view these steps as optional), your daily work routine should begin with:
- `fetch` the repo from the GitHub origin
- `pull`ing the changes to main (see https://stackoverflow.com/questions/292357/what-is-the-difference-between-git-pull-and-git-fetch to learn the difference),
- `merge` main into your branch

This means that you're never more than a day behind whatever work the other contributors have added to the codebase, which minimises the risk of merge conflicts.

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