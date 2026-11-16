# Contributing to UPSY-models

## Who can contribute?

You can! There are no rules for who can contribute - but there are rules for *how* to contribute, which are explained here.

## Before you start

- Have a look at the [UPSY Wiki](https://github.com/UPSY-group/UPSY-models/wiki)
- Set up your [development environment](https://github.com/UPSY-group/UPSY-models/wiki/Getting-started)

## Starting

### Inform the team what you're working on

We like to keep each other informed about what we're working on. We use the [GitHub Discussions](https://github.com/UPSY-group/UPSY-models/discussions) for this; if you're planning to add a new feature to the model code, start by creating a Discussion describing what you're going to build, and what you'll be using it for.

Note that you don't need to do this for every new branch/PR; just for the feature as a whole.

### Create a branch from main for your new work.

*Not sure how to do this? Take a look here: https://docs.github.com/en/desktop/overview/getting-started-with-github-desktop*

Branches should be short-lived, ideally a week at most. Of course, this is nowhere near enough to finish your project; that's not the point. This way, you're never away from main for very long, and the risk of merge conflicts between your work and that of other contributors is minimised. Also, it encourages you to chop up your big project into small, manageable chunks, which are much easier to test and review.

This inevitably means that, most of the time, you'll be merging unfinished work - which is fine, *as long as it doesn't break existing stuff*. Start with empty subroutines that are not yet called anywhere; as long as they compile, you're good. Then, start adding your new functionality. Want to see how it behaves when it's called from the actual model? Make sure it's **off by default** (best done by using config parameters). Feel like it's good to go? Add a test. Each of these steps is fine as a separate pull request.

We encourage you to use imperative branch naming. E.g., if you're planning to add a new climate model option, your branch will be called something like `add_new_climate_model_option`.

Lastly, and this really should go without saying: **don't touch someone else's branch!** ...at least, not without asking.

## Working

### Fetch/pull/merge every day

After the even more basic stuff that makes up the start of your day (I'd say getting out of bed, putting on clothes, and eating breakfast, but many programmers seem to view these steps as optional), your daily work routine should begin with:
- `fetch` the repo from the GitHub origin
- `pull` the changes to main (see https://stackoverflow.com/questions/292357/what-is-the-difference-between-git-pull-and-git-fetch to learn the difference),
- `merge` main into your branch

This means that you're never more than a day behind whatever work the other contributors have added to the codebase, which minimises the risk of merge conflicts.

### Compile often

The amount of time it costs, and the amount of frustration it causes, to trace an error in new code scales exponentially with the amount of code. The more often you compile, the smaller the amount of new code between checks, and the faster you find out if you messed up.

### Use existing debugging tools

While unfortunately, Fortran support in actual debuggers is flaky at best, UPSY offers a number of extremely useful tools to help you diagnose what's going on/going wrong in your code. Take a look at *LINK TO TO-BE-WRITTEN WIKI PAGE ABOUT UPSY DEBUGGING TOOLS* to see how these tools can make your life easier!

### Code in style

*LINK TO TO-BE-WRITTEN STYLE GUIDE*

## Merging

### Check that you're happy

The basic checklist to go over before creating a pull request is:
- does my code compile?
- does it pass whatever new tests I created?
- does it do what I intended when I created the branch?
- does it *not* do anything else?

Only when the answer to all these questions is "yes" should you create a PR.

### Create a draft pull request

*LINK TO TO-BE-WRITTEN PR TEMPLATE*
*Describe what you did and why*
*Small PR: single change (one single new feature, refactor of one module)*
*New feature: off by default, includes test*
*In case of changing output: clearly explain why, update references*

### Check the test suite

*If all tests pass: great! Un-draft your PR and ask for reviews.*
*If not: oh dear! REFER TO WIKI PAGE ABOUT TEST SUITE AND HOW TO INTERPRET IT*

### Ask for reviews

Oh dear, reviewer #2 has posted a new comment! Don't worry, code reviews are not like peer reviews on scientific literature. Yes, it's partly to check whether what you did makes sense - after all, the fact that your code compiles and runs doesn't mean it's good code. But it's also a great way for the rest of us to stay up to date on what's being added to the codebase.

## From errors to bug reports

### Ask for help!

Running into errors you don't understand is a big part of the life of a programmer. And while spending hours to solve a problem yourself is fine for the actual scientific problems you're paid to study, the simple programming problems that cause 99% of error messages have almost always already been solved by someone else at some point. So: ask! The [GitHub Discussions](https://github.com/UPSY-group/UPSY-models/discussions) are the place to go. Start with a quick search to see if someone else has already run into the problem. If not, create a new Discussion describing your problem. Chances are, you'll be sorted within a day.

### Create an issue

Not all errors are used-caused. Even automated testing and code reviews do not offer 100% protection, and bugs *do* sneak into our codebase. If you run into one and don't feel up to the task of fixing it yourself, you can [create an issue](https://github.com/UPSY-group/UPSY-models/issues) for one of our core developers to look into.

The key thing here is *reproducability*. This includes the commit hash of your source code, the versions of the external packages (all of which are printed to the terminal at the start of every simulation), and the config file. If you used input files that are not in the `UPSY-group/data` repo, check if the error persists when using those data files. If you really need external files, provide them along with your bug report using WeTransfer or a similar file-sharing service.

If possible, try to reduce the error using the smallest possible set-up. If your simulation read 10 Gigabytes of input and then crashes after an hour of running, reproducing the error is not really feasible for us. Ideally, we want something that takes less than a minute to crash, and uses only input files from the data repo. If you provide external files, make them as small as possible.

## Credit

*Rules for acknowledgements/coauthorship*