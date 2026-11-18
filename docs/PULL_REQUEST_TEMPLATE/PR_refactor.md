# Refactor

<!--
Use this template for PRs that change how the code is written, but not what it does:
restructuring a module, renaming things, moving code around, cleaning up, etc.
Text between these arrows won't show up in the PR; feel free to leave it in.
-->

## What did you refactor, and why?

<!-- Which module(s)/routine(s), and what was wrong with the old structure?
     E.g. "Split remapping_basic into submodules, because the file was 5000 lines long." -->

## Does this change model output?

<!-- A refactor should normally NOT change any results. Tick one: -->

- [ ] No: all checksums and test references are unchanged.
- [ ] Yes: I've explained why below, and I've updated the references.

<!-- If yes: explain which output changed, by how much, and why that's expected
     (e.g. a different but equivalent order of summation causing round-off differences). -->

## Checklist

- [ ] This PR refactors one thing (e.g. a single module), and doesn't add any new functionality.
- [ ] The code compiles, and the test suite passes.
- [ ] New and changed code follows the [style guide](https://github.com/UPSY-group/UPSY-models/blob/main/docs/STYLE_GUIDE.md).
- [ ] I've removed any code that has become unused.
- [ ] If I renamed or removed any config parameters, I've updated the config files in the test suite and listed the changes here.
